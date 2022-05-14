## Blinded sample size re-estimation simulator
library(survival)

## File path for enrollment table
path_name<-"/Users/jonathanluu/Dropbox/Winter Session/Simulation Code/enrollment_rate.csv"
recruitment_rate <- ceiling(read.table(path_name, header = TRUE))$enrollment # vector of # acrrued each month



##### Single simulation #####
## Variable descriptions
# recruitment_rate -> vector of how many patients have been accrued each month (cumulative)
# FU_time -> total follow up time of the study in months
# null_event_rate -> scale parameter for Weibull distribution - how often all participants 
#     will have an event in 1 years time under the null hypothesis
# lower_event_rate -> scale parameter for Weibull distribution - how often participants
#     from the treatment group will have an event under the alternative hypothesis
# upper_event_rate -> scale parameter for Weibull distribution - how often participants
#     from the control group will have an event under the alternative hypothesis
# LTFU_rate -> scale parameter for Weibull distribution - how often participants will have an event in 1 years time
# max_patient_FU_time -> if 0, there is no maximum time that patients are followed. if non-zero, maximum amount of time each patient can  be followed for 

## Result table values
# Status=1 -> Censored due to LTFU or end of study
# Status=2 -> Had an event
# Arm=1 -> Control arm
# Arm=2 -> Treatment arm

runSimulation <- function(recruitment_rate, FU_time, null_event_rate, lower_event_rate, upper_event_rate, LTFU_rate, max_patient_FU_time=-99){
  # Initialize variables and storage
  current_participant<-1
  sample_size=tail(recruitment_rate, n=1)
  recruit_time=length(recruitment_rate)

  results <- data.frame(Participant=1:sample_size, Time_to_event_null=0, Status_null=-1, Arm=-1,
                        Time_to_event=0, Status=0, Recruitment_month=0)
  
  ## Simulation code
  # iterate through each month
  for (i in 1:recruit_time){
    # obtain the number of people recruited this month
    recruited_this_month <- recruitment_rate[i]
    
    # iterate through each new patient
    results[current_participant:recruited_this_month,]<-t(apply(results[current_participant:recruited_this_month,], MARGIN=1, function(x){
      ## Null simulation
      # use a weibull RV to calculate null LTFU and event times
      event_time_null<-rweibull(1, shape=1, scale=null_event_rate)
      ltfu_time_null<-rweibull(1, shape=1, scale=LTFU_rate)
      
      # find the minimum of these three times and store result - if max_patient_FU_time is specified, also check this value
      all_times_null<-c(ltfu_time_null, FU_time-i, event_time_null, max_patient_FU_time)
      x[2]<-min(all_times_null[all_times_null>0])
      
      # set censoring status
      if (which.min(all_times_null[all_times_null>0])==3)
        x[3]<-2
      else
        x[3]<-1

    
      ## Alternative simulation
      # determine which arm the patient is randomized to
      x[4]<-sample(c(1,2),1)
      
      # assign event rate based on which arm
      if(x[4]==1)
        event_rate<-upper_event_rate
      else if (x[4]==2)
        event_rate<-lower_event_rate
      
      # use a weibull RV to calculate alternative LTFU and event times
      event_time<-rweibull(1, shape=1, scale=event_rate) 
      ltfu_time<-rweibull(1, shape=1, scale=LTFU_rate)
      
      # find the minimum of these three times and store result - if max_patient_FU_time is specified, also check this value
      all_times<-c(ltfu_time, FU_time-i, event_time, max_patient_FU_time)
      x[5]<-min(all_times[all_times>0])
      
      # set censoring status
      if (which.min(all_times[all_times>0]) == 3)
        x[6]<-2
      else
        x[6]<-1
      
      # store current month as recruitment month
      x[7]<-i
      
      return(x)
    } ))
    
    current_participant<-recruited_this_month+1
  }

  return(results)
}



## Find Weibull Scale parameters for event and LTFU rate 
# x = time (in months)
# survival_percentage - proportion of people who have not had an event at time x
findScale<- function(x, control_event_rate){
  return(-x/log(1-control_event_rate))
}

ner<-findScale(12, 0.36)
ler<-findScale(12, 0.32)
uer<-findScale(12, 0.40)
ltfur<-findScale(12, 0.1)



## Run the simulation num_rep times
set.seed(173)
num_reps<-1000

## Wrapper function to run num_rep simulations for power and type 1 error calculations
## Output variables
# null_reject - list of num_rep hypotheses values under the null -> value=1 if p<0.05, else 0 
# alt_reject - list of num_rep hypotheses values under the alternative -> value=1 if p<0.05, else 0

# num_event_null - list of num_rep event counts under the null
# num_event_ctrl - list of num_rep event counts under the alternative, assigned to control group
# num_event_trt - list of num_rep event counts under the alternative, assigned to treatment group

# mean_follow_up_null - list of num_rep mean follow up times under the null
# mean_follow_up_ctrl - list of num_rep mean follow up times under the alternative, assigned to control group
# mean_follow_up_trt - list of num_rep mean follow up times under the alternative, assigned to treatment group

runSimulations <- function(...){
  # store rejected and accepted hypotheses values
  null_reject<-alt_reject<-integer(num_reps) # Value=1 if p < 0.05, else 0
  
  # store number of events
  num_event_null<-num_event_ctrl<-num_event_trt<-integer(num_reps)

  # store mean follow up time
  mean_follow_up_null<-mean_follow_up_ctrl<-mean_follow_up_trt<-integer(num_reps)
  
  ## Scenario 1
  for (i in 1:num_reps){
    simulationResults<-runSimulation(...)
    
    # Log rank test
    null_p<-pchisq(survdiff(Surv(Time_to_event_null, Status_null)~Arm, data=simulationResults)$chisq,1,lower.tail=FALSE)
    alt_p<-pchisq(survdiff(Surv(Time_to_event, Status)~Arm, data=simulationResults)$chisq,1,lower.tail=FALSE)
    
    if (null_p < 0.05)
      null_reject[i]<-1
    else
      null_reject[i]<-0
    
    if (alt_p < 0.05)
      alt_reject[i]<-1
    else
      alt_reject[i]<-0
    
    # number of events that occurred in a trial
    num_event_null[i]<-sum(simulationResults$Status_null==2)
    num_event_ctrl[i]<-sum(simulationResults$Arm==1 & simulationResults$Status==2)
    num_event_trt[i]<-sum(simulationResults$Arm==2 & simulationResults$Status==2)
    
    # mean follow up time for a trial
    mean_follow_up_null[i]<-mean(simulationResults$Time_to_event_null)
    mean_follow_up_ctrl[i]<-mean(simulationResults[which(simulationResults$Arm==1),]$Time_to_event)
    mean_follow_up_trt[i]<-mean(simulationResults[which(simulationResults$Arm==2),]$Time_to_event)
  }
  
  return(list("null_rejections"=null_reject, "alternative_rejections"=alt_reject, "num_event_null"=num_event_null, 
              "num_event_ctrl"=num_event_ctrl, "num_event_trt"=num_event_trt, "mean_follow_up_null"=mean_follow_up_null, 
              "mean_follow_up_ctrl"=mean_follow_up_ctrl, "mean_follow_up_trt"=mean_follow_up_trt))
}



# print results given a scenario and scenario number
printResults<-function(scenario, scenario_num){
  cat("Scenario",scenario_num,"Alpha:", prop.table(table(scenario$null_rejections))[2], "\n")
  cat("Scenario",scenario_num,"Power:", prop.table(table(scenario$alternative_rejections))[2], "\n")
  cat("Scenario",scenario_num,"Mean Event Count: (null):", mean(scenario$num_event_null)/2, " (ctrl):", mean(scenario$num_event_ctrl), " (trt):", mean(scenario$num_event_trt), "\n")
  cat("Scenario",scenario_num,"Mean FU time: (null):", mean(scenario$mean_follow_up_null), " (ctrl):", mean(scenario$mean_follow_up_ctrl), " (trt):", mean(scenario$mean_follow_up_trt), "\n")
}



## Run five different scenarios
# scenario 1 - base study design - 12 months of follow up after last patient is recruited
scenario_1<-runSimulations(recruitment_rate=recruitment_rate, null_event_rate=ner, lower_event_rate=ler, upper_event_rate=uer, LTFU_rate=ltfur, FU_time=70)
printResults(scenario_1, 1)

# scenario 2 - additional year of follow-up after all patients are accrued (total of 24 months)
scenario_2<-runSimulations(recruitment_rate=recruitment_rate,null_event_rate=ner, lower_event_rate=ler, upper_event_rate=uer, LTFU_rate=ltfur, FU_time=82)
printResults(scenario_2, 2)

# scenario 3 - additional 100 participants - 12 months of follow up after last patient is recruited
recruitment_rate_s3<-c(recruitment_rate, seq(870, 970, 17)[-1], 970)
scenario_3<-runSimulations(recruitment_rate=recruitment_rate_s3, null_event_rate=ner, lower_event_rate=ler, upper_event_rate=uer, LTFU_rate=ltfur, FU_time=70)
printResults(scenario_3, 3)

# scenario 4 - only 9 months of follow up after last patient is recruited
scenario_4<-runSimulations(recruitment_rate=recruitment_rate,null_event_rate=ner, lower_event_rate=ler, upper_event_rate=uer, LTFU_rate=ltfur, FU_time=67)
printResults(scenario_4, 4)

# scenario 5 - each patient is followed for 12 months max, with 9 months of follow up after last patient is recruited
scenario_5<-runSimulations(recruitment_rate=recruitment_rate,null_event_rate=ner, lower_event_rate=ler, upper_event_rate=uer, LTFU_rate=ltfur, FU_time=67, max_patient_FU_time=12)
printResults(scenario_5, 5)