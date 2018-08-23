library(JMbayes)
library(tidyverse)
library(data.table)
library(parallel)
library(caret)
library(PRROC)
library(reshape2)

##############################################
##############################################
############## Data Preparation ##############
##############################################
##############################################

# Load the logistic dataset (same data will be used)
down_Static <- readRDS('LogisticDataD.rds')


# Reshaping the datasets for HR and SBP
HR_wide <- down_Static[,c('icustay_id', 'HR.1', 'HR.2', 'HR.3', 'HR.4')]
colnames(HR_wide) <- c('icustay_id', 0, 4, 8, 12)
HR_long <- gather(HR_wide,
                  key = obstime,
                  value = HR,
                  -'icustay_id')

SBP_wide <- down_Static[,c('icustay_id', 'SBP.1', 'SBP.2', 'SBP.3', 'SBP.4')]
colnames(SBP_wide) <- c('icustay_id', 0, 4, 8, 12)
SBP_long <- gather(SBP_wide,
                   key = obstime,
                   value = SBP,
                   -'icustay_id')

JMB_Data <- merge(HR_long, SBP_long,
               by = c('icustay_id', 'obstime'))

JMB_Data$obstime <- as.numeric(JMB_Data$obstime)


Static_info <- JM_HR_down[!duplicated(JM_HR_down$icustay_id),]
JMB_Data1 <- merge(JMB_Data,
                   Static_info[,c('icustay_id', 'T',
                                  'gender', 'race', 'service', 'age')],
                   by = 'icustay_id',
                   all.x = TRUE)

# Add case indicator
JMB_Data1$death <- 0
JMB_Data1$death[which(JMB_Data1$icustay_id %in% death72$icustay_id)] <- 1

JMB_Data1$Time <- JMB_Data1$'T'  # For ease of referencing the column

### Creating the 5 balanced subsamples

# Obtain icustay_id's of the controls and cases
N <- 169*2

controls <- NULL
for(i in 1:5){
  controls[[i]] <- as.integer(down_controls[[i]])
}

downs <- NULL
for(i in 1:5){
  downs[[i]] <- as.integer(downsample_id[[i]])
}

train.id.case <- NULL
train.id.ctrl <- NULL
train.id <- NULL
test.id <- NULL
set.seed(1)
for(i in 1:5){
  train.id.case[[i]] <- createDataPartition(cases, p=0.8)$Resample1
  train.id.ctrl[[i]] <- createDataPartition(controls[[i]], p=0.8)$Resample1
  train.id[[i]] <- c(cases[train.id.case[[i]]],
                     controls[[i]][train.id.ctrl[[i]]])
  test.id[[i]] <- setdiff(downs[[i]], train.id[[i]])
}


##############################################
##############################################
########## Joint Models (Bayesian) ###########
##############################################
##############################################

survTable <- JMB_Data1[!duplicated(JMB_Data1$icustay_id),]

CoxBasic <- NULL

set.seed(1)
for(i in 1:5){
  CoxBasic[[i]] <- coxph(Surv(Time, death) ~ age + gender + race + service,
                         data =filter(survTable,
                                      icustay_id %in% train.id[[i]]),
                         model = TRUE)
}

# MV LME model
set.seed(1)
BasicMixed <- NULL

for(i in 1:5){
  BasicMixed[[i]] <- mvglmer(list(HR ~ obstime + (obstime | icustay_id),
                                  SBP ~ obstime + (obstime | icustay_id)),
                             data = filter(JMB_Data1, icustay_id %in% train.id[[i]]),
                             families = list(gaussian, gaussian))
}


# Basic assumption that the longitudinal variables are gaussian

set.seed(1)
jointBasic <- NULL

for(i in 1:5){
  jointBasic[[i]] <- mvJointModelBayes(BasicMixed[[i]],
                                       CoxBasic[[i]],
                                       timeVar = 'obstime')
  
}

## Check predictive power

test.data <- NULL
for(i in 1:5){
  test.data[[i]] <- filter(JMB_Data1, icustay_id %in% test.id[[i]])
  test.data[[i]] <- test.data[[i]][order(test.data[[i]]$icustay_id,
                                         test.data[[i]]$obstime),]
}

pred.Basic <- NULL

set.seed(1)
for(i in 1:5){
  pred.Basic[[i]] <- survfitJM(jointBasic[[i]],
                               newdata = test.data[[i]],
                               type = "SurvProb",
                               idVar = 'icustay_id',
                               survTimes = 72,
                               seed = 1)
  
}

fit.Basic <- NULL
validation.Basic <- NULL
for(i in 1:5){
  fit.Basic[[i]] <- (pred.Basic[[i]]$summaries %>%
                       unlist)[c(FALSE, TRUE, FALSE, FALSE, FALSE)]
  
  validation.Basic[[i]] <- data.frame(icustay_id = unique(test.data[[i]]$icustay_id),
                                      caseProb = 1 - fit.Basic[[i]],
                                      case = 0)
  
  validation.Basic[[i]][which(validation.Basic[[i]]$icustay_id %in% cases),]$case <- 1
  
}

jm.PR <- NULL
jm.ROC <- NULL

for(i in 1:5){
  jm.PR[[i]] <- pr.curve(scores.class0 = validation.Basic[[i]]$caseProb,
                         weights.class0 = validation.Basic[[i]]$case,
                         curve = TRUE)
  jm.ROC[[i]] <- roc.curve(scores.class0 = validation.Basic[[i]]$caseProb,
                           weights.class0 = validation.Basic[[i]]$case,
                           curve = TRUE)
}

jm.pr.auc <- NULL
jm.roc.auc <- NULL

for(i in 1:5){
  jm.pr.auc[i] <- jm.PR[[i]]$auc.integral
  jm.roc.auc[i] <- jm.ROC[[i]]$auc
}
mean(jm.roc.auc) %>% round(4)
sd(jm.roc.auc) %>% round(4)

mean(jm.pr.auc) %>% round(4)
sd(jm.pr.auc) %>% round(4)