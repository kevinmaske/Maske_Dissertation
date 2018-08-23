##############################################
##############################################
############## Data Preparation ##############
##############################################
##############################################

library(tidyverse)
library(JM)
library(data.table)

## Function used to adhere to the long format needed by package JM

longify <- function(data.f, gran = 4, Static_Dat = logdata){
  
  # Aggregate data into hours first
  long_set <- data.f %>% group_by(icustay_id,
                                  obstime = as.numeric(floor((time/gran))*4)) %>%
    summarize(val = mean(valuenum))
  
  # Start time as just a copy of the obstime
  long_set$start <- long_set$obstime
  
  
  # Append each subjects end time,
  # the minimum of death time, discharge time, or 72
  long_set <- merge(long_set,
                    TIMES[,c('icustay_id', 'dischtime_hours', 'deathtime_hours')],
                    by = 'icustay_id', all.x = TRUE)
  
  long_set <- long_set %>% mutate(Time = pmin(dischtime_hours,
                                              deathtime_hours,
                                              72, na.rm = TRUE))
  
  # a list of the relevant times for each patient
  sub_times <- long_set %>% group_by(icustay_id) %>%
    summarise(times = list(unique(c(start,
                                    dischtime_hours,
                                    deathtime_hours,
                                    72))))  
  
  long_set <- merge(long_set,
                    sub_times,
                    by='icustay_id',
                    all.x = TRUE)
  
  end <- NULL
  
  for(i in 1:nrow(long_set)){
    ls_temp <- long_set$times[[i]]
    end[i] <- ls_temp[ls_temp > long_set$obstime[i]] %>% min(na.rm = TRUE)
  }
  
  
  long_set$end <- end
  
  case <- NULL
  
  for(i in 1:nrow(long_set)){
    if(is.na(long_set$deathtime_hours[i]) |
       long_set$deathtime_hours[i] > long_set$end[i]){
      case[i] <- 0
    } else {
      case[i] <- 1
    }
  }
  
  long_set$case <- case
  
  final_set <- merge(long_set[,c('icustay_id', 'T', 'start', 'end',
                                 'obstime', 'val', 'case')],
                     Static_Dat[,c('icustay_id', 'gender', 'race',
                                   'service', 'age')],
                     all.x = TRUE)
  
  
  return(final_set)
}



# load data
load("NA_DYNAMIC.RData")  # From the previous file
ICUSTAYS <- fread("ICUSTAYS.CSV")

# Remove outliers and invalid observations
NA_HR.f <- NA_HR %>% filter(time > 0, time < 72,
                            valuenum >= 22, valuenum <= 257,
                            icustay_id %in% logdata$icustay_id)


NA_SBP.f <- NA_SBP %>% filter(time > 0, time < 72,
                              valuenum >= 33, valuenum <= 341,
                              icustay_id %in% logdata$icustay_id)

NA_DBP.f <- NA_DBP %>% filter(time > 0, time < 72,
                              valuenum >= 5, valuenum <= 308,
                              icustay_id %in% logdata$icustay_id)

NA_RR.f <- NA_RR %>% filter(time > 0, time < 72,
                            valuenum >= 2, valuenum <= 91,
                            icustay_id %in% logdata$icustay_id)

NA_O2S.f <- NA_O2S %>% filter(time > 0, time < 72,
                              valuenum >= 70, valuenum <= 100,
                              icustay_id %in% logdata$icustay_id)

NA_HR.f <- merge(NA_HR.f,
                 TIMES[,c('icustay_id', 'dischtime_hours', 'deathtime_hours')],
                 by = 'icustay_id', all.x = TRUE)

NA_HR.f <- NA_HR.f %>% filter(time < pmin(dischtime_hours,
                                          deathtime_hours, na.rm = TRUE))


NA_SBP.f <- merge(NA_SBP.f,
                  TIMES[,c('icustay_id', 'dischtime_hours', 'deathtime_hours')],
                  by = 'icustay_id', all.x = TRUE)
NA_SBP.f <- NA_SBP.f %>% filter(time < pmin(dischtime_hours,
                                            deathtime_hours, na.rm = TRUE))

NA_DBP.f <- merge(NA_DBP.f,
                  TIMES[,c('icustay_id', 'dischtime_hours', 'deathtime_hours')],
                  by = 'icustay_id', all.x = TRUE)
NA_DBP.f <- NA_DBP.f %>% filter(time < pmin(dischtime_hours,
                                            deathtime_hours, na.rm = TRUE))

NA_RR.f <- merge(NA_RR.f,
                 TIMES[,c('icustay_id', 'dischtime_hours', 'deathtime_hours')],
                 by = 'icustay_id', all.x = TRUE)
NA_RR.f <- NA_RR.f %>% filter(time < pmin(dischtime_hours,
                                          deathtime_hours, na.rm = TRUE))

NA_O2S.f <- merge(NA_O2S.f,
                  TIMES[,c('icustay_id', 'dischtime_hours', 'deathtime_hours')],
                  by = 'icustay_id', all.x = TRUE)
NA_O2S.f <- NA_O2S.f %>% filter(time < pmin(dischtime_hours,
                                            deathtime_hours, na.rm = TRUE))


# Get subjects that have measurements per variable
subs_HR <- NA_HR.f$icustay_id %>% unique
subs_SBP <- NA_SBP.f$icustay_id %>% unique
subs_DBP <- NA_DBP.f$icustay_id %>% unique
subs_RR <- NA_RR.f$icustay_id %>% unique
subs_O2S <- NA_O2S.f$icustay_id %>% unique

# Get cohort for all
woTEMP <- Reduce(intersect, list(subs_HR,
                                 subs_SBP,
                                 subs_DBP,
                                 subs_RR,
                                 subs_O2S))

length(woTEMP)
# 14,113


# We get id's of those who are in each case group (<X hours, before dischtime)
# we get 24 and 48 as well since we are interested in checking these outcomes out
death72 <- death %>%
  filter(icustay_id %in% logdata$icustay_id,
         death > 0) %>%
  filter(deathtime_hours <= dischtime_hours) %>%
  filter(deathtime_hours <= 72)

death48 <- death %>%
  filter(icustay_id %in% logdata$icustay_id,
         death > 0) %>%
  filter(deathtime_hours <= dischtime_hours) %>%
  filter(deathtime_hours <= 48)

death24 <- death %>%
  filter(icustay_id %in% logdata$icustay_id,
         death > 0) %>%
  filter(deathtime_hours <= dischtime_hours) %>%
  filter(deathtime_hours <= 24)

# Prepare long format dataset
JM_HR <- NA_HR.f %>% filter(icustay_id %in% JM_Cohort)
JM_SBP <- NA_SBP.f %>% filter(icustay_id %in% JM_Cohort)
JM_DBP <- NA_DBP.f %>% filter(icustay_id %in% JM_Cohort)
JM_RR <- NA_RR.f %>% filter(icustay_id %in% JM_Cohort)
JM_O2S <- NA_O2S.f %>% filter(icustay_id %in% JM_Cohort)

JM_HR_long <- longify(JM_HR)
JM_SBP_long <- longify(JM_SBP)
JM_DBP_long <- longify(JM_DBP)
JM_RR_long <- longify(JM_RR)
JM_O2S_long <- longify(JM_O2S)




##############################################
##############################################
################ Joint Models ################
##############################################
##############################################

library(JM)
library(tidyverse)
library(data.table)
library(caret)
library(PRROC)
library(pROC)

# Obtain icustay_id's of the controls and cases
N <- nrow(JM_Static)

controls <- JM_Static %>% filter(!(icustay_id %in% death72$icustay_id))

# 80-20 for initial study, 10-fold CV later on
set.seed(1)
train.id.case <- createDataPartition(death72$icustay_id, p=0.8)$Resample1
train.id.ctrl <- createDataPartition(as.integer(controls$icustay_id), p=0.8)$Resample1
train.id <- c(death72$icustay_id[train.id.case],
              controls$icustay_id[train.id.ctrl])  # 80% training set
test.id <- setdiff(JM_Static$icustay_id, train.id)  # 20% development set

set.seed(2)
caseFolds <- createFolds(death72$icustay_id, k = 5)
ctrlFolds <- createFolds(controls$icustay_id, k = 5)

k = 5
random.classifier.prauc <- length(train.id.case)/length(train.id)

###################################################################
###################################################################
### To illustrate the code structure, below is the code for the ###
### heart rate joint model.                                     ###
###################################################################
###################################################################

# Attaching case or not to static dataset
case <- rep(0, nrow(JM_HR_long))
case[which(JM_HR_long$icustay_id %in% death72$icustay_id)] = 1
JM_HR_long$death <- case
JM_HR_long$Time <- JM_HR_long$'T'  # add a duplicate, better named column to avoid issues

survTable.HR <- JM_HR_long[!duplicated(JM_HR_long$icustay_id),]

# We work with the training data first, we use most default settings for most of this
set.seed(1)
lme.HR <- lme(val ~ obstime,
              random = ~ obstime | icustay_id,
              data = filter(JM_HR_long,  # Get only those in the training set
                            icustay_id %in% train.id))

# We specify all our static variables
cox.HR <- coxph(Surv(Time, death) ~ gender + race + service + age,
                data = filter(survTable.HR,
                              icustay_id %in% train.id),
                x = TRUE)

# Fit the joint model

jointHR <- jointModel(lme.HR,
                      cox.HR,
                      timeVar = 'obstime',
                      method = 'spline-PH-aGH')
# We select the smoother spline method for baseline risk.

test.HR <- filter(JM_HR_long, icustay_id %in% test.id)

HR.pred <- ((survfitJM(jointHR,
                       newdata = test.HR,
                       idVar = 'icustay_id',
                       simulate = FALSE,
                       survTimes = 72)$summaries) %>% unlist)[c(FALSE,TRUE)]


# c(FALSE, TRUE) removes odd-indexed entries, which is fine because in the unlisted
# version, the odd entries just contain 72, the survTimes. We are interested in the
# predicted survival probabilities, which are in the even entries.

validation.HR <- data.frame(icustay_id = unique(test.HR$icustay_id),
                            caseProb = 1 - HR.pred,
                            case = 0)

validation.HR$case[which(validation.HR$icustay_id %in% death72$icustay_id)] <- 1

roc.HR <- roc.curve(scores.class0 = validation.HR$caseProb,
                    weights.class0 = validation.HR$case,
                    curve = TRUE)


plot(roc.HR)

pr.HR <- pr.curve(scores.class0 = validation.HR$caseProb,
                  weights.class0 = validation.HR$case,
                  curve = TRUE)

plot(pr.HR)

## 5-fold CV
cv.HR.lme <- NULL
cv.HR.cox <- NULL
cv.HR.JM <- NULL
cv.HR.pred <- NULL
cv.HR.val <- NULL
cv.HR.roc <- NULL
cv.HR.pr <- NULL

set.seed(1)
for(fold in 1:k){
  cv.HR.lme <- lme(val ~ obstime,
                   random = ~ obstime | icustay_id,
                   data = filter(JM_HR_long,
                                 !(icustay_id %in% folds[[fold]])))
  cv.HR.cox <- coxph(Surv(Time, death) ~ gender + race +
                       service + age,
                     data = filter(survTable.HR,
                                   !(icustay_id %in% folds[[fold]])),
                     x = TRUE)
  cv.HR.JM <- jointModel(cv.HR.lme,
                         cv.HR.cox,
                         timeVar = 'obstime',
                         method = 'spline-PH-aGH')
  # Have it overwrite at each loop in order to conserve space
  
  # The model is fit above, now we make predictions.
  cv.HR.pred[[fold]] <- ((survfitJM(cv.HR.JM,
                                    newdata = filter(JM_HR_long,
                                                     icustay_id %in% folds[[fold]]),
                                    idVar = 'icustay_id',
                                    simulate = FALSE,
                                    survTimes = 72)$summaries) %>% unlist)[c(FALSE,TRUE)]
  
  
  
  cv.HR.val[[fold]] <- data.frame(icustay_id = unique(filter(JM_HR_long,
                                                             (icustay_id %in% folds[[fold]]))$icustay_id),
                                  caseProb = 1 - cv.HR.pred[[fold]],
                                  case = 0)
  
  cv.HR.val[[fold]]$case[which(cv.HR.val[[fold]]$icustay_id %in%
                                 death72$icustay_id)] <- 1
  
  # Get ROC and PR-AUC
  cv.HR.roc[fold] <- roc.curve(scores.class0 = cv.HR.val[[fold]]$caseProb,
                               weights.class0 = cv.HR.val[[fold]]$case)$auc
  cv.HR.pr[fold] <- pr.curve(scores.class0 = cv.HR.val[[fold]]$caseProb,
                             weights.class0 = cv.HR.val[[fold]]$case)$auc.integral
}

# Get CV-Metrics
mean(cv.HR.roc) %>% round(4)
sd(cv.HR.roc) %>% round(4)

mean(cv.HR.pr) %>% round(4)
sd(cv.HR.pr) %>% round(4)


####################################################################
####################################################################
### Below, the code for creating the balanced subsample is shown ###
####################################################################
####################################################################

controls <- filter(JM_Static, !(icustay_id %in% death72$icustay_id))$icustay_id

set.seed(1)
down_controls <- NULL
for(i in 1:5){
  down_controls[[i]] <- sample(controls, 169)
}

cases <- death72$icustay_id

downsample_id <- NULL
for(i in 1:5){
  downsample_id[[i]] <- c(down_controls[[i]], cases)
}

# Get total downsampleids
sample_ds <- Reduce(union, downsample_id)

JM_HR_down <- filter(JM_HR_long, icustay_id %in% sample_ds)
JM_SBP_down <- filter(JM_SBP_long, icustay_id %in% sample_ds)
down_Static <- filter(JM_Static, icustay_id %in% sample_ds)



#######################################################
#######################################################
### Below, the code for model diagnostics is shown. ###
#######################################################
#######################################################

## Longitudinal Residuals

plot(jointHR)

## Martingale Residuals (Survival Process Diagnostics)
plotResid <- function(x, y, col.loess = 'black', ...){
  plot(x, y, ...)
  lines(lowess(x,y), col = col.loess, lwd = 2)
  abline(h = 0, lty = 3, col = 'grey', lwd = 2)
}

martResDown <- residuals(jointDown, process = 'Event')
mi.tDown <- fitted(jointDown, process = 'Longitudinal',
                   type = 'EventTime')

plotResid(mi.tDown, martResDown, col.loess = 'grey62',
          ylab = 'Martingal Residuals',
          xlab = 'Fitted Values on HR',
          main = 'Martingale Residuals for Subsample')



