### libraries
library(tidyverse)
library(PRROC)
library(caret)
library(MASS)

### load data
data <- readRDS('LogisticDataD.rds')


##############################################
##############################################
############ Cumulative Approach #############
##############################################
##############################################

# get indices of the cases, and indices of the controls
cases <- which(data$case == 1)
controls <- which(data$case == 0)


set.seed(1)
train.case.idx <- createDataPartition(cases, p=0.8)$Resample1
train.ctrl.idx <- createDataPartition(controls, p=0.8)$Resample1
train.idx <- c(cases[train.case.idx],  # Combine
               controls[train.ctrl.idx]) 

## Static Only
NULL_MODEL <- glm(case ~ 1,
                  data = data[train.idx,],
                  family = "binomial")

FULL_MODEL_S <- glm(case ~ age + gender + service + race,
                    data = data[train.idx,],
                    family = "binomial")

FW_Static <- stepAIC(NULL_MODEL,
                     scope = formula(FULL_MODEL_S),
                     direction = "forward") 
BW_Static <- stepAIC(FULL_MODEL_S,
                     direction = "backward") 
STEP_STATIC <- stepAIC(NULL_MODEL,
                       scope = formula(FULL_MODEL_S),
                       direction = "both")
Model_Static <- glm(case ~ service + age,
                    data = data[train.idx,],
                    family = "binomial")

testpred_Static <- predict(Model_Static, newdata = data[-train.idx,],
                           type = 'response')
roc_Static <- roc.curve(scores.class0 = testpred_Static,
                        weights.class0 = data$case[-train.idx],
                        curve = TRUE)
plot(roc_Static)
pr_Static <- pr.curve(scores.class0 = testpred_Static,
                      weights.class0 = data$case[-train.idx],
                      curve = T)
plot(pr_Static)


### Prep for 10-fold CV
set.seed(2)
n.folds <- 10
caseFolds <- createFolds(cases, k=n.folds)
ctrlFolds <- createFolds(controls, k=n.folds)
folds <- NULL

for(i in 1:n.folds){
  folds[[i]] <- c(cases[caseFolds[[i]]],
                  controls[ctrlFolds[[i]]])
}

cv.Static <- NULL
cv.pStatic <- NULL
cv.aucStatic <- NULL
cv.paucStatic <- NULL

for(fold in 1:n.folds){
  cv.train.idx <- setdiff(1:nrow(data), folds[[fold]])  # get training indices
  cv.Static[[fold]] <- glm(case ~ service + age,  # fit model
                           data = data, subset = cv.train.idx,
                           family = "binomial")
  
  # Perform testing
  test.idx <- folds[[fold]]
  cv.pStatic[[fold]] <- predict(cv.Static[[fold]], newdata = data[test.idx,],
                                type = 'response')
  cv.aucStatic[fold] <- roc.curve(scores.class0 = cv.pStatic[[fold]],
                                  weights.class0 = data$case[test.idx])$auc
  cv.paucStatic[fold] <- pr.curve(scores.class0 = cv.pStatic[[fold]],
                                  weights.class0 = data$case[test.idx])$auc.integral
}

mean(cv.aucStatic) %>% round(4)
sd(cv.aucStatic) %>% round(4)

mean(cv.paucStatic) %>% round(4)
sd(cv.paucStatic) %>% round(4)


## S + B1
FULL_MODEL_B1 <- glm(case ~ age + gender + service + race +
                       HR.1 + SBP.1 + DBP.1 + RR.1 + O2S.1,
                     data = data[train.idx,],
                     family = "binomial")
FW_B1 <- stepAIC(NULL_MODEL,
                 scope = formula(FULL_MODEL_B1),
                 direction = "forward") 
BW_B1 <- stepAIC(FULL_MODEL_B1,
                 direction = "backward") 
STEP_B1 <- stepAIC(NULL_MODEL,
                   scope = formula(FULL_MODEL_B1),
                   direction = "both")

Model_B1 <- glm(case ~ SBP.1 + HR.1 + service + age + DBP.1 + RR.1,
                data = data[train.idx,],
                family = "binomial")

summary(Model_B1)

testpred_B1 <- predict(Model_B1, newdata = data[-train.idx,],
                       type = 'response')

roc_B1 <- roc.curve(scores.class0 = testpred_B1,
                    weights.class0 = data$case[-train.idx],
                    curve = TRUE)
plot(roc_B1)

pr_B1 <- pr.curve(scores.class0 = testpred_B1,
                  weights.class0 = data$case[-train.idx],
                  curve = T)
plot(pr_B1)

# 10-fold CV
cv.B1 <- NULL
cv.pB1 <- NULL
cv.aucB1 <- NULL
cv.paucB1 <- NULL

for(fold in 1:n.folds){
  cv.train.idx <- setdiff(1:nrow(data), folds[[fold]])  # get training indices
  cv.B1[[fold]] <- glm(formula(Model_B1),  # fit model
                       data = data, subset = cv.train.idx,
                       family = "binomial")
  
  # Perform testing
  test.idx <- folds[[fold]]
  cv.pB1[[fold]] <- predict(cv.B1[[fold]], newdata = data[test.idx,],
                            type = 'response')
  cv.aucB1[fold] <- roc.curve(scores.class0 = cv.pB1[[fold]],
                              weights.class0 = data$case[test.idx])$auc
  cv.paucB1[fold] <- pr.curve(scores.class0 = cv.pB1[[fold]],
                              weights.class0 = data$case[test.idx])$auc.integral
}

mean(cv.aucB1) %>% round(4)
cv.aucB1 %>% range %>% round(4)
sd(cv.aucB1) %>% round(4)

mean(cv.paucB1) %>% round(4)
cv.paucB1 %>% range %>% round(4)
sd(cv.paucB1) %>% round(4)


#######################################################################
#######################################################################
### All other models built using logistic regression follow an      ###
### identical code process as the two shown above, except with      ###
### different "full model" specifications, and different sets of    ###
### candidate variales. Below, it is shown how O2S was categorised. ###
#######################################################################
#######################################################################

O2S.cat <- (O2S.vals < 90)

data$O2S.cat1 <- O2S.cat[,1] %>% factor(labels = c('OK', 'LOW'))
data$O2S.cat2 <- O2S.cat[,2] %>% factor(labels = c('OK', 'LOW'))
data$O2S.cat3 <- O2S.cat[,3] %>% factor(labels = c('OK', 'LOW'))
data$O2S.cat4 <- O2S.cat[,4] %>% factor(labels = c('OK', 'LOW'))

#######################################################################
#######################################################################
### For models using categorical O2S, the variable O2S.cat1 is      ###
### specified rather than O2S.1. Similarly, O2S.cat2 instead of     ###
### O2S.2, so on and so forth.                                      ###
#######################################################################
#######################################################################
pro