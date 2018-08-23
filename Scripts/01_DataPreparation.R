#### Libraries ###
library(tidyverse)
library(data.table)
library(R.utils)

### Functions ###
# This function will assume no empty dates; we have to check that before we run it
# Returns a vector of "times"
# Default TIMES refers to the above constructed dataset
time_fmt <- "%Y-%m-%d %H:%M:%S"

TimeFn <- function(dataset, times = TIMES){
  robust_frame <- merge(dataset,
                        TIMES[,c('icustay_id', 'INTIME')],
                        by = 'icustay_id',
                        all.x = TRUE)
  
  delta_t <- difftime(as.POSIXlt(robust_frame$charttime, format = time_fmt),
                      as.POSIXlt(robust_frame$INTIME, format = time_fmt),
                      units = 'hours') %>% round(4)
  
  return(delta_t)
}

# This function aggregates measurements based on specified bin lengths
Binned <- function(NA_dataset, gran = 4){
  binned <- NA_dataset %>%
    group_by(bin = floor((time/gran)+1)) %>%
    summarise(cohort = list(unique(icustay_id)))
}

# This function is the "leniency" criteria that allows for some missing data
lenient <- function(cohort.bins, max.bins=4){
  subjects <- list()
  
  subjects[[1]] <- Reduce(union, cohort.bins$cohort[c(1:2)])
  
  for(i in 2:max.bins){  # We have an extra bin for checking last interpolation
    subjects[[i]] <- Reduce(union,
                            cohort.bins$cohort[c(i:(i+1))]) %>%
      intersect(subjects[[(i-1)]])
  }
  
  return(subjects)
}

# This is the function that interpolates between longitudinal observations
Interpol <- function(L_Data){
  L_Agg <- L_Data %>%
    group_by(bin = floor((time/4)+1), icustay_id) %>%  # Aggregate Data by bin
    summarise(val = mean(valuenum))
  
  L_tab <- xtabs(val ~ icustay_id + bin,  # Crosstabulate into serpate columns per bin
                 data = L_Agg) %>% as.data.frame.matrix
  
  L_imp <- L_tab
  # Impute for bins 2 - 4
  for(i in 2:4){
    na.idx <- which(L_tab[,i] == 0)  # get indexes of missing entries
    L_imp[na.idx,i] <- (L_tab[na.idx,(i-1)] + L_tab[na.idx,(i+1)])/2  # Interpolate
  }
  
  # Extrapolate bin 1
  na.idx1 <- which(L_tab[,1] == 0)
  L_imp[na.idx1,1] = 2*L_imp[na.idx1,2] - L_imp[na.idx1,3]  # Extrapolate
  
  # Add icustay_id at end
  L_imp[,6] <- rownames(L_imp)
  
  
  # We return only bins 1-4, exclude 5 from this point on. It was only there for imputaion.
  return(L_imp[,c(1:4,6)])
}

## Data Load
# Read the files
All_DBP <- fread('AllCohortDBP.csv')
All_HR <- fread('AllCohortHR.csv')
All_O2S <- fread('AllCohortO2S.csv')
All_PH <- fread('AllCohortPH.csv')
All_RR <- fread('AllCohortRR.csv')
All_TEMP <- fread('AllCohortTEMP.csv')
All_TIME <- fread('AllCohortTIME.csv')
All_WBC <- fread('AllCohortWBC.csv')
All_SBP <- fread('AllCohortSBP.csv')

# Load static Data Points and Death
cohort <- fread('df_cohort.csv')
death <- fread('df_death.csv')
static <- fread('df_static_data.csv')

# Label Time Data
names(All_TIME) <- c("subject_id", "icustay_id", "intime", "outtime")

# Dynamic Data labelling
names(All_HR) <- c("subject_id", "hadm_id", "icustay_id", "itemid", 
                   "charttime", "storetime", "valuenum", "valueuom")
names(All_DBP) <- c("subject_id", "hadm_id", "icustay_id", "itemid", 
                    "charttime", "storetime", "valuenum", "valueuom")
names(All_O2S) <- c("subject_id", "hadm_id", "icustay_id", "itemid", 
                    "charttime", "storetime", "valuenum", "valueuom")
names(All_PH) <- c("subject_id", "hadm_id", "icustay_id", "itemid", 
                   "charttime", "storetime", "valuenum", "valueuom")
names(All_RR) <- c("subject_id", "hadm_id", "icustay_id", "itemid", 
                   "charttime", "storetime", "valuenum", "valueuom")
names(All_SBP) <- c("subject_id", "hadm_id", "icustay_id", "itemid", 
                    "charttime", "storetime", "valuenum", "valueuom")
names(All_TEMP) <- c("subject_id", "hadm_id", "icustay_id", "itemid", 
                     "charttime", "storetime", "valuenum", "valueuom")
names(All_WBC) <- c("subject_id", "itemid", "charttime", "valuenum", "valueuom")

##############################################
##############################################
########## Cohort from G-M's Study ###########
##############################################
##############################################

nrow(cohort)  # Count total unfiltered number of recorded ICU stays
# 61,532 recorded stays

JohnsonCohort <- subset(cohort,
                        subset = (exclusion_over_15 == 0) &
                          (exclusion_valid_data == 0) &
                          (exclusion_stay_lt_4hr == 0) &
                          (exclusion_organ_donor == 0))

nrow(JohnsonCohort)  # Sample size of Johnson Cohort
# 52,085 samples


NewCohort <- JohnsonCohort[(!duplicated(JohnsonCohort$subject_id) &
                              age<89), c("subject_id","hadm_id","icustay_id","intime",
                                         "outtime", "age","gender", "ethnicity", "icu_los",
                                         "hosp_los","death_icu") ]

nrow(NewCohort)  # Samples in New Cohort
# 36,227 samples remain

# Merge with available static data
Static_Cohort <- merge(NewCohort, static, by = 'icustay_id')
nrow(Static_Cohort)
# 36, 198 samples for final static cohort


##############################################
##############################################
########## Check Units in Datasets ###########
##############################################
##############################################

unique(Static_DBP$valueuom)
# 'mmHg'

unique(Static_HR$valueuom)
# 'BPM' 'bpm'

unique(Static_O2S$valueuom)
# '%' ''

unique(Static_RR$valueuom)
# '' 'BPM' 'insp/min'

unique(Static_SBP$valueuom)
# '' 'mmHg'

# SBP
nrow(Static_SBP[valueuom == ""])
unique(Static_SBP[valueuom == ""]$valuenum)  # All Null's
# 7,944 rows, All NA's

# O2 Saturation
nrow(Static_O2S[valueuom == ""])
unique(Static_O2S[valueuom == ""]$valuenum)
# 453 rows, All NA's

# Respiration Rate
nrow(Static_RR[valueuom == ""])
unique(Static_RR[valueuom == ""]$valuenum)
# 3124 Rows, all NA's

# The prefix 'NU' is used to indicate the removal of Null Units
NU_HR <- Static_HR
NU_DBP <- Static_DBP
NU_SBP <- Static_SBP[!(Static_SBP$valueuom == ""),]  # Remove null units
NU_O2S <- Static_O2S[!(Static_O2S$valueuom == ""),]
NU_RR <- Static_RR[!(Static_RR$valueuom == ""),]



# Uniform Naming Scheme for data so far
# Considerations: All these new data-sets are taking up more memory; is there a way to retain all these annotations while minimizing memory usage?

NA_DBP <- NU_DBP
NA_HR <- NU_HR
NA_O2S <- NU_O2S
NA_SBP <- NU_SBP
NA_RR <- NU_RR

NA_DBP$valuenum <- as.numeric(NA_DBP$valuenum)
NA_HR$valuenum <- as.numeric(NA_HR$valuenum)
NA_O2S$valuenum <- as.numeric(NA_O2S$valuenum)
NA_SBP$valuenum <- as.numeric(NA_SBP$valuenum)
NA_RR$valuenum <- as.numeric(NA_RR$valuenum)


# There are some discrepancies in different tabels re. admission time;
# we specify a definitive guide to time using the ICUSTAYS table
Static_Cohort$icustay_id


TIMES <- merge(x = Static_Cohort[,c('icustay_id', 'subject_id.x')],
               y = ICUSTAYS[,c('ICUSTAY_ID', 'INTIME', 'OUTTIME')],
               by.x = 'icustay_id',
               by.y = 'ICUSTAY_ID')
TIMES <- merge(TIMES,
               death,
               by = "icustay_id")


### Outlier Removal, removal of data before and after times of interest (0-72 hours) ###
NA_HR.f <- NA_HR %>% filter(time > 0, time < 20,
                            valuenum >= 22, valuenum <= 257)


NA_SBP.f <- NA_SBP %>% filter(time > 0, time < 20,
                              valuenum >= 33, valuenum <= 341)

NA_DBP.f <- NA_DBP %>% filter(time > 0, time < 20,
                              valuenum >= 5, valuenum <= 308)

NA_RR.f <- NA_RR %>% filter(time > 0, time < 20,
                            valuenum >= 2, valuenum <= 91)

NA_O2S.f <- NA_O2S %>% filter(time > 0, time < 20,
                              valuenum >= 70, valuenum <= 100)



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


# Get cohorts for each bin of each variable
Cohort_HR <- Binned(NA_HR.f)
Cohort_RR <- Binned(NA_RR.f)
Cohort_DBP <- Binned(NA_DBP.f)
Cohort_SBP <- Binned(NA_SBP.f)
Cohort_O2S <- Binned(NA_O2S.f)

# Apply leniency criterion
Lenient_HR <- lenient(Cohort_HR)
Lenient_RR <- lenient(Cohort_RR)
Lenient_DBP <- lenient(Cohort_DBP)
Lenient_SBP <- lenient(Cohort_SBP)
Lenient_O2S <- lenient(Cohort_O2S)

# Get logistic cohort
Logistic_Cohort <- Reduce(intersect, list(Lenient_HR[[4]],
                                          Lenient_RR[[4]],
                                          Lenient_DBP[[4]],
                                          Lenient_SBP[[4]],
                                          Lenient_O2S[[4]]))



### Performing imputations ###

# Get the datasets for the Logistic Cohort
# L stands for logistic cohort
L_HR <- NA_HR.f %>%
  filter(icustay_id %in% Logistic_Cohort)
L_SBP <- NA_SBP.f %>%
  filter(icustay_id %in% Logistic_Cohort)
L_DBP <- NA_DBP.f %>%
  filter(icustay_id %in% Logistic_Cohort)
L_O2S <- NA_O2S.f %>%
  filter(icustay_id %in% Logistic_Cohort)
L_RR <- NA_O2S.f %>%
  filter(icustay_id %in% Logistic_Cohort)

# Interpolate and compile into single table
c_HR <- Interpol(L_HR)
c_DBP <- Interpol(L_DBP)
c_SBP <- Interpol(L_SBP)
c_RR <- Interpol(L_RR)
c_O2S <- Interpol(L_O2S)

colnames(c_HR) <- c('HR.1', 'HR.2', 'HR.3', 'HR.4', 'icustay_id')
colnames(c_DBP) <- c('DBP.1', 'DBP.2', 'DBP.3', 'DBP.4', 'icustay_id')
colnames(c_SBP) <- c('SBP.1', 'SBP.2', 'SBP.3', 'SBP.4', 'icustay_id')
colnames(c_RR) <- c('RR.1', 'RR.2', 'RR.3', 'RR.4', 'icustay_id')
colnames(c_O2S) <- c('O2S.1', 'O2S.2', 'O2S.3', 'O2S.4', 'icustay_id')

# Check static variables that'll be included
L_Static <- Static_Cohort %>% filter(icustay_id %in% Logistic_Cohort)

# Compile categorical vector of race
L_Race <- with(L_Static,  # arranged by icustay id occurrence in static table
               cbind(race_black,
                     race_hispanic,
                     race_asian,
                     race_other))
L_Race %>% colnames
L_Race
L_Race_v <- NULL
for(i in 1:nrow(L_Race)){
  if(sum(L_Race[i,] > 0)){
    L_Race_v[i] <- which(L_Race[i,] > 0)
  } else {
    L_Race_v[i] <- 0
  }
}

L_Race_v <- L_Race_v %>% as.factor
levels(L_Race_v) <- c("race_white", "race_black",
                      "race_hispanic", "race_asian",
                      "race_other")



# Compile categorical vector of services
Services <- with(L_Static,
                 cbind(service_med,
                       service_cmed,
                       service_nmed,
                       service_nsurg,
                       service_csurg,
                       service_surg,
                       service_traum))

max(rowSums(Services))  # nobody has 2 entries

L_Serv_v <- NULL
for(i in 1:nrow(Services)){
  if(sum(Services[i,] > 0)){
    L_Serv_v[i] <- which(Services[i,] > 0)
  } else {
    L_Serv_v[i] <- 0
  }
}

L_Serv_v <- as.factor(L_Serv_v)
levels(L_Serv_v) <- c("others",
                      "med", "cmed",
                      "nmed", "nsurg",
                      "csurg", "surg",
                      "traum")  # This is appropriate.

# build the single dataset
L_Static$gender <- L_Static$gender %>% as.factor
L_Static$race <- L_Race_v
L_Static$service <- L_Serv_v
L_Static$age <- L_Static$age.x - mean(L_Static$age.x) # note the centering



# Compiling a dataframe containing information re. being a case
case.idx <- which(TIMES$icustay_id %in% death72$icustay_id)
V_death_vec <- rep(0, nrow(TIMES))
V_death_vec[case.idx] <- 1  # These are the cases

TIMES$case <- V_death_vec

## Marge all
L_Data <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2,
                                            by = "icustay_id", all.x = TRUE),
                 list(c_HR, c_DBP, c_SBP, c_RR, c_O2S,
                      L_Static[,c('icustay_id', 'age', 'gender',
                                  'service', 'race')],
                      TIMES[,c('icustay_id', 'INTIME', 'OUTTIME',
                               'case')]))

## Save the datafile
saveRDS(L_Data, file = "LogisticDataD.rds")


