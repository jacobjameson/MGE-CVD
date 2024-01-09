################################################################################
# AUTHOR:             JACOB JAMESON
# PURPOSE:            CONSTRUCT W1,W4,W5 Analytical Dataset
################################################################################

# Load packages necessary for preparing final data
library(tidyverse)
library(haven)
library(labelled)
library(scales)

# Create and prepare wave 1 and wave 4 data
source('src/prepare data/Prepare W1.R')
source('src/prepare data/Prepare W4.R')
source('src/prepare data/Prepare W5.R')

# Get SES data
ses_path <- '~/Sue Goldie Dropbox/Jacob Jameson/Add Health/Data Upload 7.2021/Constructed Files/Constructed SES Variables'
ses.data <- read_xpt(paste0(ses_path, '/conses4.xpt'))
names(ses.data) <- tolower(names(ses.data))

# Merge wave 1, wave 4, wave 5, and SES data together
final.df <- merge(wave.1, wave.4, by=c('aid', 'psuscid', 'region'), all=TRUE)
final.df <- merge(final.df, wave.5, by=c('aid', 'psuscid', 'region'), all=TRUE)
final.df <- merge(final.df, ses.data, by='aid')

# Clear environment
rm(list=setdiff(ls(), 'final.df'))
################################################################################
# CREATE VARIABLES THAT WILL BE USED IN ANALYSIS:
################################################################################
# GE Change Variables

final.df <- final.df %>%
  mutate(
    above_school_avg = case_when(
      w1.GE_male >= school_avg_GE ~ 1,
      w1.GE_male < school_avg_GE ~ 0,
      TRUE ~ NaN
    ),
    delta_w1_w4_GE = w4.GE_male_std - w1.GE_male_std,
    increasing = ifelse(delta_w1_w4_GE > 0, 1, 0),
    school_self = case_when(
      above_school_avg == 1 & delta_w1_w4_GE > 0 ~ 1,
      above_school_avg == 1 & delta_w1_w4_GE < 0 ~ 2,
      above_school_avg == 0 & delta_w1_w4_GE > 0 ~ 3,
      above_school_avg == 0 & delta_w1_w4_GE < 0 ~ 4
    ),
    self_self = case_when(
      w1.GE_male_std >= 0 & delta_w1_w4_GE > 0 ~ 1,
      w1.GE_male_std >= 0 & delta_w1_w4_GE < 0 ~ 2,
      w1.GE_male_std < 0 & delta_w1_w4_GE > 0 ~ 3,
      w1.GE_male_std < 0 & delta_w1_w4_GE < 0 ~ 4
    )
  )

final.df$self_self <- factor(
  final.df$self_self,
  levels = c(1, 2, 3, 4),
  labels = c(
    "Above Adolescent Male Avg and Increasing",
    "Above Adolescent Male Avg and Decreasing",
    "Below Adolescent Male Avg and Increasing",
    "Below Adolescent Male Avg and Decreasing"
  )
)

################################################################################
# Determine what variables to keep for analysis

vars.keep <- c('dx_cad5', 'dx_htn5', 'dx_dm5', 'dx_hld5', 
               'w5_male', 'dx_bmi52', 'sr_bmi', 'gsw5',
               'w5_anti_hld_med_use', 'w5_anti_dm_med_use','w5_anti_htn', 
               'aid', 'region', 'psuscid', 'gswgt4_2', 'sschlcde', 
               'w1_male', 'self_self', 'delta_w1_w4_GE', 
               'w1.GE_male','w4.GE_male', 
               'w1.GE_male_std', 'w4.GE_male_std',
               'sespc_al', 'nhood1_d', 'race', 
               'w4_male', 'above_school_avg', 'w5_nhdl',
               'lipid5cat', 'w5_a1c', 'w5_bp', 'race5_d', 'race5', 'ins5',
               'edu5', 'w5biowgt')


final.df <- final.df[, vars.keep]

# Create a binary variable indicating whether a respondent is in the sample
final.df <- final.df %>%
  mutate(in_sample = ifelse(w1_male == w4_male & w1_male == w5_male & w1_male == 1, 1, 0),
         in_sample = ifelse(is.na(w1.GE_male) == F & is.na(w4.GE_male) == F, in_sample, 0),
         in_sample.5  = ifelse(in_sample == 1 & is.na(gsw5) == F, 1, 0),
         in_sample.bio = ifelse(in_sample == 1 & is.na(w5biowgt) == F, 1, 0))


#create cluster var
final.df$cluster <- paste(final.df$region,final.df$psuscid)
final.df$weights <- final.df$gsw5 / mean(final.df$gsw5, na.rm = T)
final.df$weights.bio <- final.df$w5biowgt / mean(final.df$w5biowgt, na.rm = T)

################################################################################