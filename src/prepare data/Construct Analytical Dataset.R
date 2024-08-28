################################################################################
# AUTHOR:             Jacob Jameson
#
# DESCRIPTION:        This script constructs the analytical dataset that will be 
#                     used for the analysis. The script merges the wave 1, 
#                     wave 4, wave 5, and SES data together.
#
# DEPENDENCIES:       Prepare W1.R, Prepare W4.R, Prepare W5.R
################################################################################

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
final.df <- merge(final.df, ses.data, by='aid', all=TRUE)

# Clear environment
rm(list=setdiff(ls(), 'final.df'))
################################################################################
# CREATE VARIABLES THAT WILL BE USED IN ANALYSIS:
################################################################################
# GE Change Variables

final.df <- final.df %>%
  mutate(
    delta_w1_w4_GE = w4.GE_male_std - w1.GE_male_std,
    increasing = ifelse(delta_w1_w4_GE > 0, 1, 0),
    self_self = case_when(
      w1.GE_male_std >= 0 & delta_w1_w4_GE > 0 ~ 1,
      w1.GE_male_std >= 0 & delta_w1_w4_GE < 0 ~ 2,
      w1.GE_male_std < 0 & delta_w1_w4_GE > 0 ~ 3,
      w1.GE_male_std < 0 & delta_w1_w4_GE < 0 ~ 4
    )
  )

################################################################################
# determine ages of respondents at each wave

# wave 1 ---------------------
final.df$w1_age <- final.df$iyear - final.df$h1gi1y + (final.df$imonth - final.df$h1gi1m) / 12

# wave 4 ---------------------
final.df$w4_age <- final.df$iyear4 - (1900 + final.df$h1gi1y) + (final.df$imonth4 - final.df$h1gi1m) / 12

# wave 5 ---------------------
final.df$w5_age <- final.df$iyear5 - (1900 + final.df$h1gi1y) + (final.df$imonth5 - final.df$h1gi1m) / 12

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
               'w4_male', 'w5_nhdl',
               'lipid5cat', 'w5_a1c', 'w5_bp', 'race5_d', 'race5', 'ins5',
               'edu5', 'w5biowgt', 'w1_age', 'w4_age', 'w5_age',
               'h1gh21', 'h1da5', 'h1fv5', 'h1gh28', 'h1pr4', 'h1da10', 'h1gh46', 
               'h1ee4', 'h1ed7', 'h1fs2', 'h1gh39', 'h1da11', 'h1da1', 'h1pf15', 
               'h1pr1', 'h1gh20', 'h1pf32', 'h1id5', 'h1da6', 'h1pf16', 'h1gh29', 
               'h1pf10', 'h1ee2', 'h1fs4', 'h1gh42',
               'h4to25', 'h4cj1', 'h4da16', 'h4pe5', 'h4pe9', 'h4pe2', 'h4da6', 
               'h4da23', 'h4da8', 'h4pe4', 'h4re10', 'h4da17', 'h4mi1', 'h4da4', 
               'h4mh23', 'h4pe6', 'h4mh7', 'h4pe10', 'h4pe35', 'h4da11', 'h4pe22', 
               'h4pe26')


final.df <- final.df[, vars.keep]

# Create a binary variable indicating whether a respondent is in the sample
final.df <- final.df %>%
  mutate(in_sample0 = ifelse(w1_male == 1 & w4_male == 1 & w5_male == 1, 1, 0),
         in_sample = ifelse(w1_male == 1 & w4_male == 1 & w5_male == 1, 1, 0),
         in_sample = ifelse(is.na(w1.GE_male) == F & is.na(w4.GE_male) == F, in_sample, 0),
         in_sample.5  = ifelse(in_sample == 1 & is.na(gsw5) == F, 1, 0),
         in_sample.bio = ifelse(in_sample == 1 & is.na(w5biowgt) == F, 1, 0))


#create cluster var
final.df$cluster <- paste(final.df$region,final.df$psuscid)
final.df$weights <- final.df$gsw5 / mean(final.df$gsw5, na.rm = T)
final.df$weights.bio <- final.df$w5biowgt / mean(final.df$w5biowgt, na.rm = T)

final.df$tx_dx_bp <- ifelse(
  final.df$dx_htn5 == 1 & final.df$w5_anti_htn == 1, 1, 0)

final.df$tx_dx_hld <- ifelse(
  final.df$dx_hld5 == 1 & final.df$w5_anti_hld_med_use == 1, 1, 0)

final.df$tx_dx_dm <- ifelse(
  final.df$dx_dm5 == 1 & final.df$w5_anti_dm_med_use == 1, 1, 0)

final.df$dx_bio_bp <- ifelse(
  final.df$dx_htn5 == 1 & final.df$w5_bp == 1, 1, 0)

final.df$dx_bio_hld <- ifelse(
  final.df$dx_hld5 == 1 & final.df$w5_nhdl == 1, 1, 0)

final.df$dx_bio_dm <- ifelse(
  final.df$dx_dm5 == 1 & final.df$w5_a1c == 1, 1, 0)

################################################################################