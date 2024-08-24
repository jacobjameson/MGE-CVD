################################################################################
# AUTHOR:             Jacob Jameson
#
# DESCRIPTION:        This script constructs the analytical dataset that will be
#                     used for the analysis. The script merges the wave 5 data
#                     with the biomarker data.
################################################################################

# Data paths ------------------------------------------------------------
data_path <- '~/Sue Goldie Dropbox/Jacob Jameson/Add Health/Data Update 8.2021/Core Files - Wave V'
survey_path <-  paste0(data_path, '/Wave V Mixed-Mode Survey Data')
weights_path <-  paste0(data_path, '/Wave V Mixed-Mode Survey Weights')

# Load the wave 5 data and wave 5 weights
wave.5 <- read_xpt(paste0(survey_path, '/wave5.xpt'))
weights.5 <- read_xpt(paste0(weights_path, '/weights5.xpt'))

# Merge wave 5 data with the weights
wave.5 <- merge(wave.5, weights.5, by='AID', all = T)

# Remove data no longer using
rm(weights.5)
################################################################################
### Incorporate biomarker data
bio_path <- '~/Sue Goldie Dropbox/Jacob Jameson/Add Health/Data Upload 7.2021/Wave V Biomarker Files'

# Data paths
cardio_path <-  paste0(bio_path, '/Wave V Cardiovascular Measures')
lipids_path <-  paste0(bio_path, '/Wave V Lipids')
glucose_path <-  paste0(bio_path, '/Wave V Glucose Homeostasis')
meds_path <-  paste0(bio_path, '/Wave V Medications - Home Exam')
weights_path <-  paste0(bio_path, '/Wave V Biomarker Sample Weight')

# Load Data
bmeds <- read_xpt(paste0(meds_path, '/bmeds5.xpt'))
cardio <- read_xpt(paste0(cardio_path, '/bcardio5.xpt'))
lipids <- read_xpt(paste0(lipids_path, '/blipids5.xpt'))
glucose <- read_xpt(paste0(glucose_path, '/bglua1c5.xpt'))
weights <- read_xpt(paste0(weights_path, '/bweight5.xpt'))

# Merge Data
wave.5 <- merge(wave.5, cardio, by='AID', all = T)
wave.5 <- merge(wave.5, lipids, by='AID', all = T)
wave.5 <- merge(wave.5, glucose, by='AID', all = T)
wave.5 <- merge(wave.5, weights, by='AID', all = T)


names(wave.5) <- tolower(names(wave.5))

# BMI Construction
wave.5 <- wave.5 %>%
  mutate(height_fti = ifelse(h5id2f != 98, h5id2f*12, NA),
         height_in = ifelse(h5id2i != 998, h5id2i, NA),
         height_sq = (height_fti + height_in)**2,
         weight = h5id3,
         sr_bmi = weight/height_sq*703,
         dx_bmi52 = case_when(
           sr_bmi >= 0 & sr_bmi < 30 ~ 0,
           sr_bmi >= 30 & sr_bmi < 100 ~ 1))

# Diagnoses
wave.5 <- wave.5 %>%
  mutate(dx_cad5 = h5id6e,
         dx_htn5 = h5id6c,
         dx_dm5 = h5id6d,
         dx_hld5 = h5id6b,
         w5_male = ifelse(h5od2b == 1, 1, 0))

# Medications
wave.5 <- wave.5 %>%
  mutate(
    w5_anti_htn = case_when(
      h5aht == 0 ~ 0,
      h5aht == 1 ~ 1),
    w5_anti_dm_med_use = case_when(
      h5c_med == 0 ~ 0,
      h5c_med == 1 ~ 1),
    w5_anti_hld_med_use = case_when(
      h5c_med2 == 0 ~ 0,
      h5c_med2 == 1 ~ 1))

# Biomarker data
wave.5 <- wave.5 %>%
  mutate(
    lipid5cat = case_when(
      h5tc >= 0 & h5tc < 200 ~ 0,
      h5tc >= 200 & h5tc < 1000 ~ 1),
    w5_a1c = case_when(
      h5chba1c == 1 | h5chba1c == 2  ~ 0,
      h5chba1c == 3 ~ 1),
    w5_bp = case_when(
      h5bpcls5 == 1 | h5bpcls5 == 2  ~ 0,
    #  h5bpcls5 == 1 | h5bpcls5 == 2 | h5bpcls5 == 3  ~ 0,
      h5bpcls5 == 3 | h5bpcls5 == 4 | h5bpcls5 == 5 ~ 1),
      #h5bpcls5 == 4 | h5bpcls5 == 5 ~ 1),
    w5_nhdl = case_when(
     h5nhdl >= 0 & h5nhdl < 190 ~ 0,
     h5nhdl >= 190 ~ 1 ))


wave.5 <- wave.5 %>%
  mutate(edu5 = case_when(
    h5od11 == 1 ~ 1,
    h5od11 == 2 ~ 1,
    h5od11 == 3 ~ 2,
    h5od11 == 4 ~ 2,
    h5od11 == 5 ~ 3,
    h5od11 == 6 ~ 3,
    h5od11 == 7 ~ 3,
    h5od11 == 8 ~ 3,
    h5od11 == 9 ~ 3,
    h5od11 == 10 ~ 4,
    h5od11 >= 11 & h5od11 <= 16 ~ 4),
    edu5 = factor(edu5, levels = c(1,2,3,4),
                labels = c("Some HS or Less", "HS Diploma/GED",
                           "Some College or Tech/Assoc Degree",
                           "College Degree or More")))

wave.5 <- wave.5 %>%
  mutate(ins5 = case_when(
    h5id9 == 15 ~ 1,
    h5id9 >= 1 & h5id9 <= 5 ~ 2,
    h5id9 == 7 ~ 2,
    h5id9 == 8 ~ 2,
    h5id9 == 9 ~ 3,
    h5id9 == 10 ~ 3,
    h5id9 == 6 ~ 4,
    h5id9 >= 11 & h5id9 < 14 ~ 4),
    ins5 = factor(ins5, levels = c(1,2,3,4),
                       labels = c("Uninsured", "Private or Employer-Based",
                                  "Medicaid/Medicare", "Other Government")))

wave.5 <- wave.5 %>%
  mutate(race5 = case_when(
    h5od8 == 1 ~ 1,
    h5od8 == 2 ~ 2,
    h5od8 >= 3 &  h5od8 <= 7 ~ 3,
    h5od8 >= 8 & h5od8 <= 20 ~ 4,
    h5od8 == 21 | h5od8 == 22 ~ 5,
    h5od4a == 1 ~ 1,
    h5od4b == 1 ~ 2,
    h5od4c == 1 ~ 3,
    h5od4d == 1 ~ 4,
    h5od4e == 1 ~ 4,
    h5od4f == 1 ~ 5,
    h5od4g == 1 ~ 5),
    race5 = factor(race5, levels = c(1,2,3,4,5),
                        labels = c("White", "Black/Af-Am", "Hispanic",
                                   "Asian/Pacific Islander", "Other")),
    race5_d = ifelse(race5 == 3 | race5 == 4 | race5 == 5, 3, race5),
    race5_d = factor(race5_d, levels = c(1,2,3),
                   labels = c("White", "Black/Af-Am", "Other")))


wave.5$w5_male <- factor(ifelse(wave.5$h5od2b == 1, 1, 0))
################################################################################