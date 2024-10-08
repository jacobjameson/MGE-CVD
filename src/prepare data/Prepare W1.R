################################################################################
# AUTHOR:             Jacob Jameson
#
# DESCRIPTION:        This script constructs the analytical dataset that will be
#                     used for the analysis. The script merges the wave 1 data
#                     with the wave 1 weights. The script also creates the variables
#                     that will be used in the analysis.
################################################################################

# Data paths ------------------------------------------------------------
data_path <- '~/Sue Goldie Dropbox/Jacob Jameson/Add Health/Data Upload 7.2021/Core Files - Wave I'

inhome_path <-  paste0(data_path, '/Wave I In Home Interview Data')
weights_path <-  paste0(data_path, '/Wave I In-Home Weights')

# Load the wave 1 data and wave 1 weights -------------------------------
allwave.1 <- read_xpt(paste0(inhome_path, '/allwave1.xpt'))
homeweights.1 <- read_xpt(paste0(weights_path, '/Homewt1.xpt'))

# Merge wave 1 data with the weights ------------------------------------
wave.1 <- merge(allwave.1, homeweights.1, by='AID')

# Rename variables lowercase --------------------------------------------
names(wave.1) <- tolower(names(wave.1))

# Remove data no longer using -------------------------------------------
rm(allwave.1, homeweights.1)

#-------------------------------------------------------------------------
# Create variables that will be used in the analysis
#
#   - w1.GE_male: GE score created based on Flemming et al.
#   - w1.GE_male_std: GE score by males standardized
#   - school_avg_GE: Avg GE score by (school, sex)
#-------------------------------------------------------------------------
#
# Construct pseudo.gpa variable and race ---------------------------------

wave.1 <- wave.1 %>% 
  mutate(race = case_when(
    h1gi4 == 1 ~ "Hispanic",
    h1gi6a == 1 & (h1gi6b != 1 & h1gi6c != 1 & h1gi6d != 1 & h1gi6e != 1) ~ "Non Hispanic White",
    h1gi6b == 1 & (h1gi6a != 1 & h1gi6c != 1 & h1gi6d != 1 & h1gi6e != 1) ~ "Non Hispanic Black/Af-Am",
    h1gi6c == 1 & (h1gi6a != 1 & h1gi6b != 1 & h1gi6d != 1 & h1gi6e != 1) ~ "Multiple Races or Other",
    h1gi6d == 1 & (h1gi6a != 1 & h1gi6b != 1 & h1gi6c != 1 & h1gi6e != 1)  ~ "Non Hispanic Asian/Pacific Islander",
    h1gi6e == 1 ~ "Multiple Races or Other",
    TRUE ~ "Multiple Races or Other"
  ))

wave.1$race <- factor(wave.1$race,
                      levels = c("Non Hispanic White", "Non Hispanic Black/Af-Am",
                                 "Non Hispanic Asian/Pacific Islander",
                                 "Hispanic",
                                 "Multiple Races or Other"),
                      labels = c("Non Hispanic White", "Non Hispanic Black/Af-Am",
                                 "Non Hispanic Asian/Pacific Islander",
                                 "Hispanic",
                                 "Multiple Races or Other"))

# Clean gender expression variables for regression model ----------------------
wave.1 <- wave.1 %>%
  mutate(h1gh21 = ifelse(h1gh21 > 4 | h1gh21 < 0, NA, h1gh21),
         h1da5 = ifelse(h1da5 > 3 | h1da5 < 0, NA, h1da5),
         h1fv5 = ifelse(h1fv5 > 2 | h1fv5 < 0, NA, h1fv5),
         h1gh28 = ifelse(h1gh28 > 5 | h1gh28 < 1, NA, h1gh28),
         h1pr4 = ifelse(h1pr4 > 5 | h1pr4 < 1, NA, h1pr4),
         h1da10 = ifelse(h1da10 > 99 | h1da10 < 0, NA, h1da10),
         h1gh46 = ifelse(h1gh46 > 5 | h1gh46 < 1, NA, h1gh46),
         h1ee4 = ifelse(h1ee4 > 145 | h1ee4 < 0, NA, h1ee4),
         h1ed7 = ifelse(h1ed7 > 1 | h1ed7 < 0, NA, h1ed7),
         h1fs2 = ifelse(h1fs2 > 3 | h1fs2 < 0, NA, h1fs2),
         h1gh39 = ifelse(h1gh39 > 5 | h1gh39 < 0, NA, h1gh39),
         h1da11 = ifelse(h1da11 > 99 | h1da11 < 0, NA, h1da11),
         h1da1 = ifelse(h1da1 > 3 | h1da1 < 0, NA, h1da1),
         h1pf15 = ifelse(h1pf15 > 5 | h1pf15 < 1, NA, h1pf15),
         h1pr1 = ifelse(h1pr1 > 5 | h1pr1 < 1, NA, h1pr1),
         h1gh20 = ifelse(h1gh20 > 4 | h1gh20 < 0, NA, h1gh20),
         h1pf32 = ifelse(h1pf32 > 5 | h1pf32 < 1, NA, h1pf32),
         h1id5 = ifelse(h1id5 > 1 | h1id5 < 0, NA, h1id5),
         h1da6 = ifelse(h1da6 > 3 | h1da6 < 0, NA, h1da6),
         h1pf16 = ifelse(h1pf16 > 5 | h1pf16 < 1, NA, h1pf16),
         h1gh29 = ifelse(h1gh29 > 4 | h1gh29 < 1, NA, h1gh29),
         h1pf10 = ifelse(h1pf10 > 5 | h1pf10 < 1, NA, h1pf10),
         h1ee2 = ifelse(h1ee2 > 5 | h1ee2 < 1, NA, h1ee2),
         h1fs4 = ifelse(h1fs4 > 3 | h1fs4 < 0, NA, h1fs4),
         h1gh42 = ifelse(h1gh42 > 4 | h1gh42 < 0, NA, h1gh42))
  
vars <- c("h1gh21", "h1da5", "h1fv5", "h1gh28", "h1pr4", "h1da10", "h1gh46",
          "h1ee4", "h1ed7", "h1fs2", "h1gh39", "h1da11", "h1da1", "h1pf15",
          "h1pr1", "h1gh20", "h1pf32", "h1id5", "h1da6", "h1pf16", "h1gh29",
          "h1pf10", "h1ee2", "h1fs4", "h1gh42")

# group_by bio_sex and impute missing values with median ---------------------
#wave.1 <- wave.1 %>%
#  group_by(bio_sex) %>%
#  mutate_at(vars, funs(replace(., is.na(.), median(., na.rm = TRUE))) ) %>%
#  ungroup()

# Logistic regression to predict male using GE variables ------------------------
wave.1$w1_male <- factor(ifelse(wave.1$bio_sex == 1, 1, 0))
wave.1$w1_female <- factor(ifelse(wave.1$bio_sex == 2, 1, 0))

predict_male <- glm(w1_male ~ factor(h1gh21) + factor(h1da5) + factor(h1fv5) +
                      factor(h1gh28) + factor(h1pr4) + h1da10 + factor(h1gh46) +
                      h1ee4 + factor(h1ed7) + factor(h1fs2) + factor(h1gh39) +
                      h1da11 + factor(h1da1) + factor(h1pf15) + factor(h1pr1) +
                      factor(h1gh20) + factor(h1pf32) + factor(h1id5) + 
                      factor(h1da6) + factor(h1pf16) + factor(h1gh29) + 
                      factor(h1pf10) + factor(h1ee2) + factor(h1fs4) + 
                      factor(h1gh42), data = wave.1, family = "binomial")

# Get a prediction score for male. This will be our GE score ------------------
wave.1$w1.GE_male <- predict(predict_male, wave.1, type="response")

wave.1 <- wave.1 %>%
  group_by(w1_male) %>%
  mutate(w1.GE_male_std = scale(w1.GE_male) %>% as.vector) %>% 
  ungroup()

wave.1 <- wave.1 %>%
  group_by(w1_male, sschlcde) %>%
  mutate(w1.GE_male_std_school = scale(w1.GE_male) %>% as.vector) %>% 
  ungroup()

# Create a score for school average gender expression grouped by sex ----------
wave.1 <- wave.1 %>%
  group_by(sschlcde, w1_male) %>%
  mutate(school_avg_GE = mean(w1.GE_male, na.rm=TRUE)) %>% ungroup()

# ------------------------------------------------------------------------------