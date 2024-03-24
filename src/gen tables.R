#-------------------------------------------------------------------------------
# AUTHOR:             Jacob Jameson
# PURPOSE:            Generate main regression results
#-------------------------------------------------------------------------------

rm(list = ls())

# load packages ---------------------------------------------------------------

library(lmtest)
library(sandwich)
library(broom)
library(tidyverse)
library(sjstats)
library(MASS)
library(stargazer)
library(gtsummary)
library(stargazer)
library(marginaleffects)

# load data --------------------------------------------------------------------
source('src/prepare data/Construct Analytical Dataset.R')

#-------------------------------------------------------------------------------
################################################################################
# Table 2. Marginal effects coefficients (dy/dx) estimating associations 
# between adolescent, young adult, and adolescent-to-young-adult change 
# in male gender expression (MGE) and adult awareness of hypertension, 
# diabetes, and hyperlipidemia 
################################################################################
#-------------------------------------------------------------------------------

# function to generate marginal effects coefficients (dy/dx) estimating
ape_analysis <- function(frml, df, weights){
  
  frml1 <- as.formula(frml)
  model <- glm(formula = frml1, data = df, 
            weights = weights, family = "quasibinomial")
  
  vcv_robust <- vcovHC(model, type = "HC1", cluster = ~ cluster)
  
  ape <- avg_slopes(model = model, vcov = vcv_robust) %>% 
    as.data.frame() %>%
    mutate(outcome = str_squish(word(frml,1,sep = "\\~"))) %>%
    filter(term == str_squish(word(frml,3,sep = " "))) 
  return(ape)
}

# function to generate marginal effects coefficients (dy/dx) by group
ape_analysis_by <- function(frml, df, weights, by){
  
  frml1 <- as.formula(frml)
  model <- glm(formula = frml1, data = df, 
            weights = weights, family = "quasibinomial")
  
  vcv_robust <- vcovHC(model, type = "HC1", cluster = ~ cluster)
  
  ape <- avg_slopes(model = model, vcov = vcv_robust, 
                    by = by) %>%
    as.data.frame() %>%
    mutate(outcome = str_squish(word(frml,1,sep = "\\~"))) %>%
    filter( term == str_squish(word(frml,3,sep = " "))) 

  return(ape)
}

## Model 1 --------------------------------------------------------------------

# Define the subset and weights
subset_df <- subset(final.df, in_sample.5 == 1)
subset_weights <- subset_df$weights

### Formulas

formulas_htn <- c(
  'dx_htn5 ~ w1.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'dx_htn5 ~ w4.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'dx_htn5 ~ delta_w1_w4_GE + w1.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5'
)


formulas_dm <- c(
  'dx_dm5 ~ w1.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'dx_dm5 ~ w4.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'dx_dm5 ~ delta_w1_w4_GE + w1.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5'
)

formulas_hld <- c(
  'dx_hld5 ~ w1.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'dx_hld5 ~ w4.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'dx_hld5 ~ delta_w1_w4_GE + w1.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5'
)

# Perform analysis 
results_htn <- map_dfr(
  formulas_htn, ~ ape_analysis(.x, subset_df, subset_weights)
  ) %>%
  mutate(outcome = 'Dx Hypertension')

results_dm <- map_dfr(
  formulas_dm, ~ ape_analysis(.x, subset_df, subset_weights)
) %>%
  mutate(outcome = 'Dx Diabetes')

results_hld <- map_dfr(
  formulas_hld, ~ ape_analysis(.x, subset_df, subset_weights)
) %>%
  mutate(outcome = 'Dx Hyperlipidemia')

# append data frames
results <- bind_rows(results_hld, results_htn, results_dm)

model.1 <- results %>%
  mutate(term = case_when(
    grepl("w1.GE_male_std", term) ~ "Adolescent MGE",
    grepl("w4.GE_male_std", term) ~ "Young Adult MGE",
    TRUE ~ "Change in MGE"
  ),
  outcome = factor(outcome, levels = c("Dx Hypertension", "Dx Diabetes", "Dx Hyperlipidemia")),
  term = factor(term, levels = c("Adolescent MGE", "Young Adult MGE", "Change in MGE"))
  ) %>%
  group_by(outcome, term) %>%
  summarize(
    Estimate = formatC(mean(estimate), format = "f", digits = 3),
    std.error = formatC(mean(std.error), format = "f", digits = 3),
    CI = paste0("(", formatC(mean(conf.low), format = "f", digits = 3), 
                ", ", formatC(mean(conf.high), format = "f", digits = 3), ")"),
    P_Value = formatC(mean(p.value), format = "f", digits = 3)) %>%
  ungroup()


## Model 2 --------------------------------------------------------------------

# Define the subset and weights
subset_df <- subset(final.df, in_sample.bio == 1)
subset_weights <- subset_df$weights


formulas_htn <- c(
  'dx_htn5 ~ w1.GE_male_std + w5_bp + w1.GE_male_std*w5_bp + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'dx_htn5 ~ w4.GE_male_std + w5_bp + w4.GE_male_std*w5_bp + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'dx_htn5 ~ delta_w1_w4_GE + w5_bp + delta_w1_w4_GE*w5_bp + w1.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5'
)

formulas_dm <- c(
  'dx_dm5 ~ w1.GE_male_std + w5_a1c + w1.GE_male_std*w5_a1c + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'dx_dm5 ~ w4.GE_male_std + w5_a1c + w4.GE_male_std*w5_a1c + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'dx_dm5 ~ delta_w1_w4_GE + w5_a1c + delta_w1_w4_GE*w5_a1c + w1.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5'
)

formulas_hld <- c(
  'dx_hld5 ~ w1.GE_male_std + w5_nhdl + w1.GE_male_std*w5_nhdl + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'dx_hld5 ~ w4.GE_male_std + w5_nhdl + w4.GE_male_std*w5_nhdl + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'dx_hld5 ~ delta_w1_w4_GE + w5_nhdl + delta_w1_w4_GE*w5_nhdl + w1.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5'
)


# Perform analysis 
results_htn <- map_dfr(
  formulas_htn, ~ ape_analysis_by(.x, subset_df, subset_weights, 'w5_bp')
) %>%
  mutate(outcome = 'Dx Hypertension') %>%
  filter(w5_bp == 1)

results_dm <- map_dfr(
  formulas_dm, ~ ape_analysis_by(.x, subset_df, subset_weights, 'w5_a1c')
) %>%
  mutate(outcome = 'Dx Diabetes') %>%
  filter(w5_a1c == 1)

results_hld <- map_dfr(
  formulas_hld, ~ ape_analysis_by(.x, subset_df, subset_weights, 'w5_nhdl')
) %>%
  mutate(outcome = 'Dx Hyperlipidemia') %>%
  filter(w5_nhdl == 1)

# append data frames
results <- bind_rows(results_hld, results_htn, results_dm)

model.2 <- results %>%
  mutate(term = case_when(
    grepl("w1.GE_male_std", term) ~ "Adolescent MGE",
    grepl("w4.GE_male_std", term) ~ "Young Adult MGE",
    TRUE ~ "Change in MGE"
  ),
  outcome = factor(outcome, levels = c("Dx Hypertension", "Dx Diabetes", "Dx Hyperlipidemia")),
  term = factor(term, levels = c("Adolescent MGE", "Young Adult MGE", "Change in MGE"))
  ) %>%
  group_by(outcome, term) %>%
  summarize(
    Estimate = formatC(mean(estimate), format = "f", digits = 3),
    std.error = formatC(mean(std.error), format = "f", digits = 3),
    CI = paste0("(", formatC(mean(conf.low), format = "f", digits = 3), 
                ", ", formatC(mean(conf.high), format = "f", digits = 3), ")"),
    P_Value = formatC(mean(p.value), format = "f", digits = 3)) %>%
  ungroup()

## Model 3 --------------------------------------------------------------------

# Define the subset and weights
subset_df <- subset(final.df, in_sample.5 == 1)
subset_weights <- subset_df$weights

formulas_htn <- c(
  'w5_anti_htn ~ w1.GE_male_std + w1.GE_male_std*dx_htn5 + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'w5_anti_htn ~ w4.GE_male_std + w4.GE_male_std*dx_htn5 + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'w5_anti_htn ~ delta_w1_w4_GE + delta_w1_w4_GE*dx_htn5 + w1.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5'
)

formulas_dm <- c(
  'w5_anti_dm_med_use ~ w1.GE_male_std + w1.GE_male_std*dx_dm5 + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'w5_anti_dm_med_use ~ w4.GE_male_std + w4.GE_male_std*dx_dm5 + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'w5_anti_dm_med_use ~ delta_w1_w4_GE + delta_w1_w4_GE*dx_dm5 + w1.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5'
)

formulas_hld <- c(
  'w5_anti_hld_med_use ~ w1.GE_male_std + w1.GE_male_std*dx_hld5 + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'w5_anti_hld_med_use ~ w4.GE_male_std + w4.GE_male_std*dx_hld5 + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'w5_anti_hld_med_use ~ delta_w1_w4_GE + delta_w1_w4_GE*dx_hld5 + w1.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5'
)

# Perform analysis 
results_htn <- map_dfr(
  formulas_htn, ~ ape_analysis_by(.x, subset_df, subset_weights, 'dx_htn5')
) %>%
  mutate(outcome = 'Anti Hypertension Med') %>%
  filter(dx_htn5 == 1)

results_dm <- map_dfr(
  formulas_dm, ~ ape_analysis_by(.x, subset_df, subset_weights, 'dx_dm5')
) %>%
  mutate(outcome = 'Anti Diabetes Med') %>%
  filter(dx_dm5 == 1)

results_hld <- map_dfr(
  formulas_hld, ~ ape_analysis_by(.x, subset_df, subset_weights, 'dx_hld5')
) %>%
  mutate(outcome = 'Anti Hyperlipidemia Med') %>%
  filter(dx_hld5 == 1)

# append data frames
results <- bind_rows(results_hld, results_htn, results_dm)

model.3 <- results %>%
  mutate(term = case_when(
    grepl("w1.GE_male_std", term) ~ "Adolescent MGE",
    grepl("w4.GE_male_std", term) ~ "Young Adult MGE",
    TRUE ~ "Change in MGE"
  ),
  outcome = factor(outcome, 
                   levels = c("Anti Hypertension Med", 
                              "Anti Diabetes Med",
                              "Anti Hyperlipidemia Med")),
  term = factor(term, levels = c("Adolescent MGE", "Young Adult MGE", "Change in MGE"))
  ) %>%
  group_by(outcome, term) %>%
  summarize(
    Estimate = formatC(mean(estimate), format = "f", digits = 3),
    std.error = formatC(mean(std.error), format = "f", digits = 3),
    CI = paste0("(", formatC(mean(conf.low), format = "f", digits = 3), 
                ", ", formatC(mean(conf.high), format = "f", digits = 3), ")"),
    P_Value = formatC(mean(p.value), format = "f", digits = 3)) %>%
  ungroup()


## Model 4 --------------------------------------------------------------------

# Define the subset and weights
subset_df <- subset(final.df, in_sample.bio == 1)
subset_weights <- subset_df$weights

formulas_htn <- c(
  'w5_bp ~ w1.GE_male_std + w5_anti_htn + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'w5_bp ~ w4.GE_male_std + w5_anti_htn + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'w5_bp ~ delta_w1_w4_GE + w1.GE_male_std + w5_anti_htn + race5 + sespc_al + nhood1_d + ins5 + edu5'
)

formulas_dm <- c(
  'w5_a1c ~ w1.GE_male_std + w5_anti_dm_med_use + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'w5_a1c ~ w4.GE_male_std + w5_anti_dm_med_use + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'w5_a1c ~ delta_w1_w4_GE + w1.GE_male_std + w5_anti_dm_med_use + race5 + sespc_al + nhood1_d + ins5 + edu5'
)

formulas_hld <- c(
  'w5_nhdl ~ w1.GE_male_std + w5_anti_hld_med_use + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'w5_nhdl ~ w4.GE_male_std + w5_anti_hld_med_use + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'w5_nhdl ~ delta_w1_w4_GE + w1.GE_male_std + w5_anti_hld_med_use + race5 + sespc_al + nhood1_d + ins5 + edu5'
)

# Perform analysis 
results_htn <- map_dfr(
  formulas_htn, ~ ape_analysis(.x, subset_df, subset_weights)
) %>%
  mutate(outcome = 'Biomarker Hypertension')

results_dm <- map_dfr(
  formulas_dm, ~ ape_analysis(.x, subset_df, subset_weights)
) %>%
  mutate(outcome = 'Biomarker Diabetes')

results_hld <- map_dfr(
  formulas_hld, ~ ape_analysis(.x, subset_df, subset_weights)
) %>%
  mutate(outcome = 'Biomarker Hyperlipidemia')

# append data frames
results <- bind_rows(results_hld, results_htn, results_dm)

model.4 <- results %>%
  mutate(term = case_when(
    grepl("w1.GE_male_std", term) ~ "Adolescent MGE",
    grepl("w4.GE_male_std", term) ~ "Young Adult MGE",
    TRUE ~ "Change in MGE"
  ),
  outcome = factor(outcome, levels = c("Biomarker Hypertension", 
                                       "Biomarker Diabetes", 
                                       "Biomarker Hyperlipidemia")),
  term = factor(term, levels = c("Adolescent MGE", "Young Adult MGE", "Change in MGE"))
  ) %>%
  group_by(outcome, term) %>%
  summarize(
    Estimate = formatC(mean(estimate), format = "f", digits = 3),
    std.error = formatC(mean(std.error), format = "f", digits = 3),
    CI = paste0("(", formatC(mean(conf.low), format = "f", digits = 3), 
                ", ", formatC(mean(conf.high), format = "f", digits = 3), ")"),
    P_Value = formatC(mean(p.value), format = "f", digits = 3)) %>%
  ungroup()

################################################################################
# in each model add stars to estimates depending on significance level
################################################################################
model.1 <- model.1 %>%
  mutate(
    stars = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01 ~ "**",
      P_Value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )
model.1$Estimate <- paste0(model.1$Estimate, model.1$stars)

model.2 <- model.2 %>%
  mutate(
    stars = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01 ~ "**",
      P_Value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )
model.2$Estimate <- paste0(model.2$Estimate, model.2$stars)

model.3 <- model.3 %>%
  mutate(
    stars = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01 ~ "**",
      P_Value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )
model.3$Estimate <- paste0(model.3$Estimate, model.3$stars)

model.4 <- model.4 %>%
  mutate(
    stars = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01 ~ "**",
      P_Value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )
model.4$Estimate <- paste0(model.4$Estimate, model.4$stars)


# Start the sink to capture output
sink("outputs/tables/Table 2.txt")

print('Self-Report of Adult Diagnosis (Model 1)')
print('-----------------------------------------------------------------------')

print(model.1)
print('-----------------------------------------------------------------------')
print('-----------------------------------------------------------------------')
print('-----------------------------------------------------------------------')

print('Self-Report of Adult Diagnosis Among men with elevated Biomeasure (Model 2)')
print('-----------------------------------------------------------------------')

print(model.2)

print('-----------------------------------------------------------------------')
print('-----------------------------------------------------------------------')
print('-----------------------------------------------------------------------')

print('Biomeasure Evidence of CVD Risks (Model 3)')
print('-----------------------------------------------------------------------')

print(model.3)

print('-----------------------------------------------------------------------')
print('-----------------------------------------------------------------------')
print('-----------------------------------------------------------------------')


print('Treatment of CVD Risks (Model 4)')
print('-----------------------------------------------------------------------')

print(model.4)

print('-----------------------------------------------------------------------')
print('-----------------------------------------------------------------------')
print('-----------------------------------------------------------------------')

sink()
