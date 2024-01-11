#-------------------------------------------------------------------------------
# AUTHOR:             Jacob Jameson
# PURPOSE:            Generate main figures
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

# function to generate marginal effects coefficients (dy/dx) estimating
ape_analysis <- function(frml, df, weights){
  
  frml1 <- as.formula(frml)
  model <- glm(formula = frml1, data = df, 
               weights = weights, family = "quasibinomial")
  
  vcv_robust <- vcovHC(model, type = "HC0", cluster = ~ cluster)
  
  ape <- avg_slopes(model = model, vcov = vcv_robust) %>% 
    as.data.frame() %>%
    mutate(outcome = str_squish(word(frml,1,sep = "\\~"))) %>%
    filter(grepl("male_std", term) | grepl("delta_w1_w4_GE", term)) 
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
    filter(grepl("male_std", term) | grepl("delta_w1_w4_GE", term)) 
  
  return(ape)
}

library(ggeffects)


model <- glm(dx_htn5 ~ w4.GE_male_std*w5_bp + race5 + sespc_al + nhood1_d + ins5 + edu5,
             data = subset_df, 
             weights = subset_weights, family = "quasibinomial")

vcv_robust <- vcovHC(model, type = "HC1", cluster = ~ cluster)

plot <- ggpredict(model,
            vcov.type = "HC0", 
            vcov.args = list(cluster = subset_df$cluster),
          terms = c("w4.GE_male_std [-2:2, by=0.5]", "w5_bp [0,1]")) %>%
  as.data.frame() %>%
  select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high, group)


ggplot(plot, aes(x=xvals)) +
  geom_ribbon(aes(ymin=lower,ymax=upper), alpha=0.5) +
  geom_line(aes(y=coef)) + 
  #scale_fill_nejm() +
  facet_wrap(~group, ncol = 2) +
  theme_bw() +
  geom_vline(xintercept =0, color='darkred') +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits =c(0, 0.9)) +
  labs(fill = "Wave IV Behavior",
       y = "Predicted Probability\n", 
       x='\nStandardized Gender Expression\n',
       title = 'Marginal effects of male gender expression (MGE) in adolescence and young adulthood on predicted adult substance use\n',
       caption = str_wrap('\n\nIncreases in adolescent  (Wave I) MGE and young adult (Wave IV) MGE 
                          are associated with higher predicted probabilities of young adult 
                          substance use. Note: Adolescence, Participants aged 12-18; Young 
                          adulthood, Participants aged 24-32', 140)) +
  theme(axis.text.x = element_text(color = "black", size = 18),
        axis.text.y = element_text(color = "black", size = 18),
        axis.title.x = element_text(color = 'black',size = 18),
        axis.title.y = element_text(color = 'black', size = 18),
        plot.title = element_text(color = "black", size = 18, hjust = 0),
        plot.caption = element_text(color = "black", size = 15, hjust = 0),
        legend.position = 'right',
        legend.title = element_text(size=18),
        strip.text.x = element_text(color = 'black', size = 18, face = "bold"),
        legend.text = element_text(color = "black", size = 18))


library(ggeffects)


ggpredict(model,vcv_robust, terms = c("w1.GE_male_std [-2:2]", "w5_bp [0,1]"))


# Define a function to get predictions
get_predictions <- function(model, predictor) {
  ggpredict(model, terms = c(paste(predictor, "[-2:2]")), # Define a function to get predictions
get_predictions <- function(model, predictor) {
  ggpredict(model, terms = c(paste(predictor, "[-2:2]")), vcov.type = "HC0", 
            vcov.args = list(cluster = subset(analytical_dataset, 
                                              in_sample == 1)$cluster)) %>%
    as.data.frame() %>%
    dplyr::select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high) %>%
    filter(!is.na(xvals))
}
, 
            vcov.args = list(cluster = subset(analytical_dataset, 
                                              in_sample == 1)$cluster)) %>%
    as.data.frame() %>%
    dplyr::select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high) %>%
    filter(!is.na(xvals))
}


#-------------------------------------------------------------------------------
################################################################################
# Figure 1. Associations of Male Gender Expression (MGE) with Predicted 
# Probability of Wave V CVD-Risk Awareness in Participants with elevated 
# Biomeasures
################################################################################
#-------------------------------------------------------------------------------

# Define the subset and weights
subset_df <- subset(final.df, in_sample.bio == 1)
subset_weights <- subset_df$weights


formulas_htn <- c(
  'dx_htn5 ~ w1.GE_male_std*w5_bp + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'dx_htn5 ~ w4.GE_male_std*w5_bp + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'dx_htn5 ~ delta_w1_w4_GE*w5_bp + race5 + sespc_al + nhood1_d + ins5 + edu5'
)

formulas_dm <- c(
  'dx_dm5 ~ w1.GE_male_std*w5_a1c + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'dx_dm5 ~ w4.GE_male_std*w5_a1c + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'dx_dm5 ~ delta_w1_w4_GE*w5_a1c + race5 + sespc_al + nhood1_d + ins5 + edu5'
)

formulas_hld <- c(
  'dx_hld5 ~ w1.GE_male_std*w5_nhdl + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'dx_hld5 ~ w4.GE_male_std*w5_nhdl + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'dx_hld5 ~ delta_w1_w4_GE*w5_nhdl + race5 + sespc_al + nhood1_d + ins5 + edu5'
)


# Perform analysis 
results_htn <- map_dfr(
  formulas_htn, ~ ape_analysis_by(.x, subset_df, subset_weights, 'w5_bp')
) %>%
  mutate(outcome = 'Hypertension') %>%
  filter(w5_bp == 1)

results_dm <- map_dfr(
  formulas_dm, ~ ape_analysis_by(.x, subset_df, subset_weights, 'w5_a1c')
) %>%
  mutate(outcome = 'Diabetes') %>%
  filter(w5_a1c == 1)

results_hld <- map_dfr(
  formulas_hld, ~ ape_analysis_by(.x, subset_df, subset_weights, 'w5_nhdl')
) %>%
  mutate(outcome = 'Hyperlipidemia') %>%
  filter(w5_nhdl == 1)

# append data frames
results <- bind_rows(results_hld, results_htn, results_dm)

model.2 <- results %>%
  mutate(term = case_when(
    grepl("w1.GE_male_std", term) ~ "Adolescent MGE",
    grepl("w4.GE_male_std", term) ~ "Young Adult MGE",
    TRUE ~ "Change in MGE"
  ),
  outcome = factor(outcome, levels = c("Hypertension", "Diabetes", "Hyperlipidemia")),
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


model.2





