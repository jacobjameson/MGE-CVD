#-------------------------------------------------------------------------
# AUTHOR:             Jacob Jameson
# PURPOSE:            Generate main regression results
#-------------------------------------------------------------------------
#
# load packages ----------------------------------------------------------
rm(list = ls())

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

source('src/prepare data/Construct Analytical Dataset.R')
#----------------------------------------------------------------------------------
##########################################################################
# Wave 5 Dx (Model 1)
#
# Logit models of Dx. All models use cluster robust standard errors
# all models control for race, insurance, education, SES, neighborhood
##########################################################################

# Function to perform regression and coefficient testing
run_analysis <- function(data, predictor) {
  
  # run regressions
  glm_htn <- glm(dx_htn5 ~ get(predictor) + race5 + sespc_al + 
                           nhood1_d + ins5 + edu5, 
                 data = subset(data, in_sample.5 == 1), 
                 weights = weights, family = "quasibinomial")
  
  glm_dm <- glm(dx_dm5 ~ get(predictor) + race5 + sespc_al + 
                         nhood1_d + ins5 + edu5,
                data = subset(data, in_sample.5 == 1), 
                weights = weights, family = "quasibinomial")
  
  glm_hld <- glm(dx_hld5 ~ get(predictor) + race5 + sespc_al +
                           nhood1_d + ins5 + edu5, 
                 data = subset(data, in_sample.5 == 1), 
                 weights = weights, family = "quasibinomial")
  
  # cluster robust standard errors
  vcv_robust_htn <- coeftest(
    x = glm_htn, 
    vcov = vcovCL(glm_htn, type = "HC0", cluster = ~ cluster)
    )
  
  vcv_robust_dm <- coeftest(
    x = glm_dm, 
    vcov = vcovCL(glm_dm, type = "HC0", cluster = ~ cluster)
    )
  
  vcv_robust_hld <- coeftest(
    x = glm_hld, 
    vcov = vcovCL(glm_hld, type = "HC0", cluster = ~ cluster)
    )
  
  # stargazer table of results
  stargazer(list(vcv_robust_htn, vcv_robust_dm, vcv_robust_hld), type = "text", 
            report = "vc*p", header = FALSE, 
            column.labels = c('Diagnosis HTN', 'Diagnosis DM', 'Diagnosis HLD'),
            title = paste('Dx Wave 5 with', predictor))
}

ape_analysis <- function(frml, df, weights){
  
  frml1 <- as.formula(frml)
  m0 <- glm(formula = frml1, data = df, 
            weights = weights, family = "quasibinomial")
  
  vcv_robust <- vcovHC(m0, type = "HC0", cluster = ~ cluster)
  
  ape <- avg_slopes(model = m0, vcov = vcv_robust) %>% 
    as.data.frame() %>%
    mutate(outcome = str_squish(word(frml,1,sep = "\\~"))) %>%
    filter( grepl("male_std", term) | grepl("delta_w1_w4_GE", term)) %>%
    select(outcome, term, estimate, std.error, p.value)
  
  return(ape)
}

# Start the sink to capture output
sink("outputs/tables/Wave 5 Dx Model 1.txt")

predictors <- c("w1.GE_male_std", "w4.GE_male_std", "delta_w1_w4_GE")

for (predictor in predictors) {
  run_analysis(final.df, predictor)
}


ape_analysis(
  'dx_htn5 ~ w1.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5', 
  subset(final.df, in_sample.5 == 1), 
  subset(final.df, in_sample.5 == 1)$weights
  )

ape_analysis(
  'dx_dm5 ~ w1.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5', 
  subset(final.df, in_sample.5 == 1), 
  subset(final.df, in_sample.5 == 1)$weights
  )

ape_analysis(
  'dx_hld5 ~ w1.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5', 
  subset(final.df, in_sample.5 == 1), 
  subset(final.df, in_sample.5 == 1)$weights
  )

ape_analysis(
  'dx_htn5 ~ w4.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5', 
  subset(final.df, in_sample.5 == 1), 
  subset(final.df, in_sample.5 == 1)$weights
)

ape_analysis(
  'dx_dm5 ~ w4.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5', 
  subset(final.df, in_sample.5 == 1), 
  subset(final.df, in_sample.5 == 1)$weights
)

ape_analysis(
  'dx_hld5 ~ w4.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5', 
  subset(final.df, in_sample.5 == 1), 
  subset(final.df, in_sample.5 == 1)$weights
)

ape_analysis(
  'dx_htn5 ~ delta_w1_w4_GE + race5 + sespc_al + nhood1_d + ins5 + edu5', 
  subset(final.df, in_sample.5 == 1), 
  subset(final.df, in_sample.5 == 1)$weights
)

ape_analysis(
  'dx_dm5 ~ delta_w1_w4_GE + race5 + sespc_al + nhood1_d + ins5 + edu5', 
  subset(final.df, in_sample.5 == 1), 
  subset(final.df, in_sample.5 == 1)$weights
)

ape_analysis(
  'dx_hld5 ~ delta_w1_w4_GE + race5 + sespc_al + nhood1_d + ins5 + edu5', 
  subset(final.df, in_sample.5 == 1), 
  subset(final.df, in_sample.5 == 1)$weights
)

sink()

################################################################################
# ------------------------------------------------------------------------------
# Wave 5 Dx (Model 2) stratified by biomarkers
# Logit models of Dx. All models use cluster robust standard errors
# all models control for race, insurance, education, SES, neighborhood
# ------------------------------------------------------------------------------
################################################################################

# Function to perform regression and coefficient testing
run_analysis <- function(data, predictor) {
  
  # run regressions
  glm_htn <- glm(dx_htn5 ~ get(predictor)*w5_bp + race5 + sespc_al + 
                           nhood1_d + ins5 + edu5, 
                 data = subset(data, in_sample.5 == 1), 
                 weights = weights, family = "quasibinomial")
  
  glm_dm <- glm(dx_dm5 ~ get(predictor)*w5_a1c + race5 + sespc_al + 
                         nhood1_d + ins5 + edu5, 
                data = subset(data, in_sample.5 == 1), 
                weights = weights, family = "quasibinomial")
  
  glm_hld <- glm(dx_hld5 ~ get(predictor)*w5_nhdl + race5 + sespc_al + 
                           nhood1_d + ins5 + edu5, 
                 data = subset(data, in_sample.5 == 1), 
                 weights = weights, family = "quasibinomial")
  
  # cluster robust standard errors
  vcv_robust_htn <- coeftest(
    x = glm_htn, 
    vcov = vcovCL(glm_htn, type = "HC0", cluster = ~ cluster)
    )
  
  vcv_robust_dm <- coeftest(
    x = glm_dm, 
    vcov = vcovCL(glm_dm, type = "HC0", cluster = ~ cluster)
    )
  
  vcv_robust_hld <- coeftest(
    x = glm_hld, 
    vcov = vcovCL(glm_hld, type = "HC0", cluster = ~ cluster)
    )
  
  # output results
  stargazer(list(vcv_robust_htn, vcv_robust_dm, vcv_robust_hld), type = "text", 
            report = "vc*p", header = FALSE, 
            column.labels = c('Diagnosis HTN', 'Diagnosis DM', 'Diagnosis HLD'),
            title = paste('Dx Wave 5 with', predictor))
}

ape_analysis <- function(frml, df, weights, var, by){
  
  frml1 <- as.formula(frml)
  m0 <- glm(formula = frml1, data = df, 
            weights = weights, family = "quasibinomial")
  
  vcv_robust <- vcovHC(m0, type = "HC1", cluster = ~ cluster)
  
  ape <- avg_slopes(model = m0, vcov = vcv_robust, 
                    variables = var, by = by) %>%
    as.data.frame() %>%
    mutate(outcome = str_squish(word(frml,1,sep = "\\~"))) %>%
    select(outcome, term, estimate, std.error, p.value, by)
  
  return(ape)
}

# Start the sink to capture output
sink("outputs/tables/Wave 5 Dx Controling for Bio Model 2.txt")

predictors <- c("w1.GE_male_std", "w4.GE_male_std", "delta_w1_w4_GE")
for (predictor in predictors) {
  run_analysis(final.df, predictor)
}


ape_analysis(
  'dx_htn5 ~ w1.GE_male_std*w5_bp + race5 + sespc_al + nhood1_d + ins5 + edu5', 
  subset(final.df, in_sample.5 == 1), 
  subset(final.df, in_sample.5 == 1)$weights,
  'w1.GE_male_std', 'w5_bp')

ape_analysis(
  'dx_dm5 ~ w1.GE_male_std*w5_a1c + race5 + sespc_al + nhood1_d + ins5 + edu5', 
  subset(final.df, in_sample.5 == 1), 
  subset(final.df, in_sample.5 == 1)$weights,
  'w1.GE_male_std', 'w5_a1c')

ape_analysis(
  'dx_hld5 ~ w1.GE_male_std*w5_nhdl + race5 + sespc_al + nhood1_d + ins5 + edu5', 
  subset(final.df, in_sample.5 == 1), 
  subset(final.df, in_sample.5 == 1)$weights,
  'w1.GE_male_std', 'w5_nhdl')

ape_analysis(
  'dx_htn5 ~ w4.GE_male_std*w5_bp + race5 + sespc_al + nhood1_d + ins5 + edu5', 
  subset(final.df, in_sample.5 == 1), 
  subset(final.df, in_sample.5 == 1)$weights,
  'w4.GE_male_std', 'w5_bp')

ape_analysis(
  'dx_dm5 ~ w4.GE_male_std*w5_a1c + race5 + sespc_al + nhood1_d + ins5 + edu5', 
  subset(final.df, in_sample.5 == 1), 
  subset(final.df, in_sample.5 == 1)$weights,
  'w4.GE_male_std', 'w5_a1c')

ape_analysis(
  'dx_hld5 ~ w4.GE_male_std*w5_nhdl + race5 + sespc_al + nhood1_d + ins5 + edu5', 
  subset(final.df, in_sample.5 == 1), 
  subset(final.df, in_sample.5 == 1)$weights,
  'w4.GE_male_std', 'w5_nhdl')

ape_analysis(
  'dx_htn5 ~ delta_w1_w4_GE*w5_bp + race5 + sespc_al + nhood1_d + ins5 + edu5', 
  subset(final.df, in_sample.5 == 1), 
  subset(final.df, in_sample.5 == 1)$weights,
  'delta_w1_w4_GE', 'w5_bp')

ape_analysis(
  'dx_dm5 ~ delta_w1_w4_GE*w5_a1c + race5 + sespc_al + nhood1_d + ins5 + edu5', 
  subset(final.df, in_sample.5 == 1), 
  subset(final.df, in_sample.5 == 1)$weights,
  'delta_w1_w4_GE', 'w5_a1c')

ape_analysis(
  'dx_hld5 ~ delta_w1_w4_GE*w5_nhdl + race5 + sespc_al + nhood1_d + ins5 + edu5', 
  subset(final.df, in_sample.5 == 1), 
  subset(final.df, in_sample.5 == 1)$weights,
  'delta_w1_w4_GE', 'w5_nhdl')


sink()

# Stratified -----------------------------------------------------------------------------------





suppressWarnings({
    for (i in seq_along(outcomes)) {
      for (j in seq_along(predictors)) {
        model_name <- paste0('logit', outcomes[i][[1]][1], j)
        assign(model_name, logit_analysis(
          frml = paste(outcomes[i][[1]][1], "~ race5 + sespc_al + nhood1_d + ins5 + edu5 + ", predictors[j]),
          df = filter(subset(final.df, final.df[[outcomes[[i]][2]]] == 0), 
                      in_sample.bio == 1), weights))
      }}}
)

logit.models.n <- ls()[sapply(ls(), function(x) is.data.frame(get(x))) & grepl("^logit", ls())]
logit.models.n <- do.call(rbind, mget(logit.models.n))

mod2.n <- logit.models.n %>%
  mutate(model_name = 'Model 2',
         biomarker = 'Normal Biomarker')

rm(list=setdiff(ls(), c('final.df', 'logit_analysis', 'mod2.n',
                        'outcomes', 'predictors', 'ape_analysis')))

suppressWarnings({
  for (i in seq_along(outcomes)) {
    for (j in seq_along(predictors)) {
      model_name <- paste0('logit', outcomes[i][[1]][1], j)
      assign(model_name, logit_analysis(
        frml = paste(outcomes[i][[1]][1], "~ race5 + sespc_al + nhood1_d + ins5 + edu5 + ", predictors[j]),
        df = filter(subset(final.df, final.df[[outcomes[[i]][2]]] == 1), 
                    in_sample.bio == 1), weights))
    }}}
)

logit.models.a <- ls()[sapply(ls(), function(x) is.data.frame(get(x))) & grepl("^logit", ls())]
logit.models.a <- do.call(rbind, mget(logit.models.a))

mod2.a <- logit.models.a %>%
  mutate(model_name = 'Model 2',
         biomarker = 'Abormal Biomarker')

logit.models <- rbind(mod2.n, mod2.a)

write_csv(logit.models, 'tables:figures/Logit Models 2 Stratified.csv')

# Clear environment
rm(list=setdiff(ls(), c('final.df', 'logit_analysis', 'ape_analysis')))


# Interaction Term Regression w/ Average Marginal Effects ---------------------------------------------
#------------------------------------------------------------------------------------------------------
m1 <- glm(formula = dx_htn5 ~ race5 + sespc_al + nhood1_d + ins5 + edu5 + w1.GE_male_std*w5_bp,
          data = filter(final.df, in_sample.bio == 1), 
          weights = weights.bio, family = "quasibinomial")

vcv_robust <- vcovHC(m1, type = "HC0", cluster = ~ cluster)

ape.1 <- slopes(model = m1, vcov = vcv_robust,  newdata = datagrid(w5_bp = c(0, 1))) %>% 
  as.data.frame() %>%
  mutate(outcome = 'dx_htn5') %>%
  filter( grepl("w4.GE_male_std", term) | grepl("w1.GE_male_std", term)) %>%
  dplyr::select(c("outcome", "term", biomarker_abnormal = w5_bp, "contrast", 
                  "estimate", "std.error", "statistic",
                  "p.value", "conf.low", "conf.high", "predicted", 
                  "predicted_hi", "predicted_lo"))


m1 <- glm(formula = dx_htn5 ~ race5 + sespc_al + nhood1_d + ins5 + edu5 + w4.GE_male_std*w5_bp,
          data = filter(final.df, in_sample.bio == 1), 
          weights = weights.bio, family = "quasibinomial")

vcv_robust <- vcovHC(m1, type = "HC0", cluster = ~ cluster)

ape.1b <- slopes(model = m1, vcov = vcv_robust,  newdata = datagrid(w5_bp = c(0, 1))) %>% 
  as.data.frame() %>%
  mutate(outcome = 'dx_htn5') %>%
  filter( grepl("w4.GE_male_std", term) | grepl("w1.GE_male_std", term)) %>%
  dplyr::select(c("outcome", "term", biomarker_abnormal = w5_bp, "contrast", 
                  "estimate", "std.error", "statistic",
                  "p.value", "conf.low", "conf.high", "predicted", 
                  "predicted_hi", "predicted_lo"))

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
m2 <- glm(formula = dx_dm5 ~ race5 + sespc_al + nhood1_d + ins5 + edu5 + w1.GE_male_std*w5_a1c,
          data = filter(final.df, in_sample.bio == 1), 
          weights = weights.bio, family = "quasibinomial")

vcv_robust <- vcovHC(m2, type = "HC0", cluster = ~ cluster)

ape.2 <- slopes(model = m2, vcov = vcv_robust,  newdata = datagrid(w5_a1c = c(0, 1))) %>% 
  as.data.frame() %>%
  mutate(outcome = 'dx_dm5') %>%
  filter( grepl("w4.GE_male_std", term) | grepl("w1.GE_male_std", term)) %>%
  dplyr::select(c("outcome", "term", biomarker_abnormal = w5_a1c, "contrast", 
                  "estimate", "std.error", "statistic",
                  "p.value", "conf.low", "conf.high", "predicted", 
                  "predicted_hi", "predicted_lo"))


m2 <- glm(formula = dx_dm5 ~ race5 + sespc_al + nhood1_d + ins5 + edu5 + w4.GE_male_std*w5_a1c,
          data = filter(final.df, in_sample.bio == 1), 
          weights = weights.bio, family = "quasibinomial")

vcv_robust <- vcovHC(m2, type = "HC0", cluster = ~ cluster)

ape.2b <- slopes(model = m2, vcov = vcv_robust,  newdata = datagrid(w5_a1c = c(0, 1))) %>% 
  as.data.frame() %>%
  mutate(outcome = 'dx_dm5') %>%
  filter( grepl("w4.GE_male_std", term) | grepl("w1.GE_male_std", term)) %>%
  dplyr::select(c("outcome", "term", biomarker_abnormal = w5_a1c, "contrast", 
                  "estimate", "std.error", "statistic",
                  "p.value", "conf.low", "conf.high", "predicted", 
                  "predicted_hi", "predicted_lo"))

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
m3 <- glm(formula = dx_hld5 ~ race5 + sespc_al + nhood1_d + ins5 + edu5 + w1.GE_male_std*w5_nhdl,
          data = filter(final.df, in_sample.bio == 1), 
          weights = weights.bio, family = "quasibinomial")

vcv_robust <- vcovHC(m3, type = "HC0", cluster = ~ cluster)

ape.3 <- slopes(model = m3, vcov = vcv_robust,  newdata = datagrid(w5_nhdl = c(0, 1))) %>% 
  as.data.frame() %>%
  mutate(outcome = 'dx_hld5') %>%
  filter( grepl("w4.GE_male_std", term) | grepl("w1.GE_male_std", term)) %>%
  dplyr::select(c("outcome", "term", biomarker_abnormal = w5_nhdl, "contrast", 
                  "estimate", "std.error", "statistic",
                  "p.value", "conf.low", "conf.high", "predicted", 
                  "predicted_hi", "predicted_lo"))


m3 <- glm(formula = dx_hld5 ~ race5 + sespc_al + nhood1_d + ins5 + edu5 + w4.GE_male_std*w5_nhdl,
          data = filter(final.df, in_sample.bio == 1), 
          weights = weights.bio, family = "quasibinomial")

vcv_robust <- vcovHC(m3, type = "HC0", cluster = ~ cluster)

ape.3b <- slopes(model = m3, vcov = vcv_robust,  newdata = datagrid(w5_nhdl = c(0, 1))) %>% 
  data.frame() %>%
  mutate(outcome = 'dx_hld5') %>%
  filter( grepl("w4.GE_male_std", term) | grepl("w1.GE_male_std", term)) %>%
  dplyr::select(c("outcome", "term", biomarker_abnormal = w5_nhdl, "contrast", 
                  "estimate", "std.error", "statistic",
                  "p.value", "conf.low", "conf.high", "predicted", 
                  "predicted_hi", "predicted_lo"))


ame.models <- ls()[sapply(ls(), function(x) is.data.frame(get(x))) & grepl("^ape", ls())]
ame.models <- do.call(rbind, mget(ame.models))

write_csv(ame.models, 'tables:figures/Logit Models 2 Modified AME.csv')

# Clear environment
rm(list=setdiff(ls(), c('final.df', 'logit_analysis', 'ape_analysis')))

################################################################################
# ------------------------------------------------------------------------------
# Wave 5 Biomarker (Model 3) 
# Logit models of Dx. All models use cluster robust standard errors
# all models control for race, insurance, education, SES, neighborhood
# ------------------------------------------------------------------------------
################################################################################

# Create a vector of outcome variables
outcomes <- c("w5_bp", "w5_a1c", "w5_nhdl")
meds <- c('w5_anti_htn', 'w5_anti_dm_med_use', 'w5_anti_hld_med_use')

# Create a vector of predictor variables
predictors <- c("w1.GE_male_std", "w4.GE_male_std")

suppressWarnings({
  for (i in seq_along(outcomes)) {
    for (j in seq_along(predictors)) {
      model_name <- paste0('logit', outcomes[i], j)
      assign(model_name, logit_analysis(
        frml = paste(outcomes[i], "~ race5 + sespc_al + nhood1_d + ins5 + edu5 + ", predictors[j], '+ ', meds[i]),
        df = subset(final.df, in_sample.bio == 1), weights.bio))
    }}}
)


logit.models <- ls()[sapply(ls(), function(x) is.data.frame(get(x))) & grepl("^logit", ls())]
logit.models <- do.call(rbind, mget(logit.models))

logit.models <- logit.models %>%
  mutate(model_name = 'Model 3')

write_csv(logit.models, 'tables:figures/Logit Models 3.csv')


suppressWarnings({
  for (i in seq_along(outcomes)) {
    for (j in seq_along(predictors)) {
      model_name <- paste0('AME', outcomes[i], j)
      assign(model_name, ape_analysis(
        frml = paste(outcomes[i], "~ race5 + sespc_al + nhood1_d + ins5 + edu5 + ", predictors[j], '+', meds[i]),
        df = subset(final.df, in_sample.bio == 1), weights.bio))
    }}}
)

ape.models <- ls()[sapply(ls(), function(x) is.data.frame(get(x))) & grepl("^AME", ls())]
ape.models <- do.call(rbind, mget(ape.models))

ape.models <- ape.models %>%
  mutate(model_name = 'Model 3')

write_csv(ape.models, 'tables:figures/AME Models 3.csv')

# Clear environment
rm(list=setdiff(ls(), c('final.df', 'logit_analysis', 'ape_analysis')))

################################################################################
# ------------------------------------------------------------------------------
# Wave 5 Medication Use in Men who Report (+) Diagnosis (Model 4) 
# Logit models of Dx. All models use cluster robust standard errors
# all models control for race, insurance, education, SES, neighborhood
# ------------------------------------------------------------------------------
################################################################################

# Create a vector of outcome variables
outcomes <- c('w5_anti_htn', 'w5_anti_dm_med_use', 'w5_anti_hld_med_use')

diagnoses <- c("dx_htn5", "dx_dm5", "dx_hld5")

# Create a vector of predictor variables
predictors <- c("w1.GE_male_std", "w4.GE_male_std")

suppressWarnings({
  for (i in seq_along(outcomes)) {
    for (j in seq_along(predictors)) {
      model_name <- paste0('logit', outcomes[i], j)
      assign(model_name, logit_analysis(
        frml = paste(outcomes[i], "~ race5 + sespc_al + nhood1_d + ins5 + edu5 + ", predictors[j]),
        df = filter(subset(final.df, final.df[diagnoses[i]] == 1), in_sample.bio == 1), weights.bio))
    }}}
)

logit.models <- ls()[sapply(ls(), function(x) is.data.frame(get(x))) & grepl("^logit", ls())]
logit.models <- do.call(rbind, mget(logit.models))

logit.models <- logit.models %>%
  mutate(model_name = 'Model 4')

write_csv(logit.models, 'tables:figures/Logit Models 4 Stratified dx=1.csv')

# Clear environment
rm(list=setdiff(ls(), c('final.df', 'logit_analysis', 'ape_analysis')))

# Interaction Term Regression w/ Average Marginal Effects ---------------------------------------------
#------------------------------------------------------------------------------------------------------

m1 <- glm(formula = w5_anti_htn ~ race5 + sespc_al + nhood1_d + 
            ins5 + edu5 + w1.GE_male_std*dx_htn5,
          data = filter(final.df, in_sample.bio == 1), 
          weights = weights.bio, family = "quasibinomial")

table(filter(final.df, in_sample.bio == 1)$w5_anti_htn)

vcv_robust <- vcovHC(m1, type = "HC0", cluster = ~ cluster)

summary(m1)
ape.1 <- slopes(model = m1, vcov = vcv_robust,  newdata = datagrid(dx_htn5 = c(0, 1))) %>% 
  as.data.frame() %>%
  mutate(outcome = 'w5_anti_htn') %>%
  filter( grepl("w4.GE_male_std", term) | grepl("w1.GE_male_std", term)) %>%
  dplyr::select(c("outcome", "term", diagnosis = dx_htn5, "contrast", 
                  "estimate", "std.error", "statistic",
                  "p.value", "conf.low", "conf.high", "predicted", 
                  "predicted_hi", "predicted_lo"))


m1 <- glm(formula = w5_anti_htn ~ race5 + sespc_al + nhood1_d + ins5 + edu5 + w4.GE_male_std*dx_htn5,
          data = filter(final.df, in_sample.bio == 1), 
          weights = weights.bio, family = "quasibinomial")

vcv_robust <- vcovHC(m1, type = "HC0", cluster = ~ cluster)

ape.1b <- slopes(model = m1, vcov = vcv_robust,  newdata = datagrid(dx_htn5 = c(0, 1))) %>% 
  as.data.frame() %>%
  mutate(outcome = 'w5_anti_htn') %>%
  filter( grepl("w4.GE_male_std", term) | grepl("w1.GE_male_std", term)) %>%
  dplyr::select(c("outcome", "term", diagnosis = dx_htn5, "contrast", 
                  "estimate", "std.error", "statistic",
                  "p.value", "conf.low", "conf.high", "predicted", 
                  "predicted_hi", "predicted_lo"))

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
m2 <- glm(formula = w5_anti_dm_med_use ~ race5 + sespc_al + nhood1_d + ins5 + edu5 + w1.GE_male_std*dx_dm5,
          data = filter(final.df, in_sample.bio == 1), 
          weights = weights.bio, family = "quasibinomial")

vcv_robust <- vcovHC(m2, type = "HC0", cluster = ~ cluster)

ape.2 <- slopes(model = m2, vcov = vcv_robust,  newdata = datagrid(dx_dm5 = c(0, 1))) %>% 
  as.data.frame() %>%
  mutate(outcome = 'w5_anti_dm_med_use') %>%
  filter( grepl("w4.GE_male_std", term) | grepl("w1.GE_male_std", term)) %>%
  dplyr::select(c("outcome", "term", diagnosis = dx_dm5, "contrast", 
                  "estimate", "std.error", "statistic",
                  "p.value", "conf.low", "conf.high", "predicted", 
                  "predicted_hi", "predicted_lo"))


m2 <- glm(formula = w5_anti_dm_med_use ~ race5 + sespc_al + nhood1_d + ins5 + edu5 + w4.GE_male_std*dx_dm5,
          data = filter(final.df, in_sample.bio == 1), 
          weights = weights.bio, family = "quasibinomial")

vcv_robust <- vcovHC(m2, type = "HC0", cluster = ~ cluster)

ape.2b <- slopes(model = m2, vcov = vcv_robust,  newdata = datagrid(dx_dm5 = c(0, 1))) %>% 
  as.data.frame() %>%
  mutate(outcome = 'w5_anti_dm_med_use') %>%
  filter( grepl("w4.GE_male_std", term) | grepl("w1.GE_male_std", term)) %>%
  dplyr::select(c("outcome", "term", diagnosis = dx_dm5, "contrast", 
                  "estimate", "std.error", "statistic",
                  "p.value", "conf.low", "conf.high", "predicted", 
                  "predicted_hi", "predicted_lo"))

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

m3 <- glm(formula = w5_anti_hld_med_use ~ race5 + sespc_al + nhood1_d + ins5 + edu5 + w1.GE_male_std*dx_hld5,
          data = filter(final.df, in_sample.bio == 1), 
          weights = weights.bio, family = "quasibinomial")

vcv_robust <- vcovHC(m3, type = "HC0", cluster = ~ cluster)

ape.3 <- slopes(model = m3, vcov = vcv_robust,  newdata = datagrid(dx_hld5 = c(0, 1))) %>% 
  as.data.frame() %>%
  mutate(outcome = 'w5_anti_hld_med_use') %>%
  filter( grepl("w4.GE_male_std", term) | grepl("w1.GE_male_std", term)) %>%
  dplyr::select(c("outcome", "term", diagnosis = dx_hld5, "contrast", 
                  "estimate", "std.error", "statistic",
                  "p.value", "conf.low", "conf.high", "predicted", 
                  "predicted_hi", "predicted_lo"))


m3 <- glm(formula = w5_anti_hld_med_use ~ race5 + sespc_al + nhood1_d + ins5 + edu5 + w4.GE_male_std*dx_hld5,
          data = filter(final.df, in_sample.bio == 1), 
          weights = weights.bio, family = "quasibinomial")

vcv_robust <- vcovHC(m3, type = "HC0", cluster = ~ cluster)

ape.3b <- slopes(model = m3, vcov = vcv_robust,  newdata = datagrid(dx_hld5 = c(0, 1))) %>% 
  data.frame() %>%
  mutate(outcome = 'w5_anti_hld_med_use') %>%
  filter( grepl("w4.GE_male_std", term) | grepl("w1.GE_male_std", term)) %>%
  dplyr::select(c("outcome", "term", diagnosis = dx_hld5, "contrast", 
                  "estimate", "std.error", "statistic",
                  "p.value", "conf.low", "conf.high", "predicted", 
                  "predicted_hi", "predicted_lo"))


ame.models <- ls()[sapply(ls(), function(x) is.data.frame(get(x))) & grepl("^ape", ls())]
ame.models <- do.call(rbind, mget(ame.models))

write_csv(ame.models, 'tables:figures/AME Models 4.csv')



