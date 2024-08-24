################################################################################
# AUTHOR:             Jacob Jameson
#
# DESCRIPTION:        This script generates the tables for the paper.
# 
# DEPENDENCIES:       The script requires the following dataset:
#                     -final.df, which is the final dataset constructed by 
#                      executing scr/prepare data/Construct Analytical dataset.R
################################################################################


################################################################################
#-------------------------------------------------------------------------------
# Calculate confidence intervals for the prevalence of elevated biomeasure
#-------------------------------------------------------------------------------
################################################################################

ahdsgn <- svydesign(
  id=~psuscid,
  strata=~region,
  weights=~w5biowgt,
  data=subset(
    final.df, in_sample.bio == 1 & is.na(psuscid) == F & is.na(region) == F
  ),
  nest=TRUE)

bp <- data.frame(svymean(~w5_bp, ahdsgn, na.rm=TRUE))
bp <- rename(bp, se = w5_bp)
bp$ci <- paste('[',round(bp$mean - 1.96 * bp$se, 3), ',',
               round(bp$mean + 1.96 * bp$se, 3),']')


db <- data.frame(svymean(~w5_a1c, ahdsgn, na.rm=TRUE))
db <- rename(db, se = w5_a1c)
db$ci <- paste('[',round(db$mean - 1.96 * db$se, 3), ',',
               round(db$mean + 1.96 * db$se, 3),']')


nhdl <- data.frame(svymean(~w5_nhdl, ahdsgn, na.rm=TRUE))
nhdl <- rename(nhdl, se = w5_nhdl)
nhdl$ci <- paste('[',round(nhdl$mean - 1.96 * nhdl$se, 3), ',',
                 round(nhdl$mean + 1.96 * nhdl$se, 3),']')


# put it all into a table
biomeasures <- rbind(bp, db, nhdl)
biomeasures$biomeasure <- c('Blood Pressure (weighted)', 'A1c (weighted)', 'Non-HDL Cholesterol (weighted')
biomeasures <- select(biomeasures, biomeasure, mean, se, ci)
biomeasures <- rename(biomeasures, prevalence = mean)
rownames(biomeasures) <- NULL
biomeasures

################################################################################
#-------------------------------------------------------------------------------
# Check missingness in the dataset
#-------------------------------------------------------------------------------
################################################################################

vars <- c('race5', 'sespc_al', 'nhood1_d', 'ins5', 'edu5',
          'dx_htn5', 'dx_dm5', 'dx_hld5')

vars.bio <- c('w5_anti_htn', 'w5_anti_dm_med_use', 'w5_anti_hld_med_use',
              'w5_bp', 'w5_a1c', 'w5_nhdl')

# Creating a detailed summary table
summary_table <- final.df %>%
  filter(in_sample.5 == 1) %>%
  summarise(
    Total_N = n(),  # Total number of observations in the sample
    across(all_of(vars), ~sum(is.na(.)) / Total_N * 100, .names = "perc_missing_{.col}"),  # Percentage missing for each variable
    Perc_Missing_Any = round(sum(rowSums(is.na(select(., all_of(vars)))) > 0) / Total_N * 100, 2)  # Percentage of rows missing at least one variable
  ) %>%
  pivot_longer(cols = -Total_N, 
               names_to = "Variable", 
               values_to = "Percentage_Missing") 

summary_table$Percentage_Missing <- paste(round(summary_table$Percentage_Missing, 2), '%')

summary_table.bio <- final.df %>%
  filter(in_sample.bio == 1) %>%
  summarise(
    Total_N = n(),  # Total number of observations in the sample
    across(all_of(vars.bio), ~sum(is.na(.)) / Total_N * 100, .names = "perc_missing_{.col}"),  # Percentage missing for each variable
    Perc_Missing_Any = sum(rowSums(is.na(select(., all_of(vars.bio)))) > 0) / Total_N * 100  # Percentage of rows missing at least one variable
  ) %>%
  pivot_longer(cols = -Total_N, names_to = "Variable", values_to = "Percentage_Missing")  # Reshape to long format




# Formatting and saving the table
table_output <- kable(summary_table, 
                      caption = "Summary of Missing Data in Sample") %>%
  kable_styling(font_size = 12, full_width = F) %>%
  add_header_above(c(" " = 1, "Missing Data Summary" = 2))  

# Save the table as a PDF file
save_kable(table_output, file = "outputs/tables/Missing_Data_Summary.pdf")


#-------------------------------------------------------------------------------
################################################################################
# Table 1. Demographic Characteristics of Males Included in 
# Analytic Sample (N = 4,230)
################################################################################
#-------------------------------------------------------------------------------
library(gtsummary)
theme_gtsummary_journal(journal = "jama")


final.df$group <- ifelse(final.df$w4.GE_male_std >=0 & final.df$in_sample.5 == 1,
                   'Below Average MGE in YA',
                   'Above Average MGE in YA')

vars <- c('race', 'edu5', 'ins5', 'sespc_al', 'nhood1_d', 
          'dx_htn5', 'w5_bp', 'w5_anti_htn',
          'dx_dm5', 'w5_a1c', 'w5_anti_dm_med_use',
          'dx_hld5', 'w5_nhdl', 'w5_anti_hld_med_use')

factorize <- c('dx_htn5', 'w5_bp', 'w5_anti_htn',
               'dx_dm5', 'w5_a1c', 'w5_anti_dm_med_use',
               'dx_hld5', 'w5_nhdl', 'w5_anti_hld_med_use')

final.df[factorize] <- lapply(final.df[factorize], factor)


final.df %>%
  filter(in_sample == 1) %>%
  tbl_summary(
    by = group, 
    type = all_continuous() ~ "continuous2",
    statistic = 
      all_continuous() ~ c("{median} ({min}, {max})"),
    missing = "no",
    include = c(vars, group)) %>%
  add_n() %>%
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 2),
        test = list(all_continuous() ~ "t.test",
                    all_categorical() ~ "chisq.test")) %>%
  add_overall() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_caption("**Table 1. Patient Characteristics (unweighted)**") %>%
  bold_labels()

# save table
final.df %>%
  filter(in_sample == 1) %>%
  tbl_summary(
    by = group, 
    type = all_continuous() ~ "continuous2",
    statistic = 
      all_continuous() ~ c("{median} ({min}, {max})"),
    missing = "no",
    include = c(vars, group)) %>%
  add_n() %>%
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 2),
        test = list(all_continuous() ~ "t.test",
                    all_categorical() ~ "chisq.test")) %>%
  add_overall() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_caption("**Table 1. Patient Characteristics (unweighted)**") %>%
  bold_labels() %>%
  as_kable() %>%
  kable_styling(full_width = F) %>%
  save_kable(file = "outputs/tables/Table1.pdf")

ahdsgn <- svydesign(
  id=~psuscid,
  strata=~region,
  weights=~gsw5,
  data=subset(final.df,  in_sample.5 == 1 & is.na(psuscid) == F & is.na(region) == F),
  nest=TRUE)


ahdsgn %>%
  tbl_svysummary(
    by = group, 
    type = all_continuous() ~ "continuous2",
    statistic = 
      all_continuous() ~ c("{median} ({min}, {max})"),
    missing = "always",
    include = c(ins5, group)) %>%
  # show percentage of missing data
  add_n() %>%
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 2),
        test = list(all_continuous() ~ "svy.t.test",
                    all_categorical() ~ "svy.chisq.test")) %>%
  add_overall() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_caption("**Table 1. Patient Characteristics (weighted)**") %>%
  bold_labels()




controls <- 'race5 + sespc_al + nhood1_d + ins5 + edu5'

final.df <- final.df %>%
  group_by(in_sample) %>%
  mutate_at(c('race5', 'sespc_al', 'nhood1_d', 'ins5', 'edu5'), funs(replace(., is.na(.), median(., na.rm = TRUE))) ) %>%
  ungroup()

final.df$complete_case <- ifelse(
  rowSums(
    is.na(
      final.df[, c('race5', 'sespc_al', 'nhood1_d', 'ins5', 'edu5')]
    )
  ) == 0, 1, 0)



#-------------------------------------------------------------------------------
################################################################################
# Table 2. Marginal effects coefficients (dy/dx) estimating associations 
# between adolescent, young adult, and adolescent-to-young-adult change 
# in male gender expression (MGE) and adult awareness of hypertension, 
# diabetes, and hyperlipidemia 
################################################################################
#-------------------------------------------------------------------------------

library(margins)

ame_analysis <- function(frml, design, var, outcome){
  
  frml1 <- as.formula(frml)
  model <- svyglm(formula = frml1, design = design, 
                  family = "quasibinomial")
  
  ame <-
    summary(
      margins(model, method = "dydx", variables = c(var), design = design)
    ) %>%
    as.data.frame() %>% 
    mutate(
      outcome = outcome,
      AME = sprintf("%.3f", AME),
      SE = sprintf("%.3f", SE),
      CI = sprintf("(%.3f, %.3f)", lower, upper),
      P_value = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
    ) %>%
    select(Outcome = outcome, Term = factor, AME, `AME SE` = SE, 
           `AME P_value`= P_value, `AME 95% CI` = CI)
  
  return(ame)
}

################################################################################
## Model 1 --------------------------------------------------------------------
### (DX ~ MGE  + controls)
################################################################################

ahdsgn <- svydesign(
  id=~psuscid,
  strata=~region,
  weights=~gsw5,
  data=subset(final.df, is.na(gsw5) == F & is.na(region) == F),
  nest=TRUE)

# hypertension -----------------------------------------------------------------

htn.1a <- ame_analysis(paste('dx_htn5 ~ w1.GE_male_std +', controls),
                       subset(ahdsgn, in_sample.5 == 1 & complete_case == 1), 
                       'w1.GE_male_std', 'dx_htn5')

htn.1b <- ame_analysis(paste('dx_htn5 ~ w4.GE_male_std +', controls),
                       subset(ahdsgn, in_sample.5 == 1 & complete_case == 1),
                       'w4.GE_male_std', 'dx_htn5')

htn.1c <- ame_analysis(paste('dx_htn5 ~ delta_w1_w4_GE + w1.GE_male_std +', controls),
                       subset(ahdsgn, in_sample.5 == 1 & complete_case == 1),
                       'delta_w1_w4_GE', 'dx_htn5')

htn.1 <- rbind(htn.1a, htn.1b, htn.1c)

# diabetes ---------------------------------------------------------------------
dm.1a <- ame_analysis(paste('dx_dm5 ~ w1.GE_male_std + ', controls),
                      subset(ahdsgn, in_sample.5 == 1 & complete_case == 1), 
                      'w1.GE_male_std', 'dx_dm5')

dm.1b <- ame_analysis(paste('dx_dm5 ~ w4.GE_male_std + ', controls),
                      subset(ahdsgn, in_sample.5 == 1 & complete_case == 1), 
                      'w4.GE_male_std', 'dx_dm5')

dm.1c <- ame_analysis(paste('dx_dm5 ~ delta_w1_w4_GE + w1.GE_male_std + ', controls),
                      subset(ahdsgn, in_sample.5 == 1 & complete_case == 1), 
                      'delta_w1_w4_GE', 'dx_dm5')

dm.1 <- rbind(dm.1a, dm.1b, dm.1c)

# hyperlipidemia ---------------------------------------------------------------
hld.1a <- ame_analysis(paste('dx_hld5 ~ w1.GE_male_std + ', controls),
                      subset(ahdsgn, in_sample.5 == 1 & complete_case == 1), 
                      'w1.GE_male_std', 'dx_hld5')

hld.1b <- ame_analysis(paste('dx_hld5 ~ w4.GE_male_std + ', controls),
                       subset(ahdsgn, in_sample.5 == 1 & complete_case == 1),
                      'w4.GE_male_std', 'dx_hld5')

hld.1c <- ame_analysis(paste('dx_hld5 ~ delta_w1_w4_GE + w1.GE_male_std + ', controls),
                       subset(ahdsgn, in_sample.5 == 1 & complete_case == 1), 
                       'delta_w1_w4_GE', 'dx_hld5')

hld.1 <- rbind(hld.1a, hld.1b, hld.1c)

model.1 <- rbind(htn.1, dm.1, hld.1) %>%
  mutate(stars = case_when(
    `AME P_value` < 0.01 ~ "***",
    `AME P_value` < 0.05 ~ "**",
    `AME P_value` < 0.1 ~ "*",
    TRUE ~ ""),
    values = paste(paste0(AME,stars), `AME 95% CI`, 'P Value =', `AME P_value`)) %>%
  select(Outcome, Term, values) %>%
  pivot_wider(names_from = Term, values_from = values) %>%
  mutate(Outcome = case_when(
    Outcome == "dx_htn5" ~ "Dx of hypertension",
    Outcome == "dx_dm5"  ~ "Dx of diabetes",
    Outcome == "dx_hld5" ~ "Dx of hyperlipidemia"
  ),
  Model = 'Model 1') %>%
  rename(`Adolescent MGE` = w1.GE_male_std, 
         `Young adult MGE` = w4.GE_male_std,
         `Change in MGE` = delta_w1_w4_GE) 

print(model.1)

################################################################################
## Model 2 --------------------------------------------------------------------
### (DX ~ MGE##Bio-M + controls
################################################################################

ahdsgnbio <- svydesign(
  id=~psuscid,
  strata=~region,
  weights=~weights.bio,
  data=subset(final.df, in_sample.bio == 1 & is.na(psuscid) == F & is.na(region) == F),
  nest=TRUE)

ame_by_bio_analysis <- function(frml, design, var, outcome, bio_var){
  
  frml1 <- as.formula(frml)
  model <- svyglm(formula = frml1, design = design, 
                  family = "quasibinomial")
  
  if (bio_var == 'w5_bp'){
    ame <-
      summary(
        margins(model, method = "dydx", design = design,
                variables = c(var), 
                at = list(w5_bp = c(1)))
      )
  }
  
  if (bio_var == 'w5_a1c'){
    ame <-
      summary(
        margins(model, method = "dydx", design = design,
                variables = c(var), 
                at = list(w5_a1c = c(1)))
      )
  }
  
  if (bio_var == 'w5_nhdl'){
    ame <-
      summary(
        margins(model, method = "dydx", design = design,
                variables = c(var), 
                at = list(w5_nhdl = c(1)))
      )
  }

  ame <- ame %>%
    as.data.frame() %>% 
    mutate(
      outcome = outcome,
      AME = sprintf("%.3f", AME),
      SE = sprintf("%.3f", SE),
      CI = sprintf("(%.3f, %.3f)", lower, upper),
      P_value = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
    ) %>%
    select(Outcome = outcome, Term = factor, AME, `AME SE` = SE, 
           `AME P_value`= P_value, `AME 95% CI` = CI)
  
  return(ame)
}

# hypertension -----------------------------------------------------------------
htn.2a <- ame_by_bio_analysis(paste('dx_htn5 ~ w1.GE_male_std*w5_bp +', controls),
                              subset(ahdsgnbio, complete_case == 1), 
                              'w1.GE_male_std', 'dx_htn5', 'w5_bp')

htn.2b <- ame_by_bio_analysis(paste('dx_htn5 ~ w4.GE_male_std*w5_bp +', controls),
                              subset(ahdsgnbio, complete_case == 1), 
                              'w4.GE_male_std', 'dx_htn5', 'w5_bp')

htn.2c <- ame_by_bio_analysis(paste('dx_htn5 ~ delta_w1_w4_GE*w5_bp + w1.GE_male_std +', controls),
                              subset(ahdsgnbio, complete_case == 1), 
                              'delta_w1_w4_GE', 'dx_htn5', 'w5_bp')

htn.2 <- rbind(htn.2a, htn.2b, htn.2c)

# diabetes ---------------------------------------------------------------------
dm.2a <- ame_by_bio_analysis(paste('dx_dm5 ~ w1.GE_male_std*w5_a1c +', controls),
                          subset(ahdsgnbio, complete_case == 1), 
                          'w1.GE_male_std', 'dx_dm5', 'w5_a1c')

dm.2b <- ame_by_bio_analysis(paste('dx_dm5 ~ w4.GE_male_std*w5_a1c +', controls),
                          subset(ahdsgnbio, complete_case == 1), 
                          'w4.GE_male_std', 'dx_dm5', 'w5_a1c')

dm.2c <- ame_by_bio_analysis(paste('dx_dm5 ~ delta_w1_w4_GE*w5_a1c + w1.GE_male_std +', controls),
                             subset(ahdsgnbio, complete_case == 1), 
                             'delta_w1_w4_GE', 'dx_dm5', 'w5_a1c')

dm.2 <- rbind(dm.2a, dm.2b, dm.2c)

# hyperlipidemia --------------------------------------------------------------
hld.2a <- ame_by_bio_analysis(paste('dx_hld5 ~ w1.GE_male_std*w5_nhdl +', controls),
                              subset(ahdsgnbio, complete_case == 1),
                              'w1.GE_male_std', 'dx_hld5', 'w5_nhdl')

hld.2b <- ame_by_bio_analysis(paste('dx_hld5 ~ w4.GE_male_std*w5_nhdl +', controls),
                              subset(ahdsgnbio, complete_case == 1),
                              'w4.GE_male_std', 'dx_hld5', 'w5_nhdl')

hld.2c <- ame_by_bio_analysis(paste('dx_hld5 ~ delta_w1_w4_GE*w5_nhdl + w1.GE_male_std +', controls),
                              subset(ahdsgnbio, complete_case == 1),
                              'delta_w1_w4_GE', 'dx_hld5', 'w5_nhdl')

hld.2 <- rbind(hld.2a, hld.2b, hld.2c)

model.2 <- rbind(htn.2, dm.2, hld.2) %>%
  mutate(stars = case_when(
    `AME P_value` < 0.01 ~ "***",
    `AME P_value` < 0.05 ~ "**",
    `AME P_value` < 0.1 ~ "*",
    TRUE ~ ""),
    values = paste(paste0(AME,stars), `AME 95% CI`, 'P Value =', `AME P_value`)) %>%
  select(Outcome, Term, values) %>%
  pivot_wider(names_from = Term, values_from = values) %>%
  mutate(Outcome = case_when(
    Outcome == "dx_htn5" ~ "Dx of hypertension",
    Outcome == "dx_dm5"  ~ "Dx of diabetes",
    Outcome == "dx_hld5" ~ "Dx of hyperlipidemia"
  ),
  Model = 'Model 2') %>%
  rename(`Adolescent MGE` = w1.GE_male_std, 
         `Young adult MGE` = w4.GE_male_std,
         `Change in MGE` = delta_w1_w4_GE) 

print(model.2)

################################################################################
# Model 3 ----------------------------------------------------------------------
### (TX ~ MGE##DX + controls)  
################################################################################

ame_by_dx_analysis <- function(frml, design, var, outcome, dx){
  
  frml1 <- as.formula(frml)
  model <- svyglm(formula = frml1, design = design, 
                  family = "quasibinomial")
  
  if (dx == 'dx_htn5'){
    ame <-
      summary(
        margins(model, method = "dydx", design = design,
                variables = c(var), 
                at = list(dx_htn5 = c(1)))
      )
  }
  
  if (dx == 'dx_dm5'){
    ame <-
      summary(
        margins(model, method = "dydx", design = design,
                variables = c(var), 
                at = list(dx_dm5 = c(1)))
      )
  }
  
  if (dx == 'dx_hld5'){
    ame <-
      summary(
        margins(model, method = "dydx", design = design,
                variables = c(var), 
                at = list(dx_hld5 = c(1)))
      )
  }
  
  ame <- ame %>%
    as.data.frame() %>% 
    mutate(
      outcome = outcome,
      AME = sprintf("%.3f", AME),
      SE = sprintf("%.3f", SE),
      CI = sprintf("(%.3f, %.3f)", lower, upper),
      P_value = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
    ) %>%
    select(Outcome = outcome, Term = factor, AME, `AME SE` = SE, 
           `AME P_value`= P_value, `AME 95% CI` = CI)
  
  return(ame)
}


# hypertension -----------------------------------------------------------------

htn.3a <- ame_by_dx_analysis(paste('w5_anti_htn ~ w1.GE_male_std*dx_htn5 +', controls),
                              subset(ahdsgnbio, complete_case == 1),
                             'w1.GE_male_std', 'w5_anti_htn', 'dx_htn5')

htn.3b <- ame_by_dx_analysis(paste('w5_anti_htn ~ w4.GE_male_std*dx_htn5 +', controls),
                             subset(ahdsgnbio, complete_case == 1),
                             'w4.GE_male_std', 'w5_anti_htn', 'dx_htn5')

htn.3c <- ame_by_dx_analysis(paste('w5_anti_htn ~ delta_w1_w4_GE*dx_htn5 + w1.GE_male_std +', controls),
                             subset(ahdsgnbio, complete_case == 1),
                             'delta_w1_w4_GE', 'w5_anti_htn', 'dx_htn5')

htn.3 <- rbind(htn.3a, htn.3b, htn.3c)

# diabetes ---------------------------------------------------------------------

dm.3a <- ame_by_dx_analysis(paste('w5_anti_dm_med_use ~ w1.GE_male_std*dx_dm5 +', controls),
                            subset(ahdsgnbio, complete_case == 1),
                            'w1.GE_male_std', 'w5_anti_dm_med_use', 'dx_dm5')
                  
dm.3b <- ame_by_dx_analysis(paste('w5_anti_dm_med_use ~ w4.GE_male_std*dx_dm5 +', controls),
                            subset(ahdsgnbio, complete_case == 1),
                            'w4.GE_male_std', 'w5_anti_dm_med_use', 'dx_dm5')

dm.3c <- ame_by_dx_analysis(paste('w5_anti_dm_med_use ~ delta_w1_w4_GE*dx_dm5 + w1.GE_male_std +', controls),
                            subset(ahdsgnbio, complete_case == 1),
                            'delta_w1_w4_GE', 'w5_anti_dm_med_use', 'dx_dm5')

dm.3 <- rbind(dm.3a, dm.3b, dm.3c)


# hyperlipidemia ---------------------------------------------------------------

hld.3a <- ame_by_dx_analysis(paste('w5_anti_hld_med_use ~ w1.GE_male_std*dx_hld5 +', controls),
                             subset(ahdsgnbio, complete_case == 1),
                             'w1.GE_male_std', 'w5_anti_hld_med_use', 'dx_hld5')

hld.3b <- ame_by_dx_analysis(paste('w5_anti_hld_med_use ~ w4.GE_male_std*dx_hld5 +', controls),
                             subset(ahdsgnbio, complete_case == 1),
                             'w4.GE_male_std', 'w5_anti_hld_med_use', 'dx_hld5')

hld.3c <- ame_by_dx_analysis(paste('w5_anti_hld_med_use ~ delta_w1_w4_GE*dx_hld5 + w1.GE_male_std +', controls),
                             subset(ahdsgnbio, complete_case == 1),
                             'delta_w1_w4_GE', 'w5_anti_hld_med_use', 'dx_hld5')

hld.3 <- rbind(hld.3a, hld.3b, hld.3c)

# combine ----------------------------------------------------------------------
model.3 <- rbind(htn.3, dm.3, hld.3) %>%
  mutate(stars = case_when(
    `AME P_value` < 0.01 ~ "***",
    `AME P_value` < 0.05 ~ "**",
    `AME P_value` < 0.1 ~ "*",
    TRUE ~ ""),
    values = paste(paste0(AME,stars), `AME 95% CI`, 'P Value =', `AME P_value`)) %>%
  select(Outcome, Term, values) %>%
  pivot_wider(names_from = Term, values_from = values) %>%
  mutate(Outcome = case_when(
    Outcome == "w5_anti_htn" ~ "Tx of hypertension",
    Outcome == "w5_anti_dm_med_use"  ~ "Tx of diabetes",
    Outcome == "w5_anti_hld_med_use" ~ "Tx of hyperlipidemia"
  ),
  Model = 'Model 3') %>%
  rename(`Adolescent MGE` = w1.GE_male_std, 
         `Young adult MGE` = w4.GE_male_std,
         `Change in MGE` = delta_w1_w4_GE) 

print(model.3)

################################################################################
# Model 4 ----------------------------------------------------------------------
### (TX ~ MGE##DX + Bio-M + controls)  
################################################################################

ahdsgnbp <- svydesign(
  id=~psuscid,
  strata=~region,
  weights=~weights.bio,
  data=subset(final.df, in_sample.bio == 1 & is.na(psuscid) == F & is.na(region) == F & is.na(w5_bp) == F),
  nest=TRUE)

# hypertension -----------------------------------------------------------------

htn.4a <- ame_by_dx_analysis(paste('w5_anti_htn ~ w5_bp + w1.GE_male_std*dx_htn5 +', controls),
                             subset(ahdsgnbp, complete_case == 1), 
                             'w1.GE_male_std', 'w5_anti_htn', 'dx_htn5')

htn.4b <- ame_by_dx_analysis(paste('w5_anti_htn ~ w5_bp + w4.GE_male_std*dx_htn5 +', controls),
                             subset(ahdsgnbp, complete_case == 1),
                             'w4.GE_male_std', 'w5_anti_htn', 'dx_htn5')

htn.4c <- ame_by_dx_analysis(paste('w5_anti_htn ~ w5_bp + delta_w1_w4_GE*dx_htn5 + w1.GE_male_std +', controls),
                             subset(ahdsgnbp, complete_case == 1),
                             'delta_w1_w4_GE', 'w5_anti_htn', 'dx_htn5')

htn.4 <- rbind(htn.4a, htn.4b, htn.4c)

# diabetes ---------------------------------------------------------------------

ahdsgndm <- svydesign(
  id=~psuscid,
  strata=~region,
  weights=~weights.bio,
  data=subset(final.df, in_sample.bio == 1 & is.na(psuscid) == F & is.na(region) == F & is.na(w5_a1c) == F),
  nest=TRUE)


dm.4a <- ame_by_dx_analysis(paste('w5_anti_dm_med_use ~ w5_a1c + w1.GE_male_std*dx_dm5 +', controls),
                            subset(ahdsgndm, complete_case == 1),
                            'w1.GE_male_std', 'w5_anti_dm_med_use', 'dx_dm5')

dm.4b <- ame_by_dx_analysis(paste('w5_anti_dm_med_use ~ w5_a1c + w4.GE_male_std*dx_dm5 +', controls),
                            subset(ahdsgndm, complete_case == 1),
                            'w4.GE_male_std', 'w5_anti_dm_med_use', 'dx_dm5')

dm.4c <- ame_by_dx_analysis(paste('w5_anti_dm_med_use ~ w5_a1c + delta_w1_w4_GE*dx_dm5 + w1.GE_male_std +', controls),
                            subset(ahdsgndm, complete_case == 1),
                            'delta_w1_w4_GE', 'w5_anti_dm_med_use', 'dx_dm5')

dm.4 <- rbind(dm.4a, dm.4b, dm.4c)


# hyperlipidemia ---------------------------------------------------------------

ahdsgnhld <- svydesign(
  id=~psuscid,
  strata=~region,
  weights=~weights.bio,
  data=subset(final.df, in_sample.bio == 1 & is.na(psuscid) == F & is.na(region) == F & is.na(w5_nhdl) == F),
  nest=TRUE)

hld.4a <- ame_by_dx_analysis(paste('w5_anti_hld_med_use ~ w5_nhdl + w1.GE_male_std*dx_hld5 +', controls),
                             subset(ahdsgnhld, complete_case == 1), 
                             'w1.GE_male_std', 'w5_anti_hld_med_use', 'dx_hld5')

hld.4b <- ame_by_dx_analysis(paste('w5_anti_hld_med_use ~ w5_nhdl + w4.GE_male_std*dx_hld5 +', controls),
                             subset(ahdsgnhld, complete_case == 1), 
                             'w4.GE_male_std', 'w5_anti_hld_med_use', 'dx_hld5')

hld.4c <- ame_by_dx_analysis(paste('w5_anti_hld_med_use ~ w5_nhdl + delta_w1_w4_GE*dx_hld5 + w1.GE_male_std +', controls),
                             subset(ahdsgnhld, complete_case == 1), 
                             'delta_w1_w4_GE', 'w5_anti_hld_med_use', 'dx_hld5')

hld.4 <- rbind(hld.4a, hld.4b, hld.4c)


model.4 <- rbind(htn.4, dm.4, hld.4) %>%
  mutate(stars = case_when(
    `AME P_value` < 0.01 ~ "***",
    `AME P_value` < 0.05 ~ "**",
    `AME P_value` < 0.1 ~ "*",
    TRUE ~ ""),
    values = paste(paste0(AME,stars), `AME 95% CI`, 'P Value =', `AME P_value`)) %>%
  select(Outcome, Term, values) %>%
  pivot_wider(names_from = Term, values_from = values) %>%
  mutate(Outcome = case_when(
    Outcome == "w5_anti_htn" ~ "Tx of hypertension",
    Outcome == "w5_anti_dm_med_use"  ~ "Tx of diabetes",
    Outcome == "w5_anti_hld_med_use" ~ "Tx of hyperlipidemia"
  ),
  Model = 'Model 4') %>%
  rename(`Adolescent MGE` = w1.GE_male_std, 
         `Young adult MGE` = w4.GE_male_std,
         `Change in MGE` = delta_w1_w4_GE) 

print(model.4)

################################################################################
# Model 5 ----------------------------------------------------------------------
### (Bio-M ~ MGE + controls)
################################################################################

ahdsgn <- svydesign(
  id=~psuscid,
  strata=~region,
  weights=~w5biowgt,
  data=subset(
    final.df, in_sample.bio == 1 & is.na(psuscid) == F & is.na(region) == F
  ),
  nest=TRUE)

# hypertension -----------------------------------------------------------------

ahdsgnbp <- svydesign(
  id=~psuscid,
  strata=~region,
  weights=~weights.bio,
  data=subset(final.df, in_sample.bio == 1 & is.na(psuscid) == F & is.na(region) == F & is.na(w5_anti_htn) == F),
  nest=TRUE)


htn.5a <- ame_analysis(paste('w5_bp ~ w1.GE_male_std + w5_anti_htn +', controls),
                       subset(ahdsgnbp, complete_case == 1),
                       'w1.GE_male_std', 'w5_bp')

htn.5b <- ame_analysis(paste('w5_bp ~ w4.GE_male_std + w5_anti_htn +', controls),
                       subset(ahdsgnbp, complete_case == 1),
                       'w4.GE_male_std', 'w5_bp')

htn.5c <- ame_analysis(paste('w5_bp ~ delta_w1_w4_GE + w1.GE_male_std + w5_anti_htn +', controls),
                       subset(ahdsgnbp, complete_case == 1),
                       'delta_w1_w4_GE', 'w5_bp')

htn.5 <- rbind(htn.5a, htn.5b, htn.5c)

# diabetes ---------------------------------------------------------------------

ahdsgndm <- svydesign(
  id=~psuscid,
  strata=~region,
  weights=~weights.bio,
  data=subset(final.df, in_sample.bio == 1 & is.na(psuscid) == F & is.na(region) == F & is.na(w5_anti_dm_med_use) == F),
  nest=TRUE)

dm.5a <- ame_analysis(paste('w5_a1c ~ w1.GE_male_std + w5_anti_dm_med_use +', controls),
                      subset(ahdsgndm, complete_case == 1),
                      'w1.GE_male_std', 'w5_a1c')

dm.5b <- ame_analysis(paste('w5_a1c ~ w4.GE_male_std + w5_anti_dm_med_use +', controls),
                      subset(ahdsgndm, complete_case == 1),
                      'w4.GE_male_std', 'w5_a1c')

dm.5c <- ame_analysis(paste('w5_a1c ~ delta_w1_w4_GE + w1.GE_male_std + w5_anti_dm_med_use +', controls),
                      subset(ahdsgndm, complete_case == 1),
                      'delta_w1_w4_GE', 'w5_a1c')

dm.5 <- rbind(dm.5a, dm.5b, dm.5c)

# hyperlipidemia ---------------------------------------------------------------

ahdsgnhld <- svydesign(
  id=~psuscid,
  strata=~region,
  weights=~weights.bio,
  data=subset(final.df, in_sample.bio == 1 & is.na(psuscid) == F & is.na(region) == F & is.na(w5_anti_hld_med_use) == F),
  nest=TRUE)

hld.5a <- ame_analysis(paste('w5_nhdl ~ w1.GE_male_std + w5_anti_hld_med_use +', controls),
                       subset(ahdsgnhld, complete_case == 1),
                       'w1.GE_male_std', 'w5_nhdl')

hld.5b <- ame_analysis(paste('w5_nhdl ~ w4.GE_male_std + w5_anti_hld_med_use +', controls),
                       subset(ahdsgnhld, complete_case == 1),
                       'w4.GE_male_std', 'w5_nhdl')

hld.5c <- ame_analysis(paste('w5_nhdl ~ delta_w1_w4_GE + w1.GE_male_std + w5_anti_hld_med_use +', controls),
                       subset(ahdsgnhld, complete_case == 1),
                       'delta_w1_w4_GE', 'w5_nhdl')

hld.5 <- rbind(hld.5a, hld.5b, hld.5c)

model.5 <- rbind(htn.5, dm.5, hld.5) %>%
  mutate(stars = case_when(
    `AME P_value` < 0.01 ~ "***",
    `AME P_value` < 0.05 ~ "**",
    `AME P_value` < 0.1 ~ "*",
    TRUE ~ ""),
    values = paste(paste0(AME,stars), `AME 95% CI`, 'P Value =', `AME P_value`)) %>%
  select(Outcome, Term, values) %>%
  pivot_wider(names_from = Term, values_from = values) %>%
  mutate(Outcome = case_when(
    Outcome == "w5_bp" ~ "Bio-M of hypertension",
    Outcome == "w5_a1c"  ~ "Bio-M of diabetes",
    Outcome == "w5_nhdl" ~ "Bio-M of hyperlipidemia"
  ),
  Model = 'Model 5') %>%
  rename(`Adolescent MGE` = w1.GE_male_std, 
         `Young adult MGE` = w4.GE_male_std,
         `Change in MGE` = delta_w1_w4_GE) 

print(model.5)

################################################################################
# Create Table 2 --------------------------------------------------------------
################################################################################

table2 <- rbind(model.1, model.2, model.3, model.4, model.5)

table2 <- table2 %>%
kable(caption = "Table 2. Coefficients (dy/dx) Estimating Associations 
      between MGE and Adult Awareness, Treatment, 
      and Biomeasure Evidence of Adult CVD Risks",
      col.names = c("Outcome", 
                    "Adolescent MGE", "Young Adult MGE", "Change in MGE", "Model")) %>%
  kable_styling(font_size = 12, full_width = F) %>%
  add_header_above(c(" " = 5)) 

save_kable(table2, file = "outputs/tables/Table 2 with Imputation All Missing.pdf")


