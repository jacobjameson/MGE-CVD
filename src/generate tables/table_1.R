#-------------------------------------------------------------------------
# AUTHOR:             Jacob Jameson
# PURPOSE:            Generate table 1
#-------------------------------------------------------------------------
#
# load packages ----------------------------------------------------------
rm(list = ls())

libs <- c("tidyverse", "haven", 'scales', 'survey', 'gridExtra',
          'gtsummary', 'labelled', 'broom', 'ggpubr', 'sandwich',
          'lmtest', 'marginaleffects', 'webshot2', 'gt',
          'kableExtra')

installed_libs <- libs %in% rownames (installed.packages ())
if (any (installed_libs == F)) {
  install.packages (libs[!installed_libs])
}

invisible(lapply (libs, library, character.only = T))

theme_gtsummary_journal(journal = "jama")
theme_gtsummary_compact()

source('src/prepare data/Construct Analytical Dataset.R')
#-------------------------------------------------------------------------
# Table 1 ----------------------------------------------------------------

# Unweighted -------------------------------------------------------------

t1.vars <- c("race5", 'edu5', 'ins5', 'sespc_al', 'nhood1_d', 
             "dx_htn5", "w5_bp", "w5_anti_htn",
             "dx_hld5", "w5_anti_hld_med_use",
             "lipid5cat", "w5_anti_htn", "dx_dm5", "w5_a1c", 
             "w5_anti_dm_med_use")

final.df$increasing <- ifelse(final.df$w4.GE_male_std > 0, 
                              'Wave IV GE Above Avg', 'Wave IV GE Below Avg')

sum <- final.df %>% 
  filter(in_sample.5 == 1) %>%
  select(one_of(t1.vars), increasing) %>%
  tbl_summary(by = increasing,
              type = all_continuous() ~ "continuous2",
              statistic = all_continuous() ~ c("{mean} ({sd})"),
              missing_text = "(Missing)") %>% 
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 2),
        test = list(all_continuous() ~ "aov",
                    all_categorical() ~ "chisq.test")) %>%
  add_overall() %>% modify_header(label ~ "**Variable**") %>%
  modify_caption("**Patient Characteristics (unweighted)**") %>%
  bold_labels() %>% as_gt()

print(sum)
