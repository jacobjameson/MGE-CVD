################################################################################
# AUTHOR:             Jacob Jameson
#
# DESCRIPTION:        This script generates the figures for the paper.
#
# DEPENDENCIES:       The script requires the following dataset:
#                     -final.df, which is the final dataset constructed by 
#                      executing scr/prepare data/Construct Analytical dataset.R
################################################################################

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


plot_data_frames <- list() 

for (FORMULA in c(formulas_htn, formulas_dm, formulas_hld)) {
  model <- glm(FORMULA,
               data = subset_df, 
               weights = subset_weights, 
               family = "quasibinomial")
  
  vcv_robust <- vcovHC(model, type = "HC0", cluster = ~cluster)
  
  term.1 <- str_squish(word(FORMULA,3,sep = " "))
  term.2 <- str_squish(word(FORMULA,5,sep = " "))
  

  plot_df <- ggpredict(model,
                       vcov.type = "HC0", 
                       vcov.args = list(cluster = subset_df$cluster),
                       terms = c(paste(term.1, '[-2:2, by=0.1]'), paste(term.2, " [0,1]"))) %>%
    as.data.frame() %>%
    select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high, group) %>%
    filter(group == 1)
  
  # Append the model identifier for clarity in the combined plot data frame
  plot_df$model_type <- term.1
  plot_df$model_outcome <- str_squish(word(FORMULA,1,sep = "\\~"))
  
  plot_data_frames[[length(plot_data_frames) + 1]] <- plot_df
}

# Combine all plot data frames into one
combined_plot_data <- do.call(rbind, plot_data_frames)

# rename the model type if w1.GE_male_std then Adolescent MGE
combined_plot_data <- combined_plot_data %>%
  mutate(model_type = case_when(
    model_type == 'w1.GE_male_std' ~ 'Adolescent MGE',
    model_type == 'w4.GE_male_std' ~ 'Young Adult MGE',
    model_type == 'delta_w1_w4_GE' ~ 'Adolescent-to-Young Adult \nChange in MGE'
  ),
  model_outcome = case_when(
    model_outcome == 'dx_htn5' ~ 'Hypertension Awareness',
    model_outcome == 'dx_dm5' ~ 'Diabetes Awareness',
    model_outcome == 'dx_hld5' ~ 'Hyperlipidemia Awareness'
  ))

combined_plot_data$model_outcome <- factor(
  combined_plot_data$model_outcome, levels = c('Hypertension Awareness',
                                               'Diabetes Awareness', 
                                               'Hyperlipidemia Awareness'))

combined_plot_data$model_type <- factor(
  combined_plot_data$model_type, 
  levels = c('Adolescent MGE', 'Young Adult MGE', 'Adolescent-to-Young Adult \nChange in MGE'))


ggplot(combined_plot_data, aes(x=xvals)) +
  geom_ribbon(aes(ymin=lower,ymax=upper, fill= model_outcome), alpha=0.65) +
  geom_line(aes(y=coef)) + 
  facet_grid(~model_outcome ~ model_type) +
  scale_fill_manual(values = 
                      c('Diabetes Awareness' = '#6C3E93', 
                        'Hypertension Awareness' = '#C04848', 
                        'Hyperlipidemia Awareness' = '#FFC857')) +
  theme_bw(base_family = "Times New Roman") +
  geom_vline(xintercept =0, color='darkred') +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = c(0,0.25,0.5,0.75,1)) +
  labs(fill = "MGE Measure",
       y = "Predicted Probability\n", 
       x='\nStandardized Male Gender Expressivity\n',
       title = 'Associations of Male Gender Expressivity (MGE) with Predicted Probability of \nHypertension, Diabetes, or Hyperlipidemia Diagnoses Awareness in Men with \nBiomeasure Evidence of Disease\n',
       caption = str_wrap('\n\nIncreased male gender expressivity in adolescence 
                          is significantly associated with decreased adult awareness 
                          of diabetes among men with elevated hemoglobin a1c levels > 6.5%. 
                          Increased male gender expression in young adulthood is significantly
                          associated with decreased adult awareness of hypertension among 
                          men with elevated blood pressure readings (> 130 mm Hg systolic 
                          and/or > 80 mm Hg diastolic)', 125)) +
  theme_economist() +
  theme(strip.text.y = element_blank(),
      plot.background = element_rect(fill = 'white', color = NA),
        panel.grid.major = element_line(color = 'grey85', size = 0.5),
        axis.text.x = element_text(color = "black", size = 18),
        axis.text.y = element_text(color = "black", size = 18),
        axis.title.x = element_text(color = 'black',size = 20),
        panel.spacing = unit(1, "lines"),
        axis.line.x = element_line(color = 'black', size = 0.5),
        axis.line.y = element_line(color = 'black', size = 0.5),
        axis.title.y = element_text(color = 'black', size = 20),
        plot.title = element_text(color = "black", size = 20, hjust = 0, face = "bold"),
        plot.caption = element_text(color = "black", size = 15, hjust = 0),
        legend.position = 'top',
        legend.title = element_blank(),
        strip.text.x = element_text(color = 'black', size = 18, face = "bold",
                                    margin = margin(0.1,0,0.1,0, "cm")),
        legend.text = element_text(color = "black", size = 18)) 

# save plot with white background
ggsave('outputs/figures/figure.1a.png', width = 13, height = 11)

################################################################################

# Define the subset and weights
subset_df <- subset(final.df, in_sample.5 == 1)
subset_weights <- subset_df$weights

formulas_htn <- c(
  'w5_anti_htn ~ w1.GE_male_std + dx_htn5 + w1.GE_male_std*dx_htn5 + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'w5_anti_htn ~ w4.GE_male_std + dx_htn5 + w4.GE_male_std*dx_htn5 + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'w5_anti_htn ~ delta_w1_w4_GE + dx_htn5 + delta_w1_w4_GE*dx_htn5 + w1.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5'
)

formulas_dm <- c(
  'w5_anti_dm_med_use ~ w1.GE_male_std + dx_dm5 + w1.GE_male_std*dx_dm5 + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'w5_anti_dm_med_use ~ w4.GE_male_std + dx_dm5 + w4.GE_male_std*dx_dm5 + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'w5_anti_dm_med_use ~ delta_w1_w4_GE + dx_dm5 + delta_w1_w4_GE*dx_dm5 + w1.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5'
)

formulas_hld <- c(
  'w5_anti_hld_med_use ~ w1.GE_male_std + dx_hld5 + w1.GE_male_std*dx_hld5 + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'w5_anti_hld_med_use ~ w4.GE_male_std + dx_hld5 + w4.GE_male_std*dx_hld5 + race5 + sespc_al + nhood1_d + ins5 + edu5',
  'w5_anti_hld_med_use ~ delta_w1_w4_GE + dx_hld5 + delta_w1_w4_GE*dx_hld5 + w1.GE_male_std + race5 + sespc_al + nhood1_d + ins5 + edu5'
)


plot_data_frames <- list() 

for (FORMULA in c(formulas_htn, formulas_dm, formulas_hld)) {
  model <- glm(FORMULA,
               data = subset_df, 
               weights = subset_weights, 
               family = "quasibinomial")
  
  vcv_robust <- vcovHC(model, type = "HC0", cluster = ~cluster)
  
  term.1 <- str_squish(word(FORMULA,3,sep = " "))
  term.2 <- str_squish(word(FORMULA,5,sep = " "))
  
  
  plot_df <- ggpredict(model,
                       vcov.type = "HC0", 
                       vcov.args = list(cluster = subset_df$cluster),
                       terms = c(paste(term.1, '[-2:2, by=0.1]'), paste(term.2, " [0,1]"))) %>%
    as.data.frame() %>%
    select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high, group) %>%
    filter(group == 1)
  
  # Append the model identifier for clarity in the combined plot data frame
  plot_df$model_type <- term.1
  plot_df$model_outcome <- str_squish(word(FORMULA,1,sep = "\\~"))
  
  plot_data_frames[[length(plot_data_frames) + 1]] <- plot_df
}

# Combine all plot data frames into one
combined_plot_data <- do.call(rbind, plot_data_frames)

# rename the model type if w1.GE_male_std then Adolescent MGE
combined_plot_data <- combined_plot_data %>%
  mutate(model_type = case_when(
    model_type == 'w1.GE_male_std' ~ 'Adolescent MGE',
    model_type == 'w4.GE_male_std' ~ 'Young Adult MGE',
    model_type == 'delta_w1_w4_GE' ~ 'Adolescent-to-Young Adult Change in MGE'
  ),
  model_outcome = case_when(
    model_outcome == 'w5_anti_htn' ~ 'Treatment of Hypertension',
    model_outcome == 'w5_anti_dm_med_use' ~ 'Treatment of Diabetes',
    model_outcome == 'w5_anti_hld_med_use' ~ 'Treatment of Hyperlipidemia'
  ))

combined_plot_data$model_outcome <- factor(
  combined_plot_data$model_outcome, levels = c('Treatment of Hypertension',
                                               'Treatment of Diabetes', 
                                               'Treatment of Hyperlipidemia'))

combined_plot_data$model_type <- factor(
  combined_plot_data$model_type, 
  levels = c('Adolescent MGE', 'Young Adult MGE', 'Adolescent-to-Young Adult Change in MGE'))


ggplot(combined_plot_data, aes(x=xvals)) +
  geom_ribbon(aes(ymin=lower,ymax=upper, fill= model_type), alpha=0.5) +
  geom_line(aes(y=coef)) + 
  facet_grid(~model_type ~ model_outcome) +
  scale_fill_manual(values = 
                      c('Adolescent MGE' = '#DF8F44FF', 
                        'Young Adult MGE' = '#00A1D5FF', 
                        'Adolescent-to-Young Adult Change in MGE' = '#B24745FF')) +
  theme_bw(base_family = "Times New Roman") +
  geom_vline(xintercept =0, color='darkred') +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = c(0,0.25,0.5,0.75,1)) +
  labs(fill = "MGE Measure\n",
       y = "Predicted Probability\n", 
       x='\nStandardized Male Gender Expressivity\n',
       title = 'Associations of Male Gender Expressivity (MGE) with Predicted Probability of \nAdult Treatment of Hypertension, Diabetes, or Hyperlipidemia in Men who \nReport Having These Conditions. \n',
       caption = str_wrap('\n\n\nIncreased male gender expressivity in adolescence 
                          is significantly associated with decreased adult 
                          treatment of hypertension among men who report 
                          awareness of hypertension diagnoses. Increased 
                          male gender expression in young adulthood is 
                          significantly associated with decreased adult 
                          treatment of hypertension and diabetes among men 
                          who report awareness of hypertension and diabetes respectively)', 125)) +
  theme(strip.text.y = element_blank(),
        axis.text.x = element_text(color = "black", size = 18),
        axis.text.y = element_text(color = "black", size = 18),
        axis.title.x = element_text(color = 'black',size = 20),
        axis.title.y = element_text(color = 'black', size = 20),
        plot.title = element_text(color = "black", size = 20, hjust = 0, face = "bold"),
        plot.caption = element_text(color = "black", size = 15, hjust = 0),
        legend.position = 'top',
        legend.title = element_blank(),
        strip.text.x = element_text(color = 'black', size = 18, face = "bold"),
        legend.text = element_text(color = "black", size = 18))

ggsave('outputs/figures/figure.2.png', width = 12, height = 11)

