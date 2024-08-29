################################################################################
# AUTHOR:             Jacob Jameson
#
# DESCRIPTION:        This script generates the figures for the paper.
#
# DEPENDENCIES:       The script requires the following dataset:
#                     -final.df, which is the final dataset constructed by 
#                      executing scr/prepare data/Construct Analytical dataset.R
################################################################################

final.df$complete_case <- ifelse(
  rowSums(
    is.na(
      final.df[, c('race5', 'sespc_al', 'nhood1_d', 'ins5', 'edu5')]
    )
  ) == 0, 1, 0)


ahdsgnbio <- svydesign(
  id=~psuscid,
  strata=~region,
  weights=~weights.bio,
  data=subset(final.df, in_sample.bio == 1 & is.na(psuscid) == F & is.na(region) == F),
  nest=TRUE)

# hypertension -----------------------------------------------------------------
htn.2a <- svyglm(dx_htn5 ~ w1.GE_male_std*w5_bp + 
                 race5 + sespc_al + nhood1_d + ins5 + edu5,
                 subset(ahdsgnbio, complete_case == 1),
                 family = 'quasibinomial')

htn.2a <- ggpredict(htn.2a, terms = c('w1.GE_male_std [-2:2, by=0.5]', "w5_bp [0,1]")) %>%
  as.data.frame() %>%
  select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high, group) %>%
  filter(group == 1) %>%
  mutate(model_type = 'w1.GE_male_std')

htn.2b <- svyglm(dx_htn5 ~ w4.GE_male_std*w5_bp + 
                   race5 + sespc_al + nhood1_d + ins5 + edu5,
                 subset(ahdsgnbio, complete_case == 1),
                 family = 'quasibinomial')

htn.2b <- ggpredict(htn.2b, terms = c('w4.GE_male_std [-2:2, by=0.5]', "w5_bp [0,1]")) %>%
  as.data.frame() %>%
  select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high, group) %>%
  filter(group == 1)  %>%
  mutate(model_type = 'w4.GE_male_std')


htn.2c <- svyglm(dx_htn5 ~ delta_w1_w4_GE*w5_bp + w1.GE_male_std +
                   race5 + sespc_al + nhood1_d + ins5 + edu5,
                 subset(ahdsgnbio, complete_case == 1),
                 family = 'quasibinomial')

htn.2c <- ggpredict(htn.2c, terms = c('delta_w1_w4_GE [-2:2, by=0.5]', "w5_bp [0,1]")) %>%
  as.data.frame() %>%
  select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high, group) %>%
  filter(group == 1) %>%
  mutate(model_type = 'delta_w1_w4_GE')

htn.2 <- rbind(htn.2a, htn.2b, htn.2c)
htn.2$model_outcome <- 'dx_htn5'
# diabetes ---------------------------------------------------------------------

dm.2a <- svyglm(dx_dm5 ~ w1.GE_male_std*w5_a1c + 
                   race5 + sespc_al + nhood1_d + ins5 + edu5,
                 subset(ahdsgnbio, complete_case == 1),
                 family = 'quasibinomial')

dm.2a <- ggpredict(dm.2a, terms = c('w1.GE_male_std [-2:2, by=0.5]', "w5_a1c [0,1]")) %>%
  as.data.frame() %>%
  select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high, group) %>%
  filter(group == 1) %>%
  mutate(model_type = 'w1.GE_male_std')


dm.2b <- svyglm(dx_dm5 ~ w4.GE_male_std*w5_a1c + 
                   race5 + sespc_al + nhood1_d + ins5 + edu5,
                 subset(ahdsgnbio, complete_case == 1),
                 family = 'quasibinomial')

dm.2b <- ggpredict(dm.2b, terms = c('w4.GE_male_std [-2:2, by=0.5]', "w5_a1c [0,1]")) %>%
  as.data.frame() %>%
  select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high, group) %>%
  filter(group == 1)  %>%
  mutate(model_type = 'w4.GE_male_std')


dm.2c <- svyglm(dx_dm5 ~ delta_w1_w4_GE*w5_a1c + w1.GE_male_std +
                   race5 + sespc_al + nhood1_d + ins5 + edu5,
                 subset(ahdsgnbio, complete_case == 1),
                 family = 'quasibinomial')

dm.2c <- ggpredict(dm.2c, terms = c('delta_w1_w4_GE [-2:2, by=0.5]', "w5_a1c [0,1]")) %>%
  as.data.frame() %>%
  select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high, group) %>%
  filter(group == 1)  %>%
  mutate(model_type = 'delta_w1_w4_GE')

dm.2 <- rbind(dm.2a, dm.2b, dm.2c)
dm.2$model_outcome <- 'dx_dm5'

# hyperlipidemia --------------------------------------------------------------

hld.2a <- svyglm(dx_hld5 ~ w1.GE_male_std*w5_nhdl + 
                  race5 + sespc_al + nhood1_d + ins5 + edu5,
                subset(ahdsgnbio, complete_case == 1),
                family = 'quasibinomial')

hld.2a <- ggpredict(hld.2a, terms = c('w1.GE_male_std [-2:2, by=0.5]', "w5_nhdl [0,1]")) %>%
  as.data.frame() %>%
  select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high, group) %>%
  filter(group == 1) %>%
  mutate(model_type = 'w1.GE_male_std')

hld.2b <- svyglm(dx_hld5 ~ w4.GE_male_std*w5_nhdl + 
                  race5 + sespc_al + nhood1_d + ins5 + edu5,
                subset(ahdsgnbio, complete_case == 1),
                family = 'quasibinomial')

hld.2b <- ggpredict(hld.2b, terms = c('w4.GE_male_std [-2:2, by=0.5]', "w5_nhdl [0,1]")) %>%
  as.data.frame() %>%
  select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high, group) %>%
  filter(group == 1) %>%
  mutate(model_type = 'w4.GE_male_std')


hld.2c <- svyglm(dx_hld5 ~ delta_w1_w4_GE*w5_nhdl + w1.GE_male_std +
                  race5 + sespc_al + nhood1_d + ins5 + edu5,
                subset(ahdsgnbio, complete_case == 1),
                family = 'quasibinomial')

hld.2c <- ggpredict(hld.2c, terms = c('delta_w1_w4_GE [-2:2, by=0.5]', "w5_nhdl [0,1]")) %>%
  as.data.frame() %>%
  select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high, group) %>%
  filter(group == 1) %>%
  mutate(model_type = 'delta_w1_w4_GE')

hld.2 <- rbind(hld.2a, hld.2b, hld.2c)
hld.2$model_outcome <- 'dx_hld5'

# Combine all plot data frames into one ----------------------------------------
combined_plot_data <- rbind(htn.2, dm.2, hld.2)

# rename the model type if w1.GE_male_std then Adolescent MGE
combined_plot_data <- combined_plot_data %>%
  mutate(model_type = case_when(
    model_type == 'w1.GE_male_std' ~ 'Adolescent MGE',
    model_type == 'w4.GE_male_std' ~ 'Young Adult MGE',
    model_type == 'delta_w1_w4_GE' ~ 'Adolescent-to-Young Adult Change in MGE'
  ),
  model_outcome = case_when(
    model_outcome == 'dx_htn5' ~ 'Hypertension Diagnosis',
    model_outcome == 'dx_dm5' ~ 'Diabetes Diagnosis',
    model_outcome == 'dx_hld5' ~ 'Hyperlipidemia Diagnosis'
  ))

combined_plot_data$model_outcome <- factor(
  combined_plot_data$model_outcome, levels = c('Hypertension Diagnosis',
                                               'Diabetes Diagnosis', 
                                               'Hyperlipidemia Diagnosis'))

combined_plot_data$model_type <- factor(
  combined_plot_data$model_type, 
  levels = c('Adolescent MGE', 'Young Adult MGE', 'Adolescent-to-Young Adult Change in MGE'))



ggplot(combined_plot_data, aes(x=xvals)) +
  geom_ribbon(aes(ymin=lower,ymax=upper, fill= model_type), alpha=0.45) +
  geom_line(aes(y=coef), size=1) + 
  facet_grid(~model_type ~ model_outcome) +
  geom_vline(xintercept =0, color='darkred', linetype = 'dotted') +
  scale_fill_manual(values = c('Adolescent MGE' = '#DF8F44FF', 
                               'Young Adult MGE' = '#00A1D5FF', 
                               'Adolescent-to-Young Adult Change in MGE' = '#B24745FF')) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = c(0,0.25,0.5,0.75,1)) +
  labs(fill = "",
       y = "Predicted Probability\n", 
       x='\nStandardized Male Gender Expressivity\n',
       title = 'Figure 1: Associations of Male Gender Expressivity (MGE) with \nAdult Diagnosis of Hypertension, Diabetes, and Hyperlipidemia\n',
       caption = str_wrap('\n\nCaption: Higher male gender expressivity (MGE) in adolescence 
                          was significantly associated with a lower predicted probability of adult  
                          diabetes diagnosis among men with elevated hemoglobin a1c levels ≥ 6.5%. 
                          Higher MGE in young adulthood was significantly
                          associated with a lower predicted probability of adult hypertension diagnosis among 
                          men with elevated blood pressure (≥ 130 mm Hg systolic 
                          and/or ≥ 80 mm Hg diastolic).', 125)) +
  theme_calc(base_family = "Arial") +
  theme(plot.title = element_text(size = rel(1.2)),
        text = element_text(size = 20),
        strip.text.y = element_blank(),
        panel.background = element_rect(fill = 'white'),
        plot.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = NA),
        axis.title = element_text(size = rel(1)),
        axis.title.y = element_text(angle=90,vjust =2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(), 
        plot.caption = element_text(color = "black", size = 15, hjust = 0),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(),
        panel.spacing = unit(1, "lines"),
        panel.grid.major = element_line(colour="grey90"),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(colour = NA),
        legend.background = element_rect(fill = 'white'),
        legend.position = "top",
        legend.title.position = 'top',
        legend.direction = "horizontal",
        legend.key.size= unit(0.8, "cm"),
        legend.title = element_text(face="italic", hjust = 0.5),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="black",fill="#374E55FF"),
        strip.text = element_text(face="bold", color="white", size=rel(1))
  )

# save plot with white background
ggsave('outputs/figures/figure 1.png', width = 13.5, height = 11)
ggsave('outputs/figures/figure 1.svg', width = 13.5, height = 11)

# hypertension -----------------------------------------------------------------
htn.3a <- svyglm(w5_anti_htn ~ w1.GE_male_std*dx_htn5 + 
                   race5 + sespc_al + nhood1_d + ins5 + edu5,
                 subset(ahdsgnbio, complete_case == 1),
                 family = 'quasibinomial')

htn.3a <- ggpredict(htn.3a, terms = c('w1.GE_male_std [-2:2, by=0.5]', "dx_htn5 [0,1]")) %>%
  as.data.frame() %>%
  select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high, group) %>%
  filter(group == 1) %>%
  mutate(model_type = 'w1.GE_male_std')

htn.3b <- svyglm(w5_anti_htn ~ w4.GE_male_std*dx_htn5 + 
                   race5 + sespc_al + nhood1_d + ins5 + edu5,
                 subset(ahdsgnbio, complete_case == 1),
                 family = 'quasibinomial')

htn.3b <- ggpredict(htn.3b, terms = c('w4.GE_male_std [-2:2, by=0.5]', "dx_htn5 [0,1]")) %>%
  as.data.frame() %>%
  select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high, group) %>%
  filter(group == 1)  %>%
  mutate(model_type = 'w4.GE_male_std')


htn.3c <- svyglm(w5_anti_htn ~ delta_w1_w4_GE*dx_htn5 + w1.GE_male_std +
                   race5 + sespc_al + nhood1_d + ins5 + edu5,
                 subset(ahdsgnbio, complete_case == 1),
                 family = 'quasibinomial')

htn.3c <- ggpredict(htn.3c, terms = c('delta_w1_w4_GE [-2:2, by=0.5]', "dx_htn5 [0,1]")) %>%
  as.data.frame() %>%
  select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high, group) %>%
  filter(group == 1) %>%
  mutate(model_type = 'delta_w1_w4_GE')

htn.3 <- rbind(htn.3a, htn.3b, htn.3c)
htn.3$model_outcome <- 'w5_anti_htn'
# diabetes ---------------------------------------------------------------------

dm.3a <- svyglm(w5_anti_dm_med_use ~ w1.GE_male_std*dx_dm5 + 
                  race5 + sespc_al + nhood1_d + ins5 + edu5,
                subset(ahdsgnbio, complete_case == 1),
                family = 'quasibinomial')

dm.3a <- ggpredict(dm.3a, terms = c('w1.GE_male_std [-2:2, by=0.5]', "dx_dm5 [0,1]")) %>%
  as.data.frame() %>%
  select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high, group) %>%
  filter(group == 1) %>%
  mutate(model_type = 'w1.GE_male_std')


dm.3b <- svyglm(w5_anti_dm_med_use ~ w4.GE_male_std*dx_dm5 + 
                  race5 + sespc_al + nhood1_d + ins5 + edu5,
                subset(ahdsgnbio, complete_case == 1),
                family = 'quasibinomial')

dm.3b <- ggpredict(dm.3b, terms = c('w4.GE_male_std [-2:2, by=0.5]', "dx_dm5 [0,1]")) %>%
  as.data.frame() %>%
  select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high, group) %>%
  filter(group == 1)  %>%
  mutate(model_type = 'w4.GE_male_std')


dm.3c <- svyglm(w5_anti_dm_med_use ~ delta_w1_w4_GE*dx_dm5 + w1.GE_male_std +
                  race5 + sespc_al + nhood1_d + ins5 + edu5,
                subset(ahdsgnbio, complete_case == 1),
                family = 'quasibinomial')

dm.3c <- ggpredict(dm.3c, terms = c('delta_w1_w4_GE [-2:2, by=0.5]', "dx_dm5 [0,1]")) %>%
  as.data.frame() %>%
  select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high, group) %>%
  filter(group == 1)  %>%
  mutate(model_type = 'delta_w1_w4_GE')

dm.3 <- rbind(dm.3a, dm.3b, dm.3c)
dm.3$model_outcome <- 'w5_anti_dm_med_use'

# hyperlipidemia --------------------------------------------------------------

hld.3a <- svyglm(w5_anti_hld_med_use ~ w1.GE_male_std*dx_hld5 + 
                   race5 + sespc_al + nhood1_d + ins5 + edu5,
                 subset(ahdsgnbio, complete_case == 1),
                 family = 'quasibinomial')

hld.3a <- ggpredict(hld.3a, terms = c('w1.GE_male_std [-2:2, by=0.5]', "dx_hld5 [0,1]")) %>%
  as.data.frame() %>%
  select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high, group) %>%
  filter(group == 1) %>%
  mutate(model_type = 'w1.GE_male_std')

hld.3b <- svyglm(w5_anti_hld_med_use ~ w4.GE_male_std*dx_hld5 + 
                   race5 + sespc_al + nhood1_d + ins5 + edu5,
                 subset(ahdsgnbio, complete_case == 1),
                 family = 'quasibinomial')

hld.3b <- ggpredict(hld.3b, terms = c('w4.GE_male_std [-2:2, by=0.5]', "dx_hld5 [0,1]")) %>%
  as.data.frame() %>%
  select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high, group) %>%
  filter(group == 1) %>%
  mutate(model_type = 'w4.GE_male_std')


hld.3c <- svyglm(w5_anti_hld_med_use ~ delta_w1_w4_GE*dx_hld5 + w1.GE_male_std +
                   race5 + sespc_al + nhood1_d + ins5 + edu5,
                 subset(ahdsgnbio, complete_case == 1),
                 family = 'quasibinomial')

hld.3c <- ggpredict(hld.3c, terms = c('delta_w1_w4_GE [-2:2, by=0.5]', "dx_hld5 [0,1]")) %>%
  as.data.frame() %>%
  select(xvals = x, coef = predicted, lower = conf.low, upper = conf.high, group) %>%
  filter(group == 1) %>%
  mutate(model_type = 'delta_w1_w4_GE')

hld.3 <- rbind(hld.3a, hld.3b, hld.3c)
hld.3$model_outcome <- 'w5_anti_hld_med_use'

# Combine all plot data frames into one ----------------------------------------
combined_plot_data <- rbind(htn.3, dm.3, hld.3)

# rename the model type if w1.GE_male_std then Adolescent MGE
combined_plot_data <- combined_plot_data %>%
  mutate(model_type = case_when(
    model_type == 'w1.GE_male_std' ~ 'Adolescent MGE',
    model_type == 'w4.GE_male_std' ~ 'Young Adult MGE',
    model_type == 'delta_w1_w4_GE' ~ 'Adolescent-to-Young Adult Change in MGE'
  ),
  model_outcome = case_when(
    model_outcome == 'w5_anti_htn' ~ 'Hypertension Treatment',
    model_outcome == 'w5_anti_dm_med_use' ~ 'Diabetes Treatment',
    model_outcome == 'w5_anti_hld_med_use' ~ 'Hyperlipidemia Treatment'
  ))

combined_plot_data$model_outcome <- factor(
  combined_plot_data$model_outcome, levels = c('Hypertension Treatment',
                                               'Diabetes Treatment', 
                                               'Hyperlipidemia Treatment'))

combined_plot_data$model_type <- factor(
  combined_plot_data$model_type, 
  levels = c('Adolescent MGE', 'Young Adult MGE', 
             'Adolescent-to-Young Adult Change in MGE'))




ggplot(combined_plot_data, aes(x=xvals)) +
  geom_ribbon(aes(ymin=lower,ymax=upper, fill= model_type), alpha=0.45) +
  geom_line(aes(y=coef), size=1) + 
  facet_grid(~model_type ~ model_outcome) +
  geom_vline(xintercept =0, color='darkred', linetype = 'dotted') +
  scale_fill_manual(values = c('Adolescent MGE' = '#DF8F44FF', 
                               'Young Adult MGE' = '#00A1D5FF', 
                               'Adolescent-to-Young Adult Change in MGE' = '#B24745FF')) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = c(0,0.25,0.5,0.75,1)) +
  labs(fill = "",
       y = "Predicted Probability\n", 
       x='\nStandardized Male Gender Expressivity\n',
       title = 'Figure 2: Associations of Male Gender Expressivity (MGE) with \nAdult Treatment of Hypertension, Diabetes, and Hyperlipidemia\n',
       caption = str_wrap('\n\n\n Caption: Higher male gender expressivity (MGE) in adolescence 
                          was significantly associated with a lower predicted probability of adult 
                          treatment of hypertension among men who reported 
                          hypertension diagnoses. Higher 
                          MGE in young adulthood was 
                          significantly associated with a lower predicted probability of adult 
                          treatment of hypertension and diabetes among men 
                          who reported diagnoses of hypertension and diabetes respectively.', 125)) +
  theme_calc(base_family = "Arial") +
  theme(plot.title = element_text(size = rel(1.2)),
        text = element_text(size = 20),
        strip.text.y = element_blank(),
        panel.background = element_rect(fill = 'white'),
        plot.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = NA),
        axis.title = element_text(size = rel(1)),
        axis.title.y = element_text(angle=90,vjust =2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(), 
        plot.caption = element_text(color = "black", size = 15, hjust = 0),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(),
        panel.spacing = unit(1, "lines"),
        panel.grid.major = element_line(colour="grey90"),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(colour = NA),
        legend.background = element_rect(fill = 'white'),
        legend.position = "top",
        legend.title.position = 'top',
        legend.direction = "horizontal",
        legend.key.size= unit(0.8, "cm"),
        legend.title = element_text(face="italic", hjust = 0.5),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="black",fill="#374E55FF"),
        strip.text = element_text(face="bold", color="white", size=rel(1))
  )



ggsave('outputs/figures/figure 2.png', width = 13.5, height = 11)
ggsave('outputs/figures/figure 2.svg', width = 13.5, height = 11)

