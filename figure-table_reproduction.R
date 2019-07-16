########### BOOTSTRAP VARIABLE EFFECTS PLOT

# Joining and mutating the bootstrap results
plot_boot <- mutate(bootStepImp, Model = "Stepwise") %>% 
  bind_rows(mutate(bootLassoImp, Model = "Lasso")) %>% 
  bind_rows(mutate(bootElnetImp, Model = "Elastic-Net")) %>% 
  mutate(Method = ifelse(Model == "Stepwise", "Least Squares", "Regularized Least Squares"),
         Model = factor(Model, levels = c("Stepwise", "Lasso", "Elastic-Net"))) %>% 
  select(-cutoff, -AUC, -ACC, -SNS, -SPC) %>% 
  gather(-Model,-Method, key = "Coefficient", value = "Estimate")

# Creating confidence intervals
plot_boot_ci <- plot_boot %>% 
  filter(Coefficient != "(Intercept)") %>% 
  group_by(Coefficient, Model, Method) %>% 
  summarise(LInf = quantile(Estimate, 0.05, na.rm = T),
            LSup = quantile(Estimate, 0.95, na.rm = T),
            Estimate = mean(Estimate, na.rm = T))
plot_boot$Estimate <- plot_boot$Estimate %>% replace_na(0)

# Creating the plot
plot_boot %>% 
  filter(Coefficient != "(Intercept)") %>% 
  ggplot(aes(x = Coefficient, y = Estimate)) + 
  geom_violin(aes(color = Model, fill = Model), color = FALSE, alpha = 0.5) + 
  geom_pointrange(aes(x = Coefficient, y = Estimate, ymin = LInf, ymax = LSup, group = Model), 
                  data = plot_boot_ci, position = position_dodge(width = 0.91), size = 0.30) +
  facet_wrap(~Model, scales = "free") + guides(fill = FALSE) +
  scale_fill_viridis_d() + scale_color_viridis_d() + 
  geom_hline(yintercept = 0) + labs(x = NULL, y = NULL) +
  scale_x_discrete(labels = c("Chol", "HDL", "LDL", "Potas", "Sod", "Trig")) +
  theme(panel.grid.major.x = element_line(color = "black")) + 
  ggpubr::theme_pubr() + theme(text = element_text(size = 16))

########### COEFFICIENTS TABLE

# Joining all coeficients together
full_table <- logist_coefs %>% 
  full_join(lasso_coefs) %>% 
  full_join(elnet_coefs)

# Creating the table
coef_table <- full_table %>% 
  transmute(Names = Names,
            StepIMP = StepIMP,
            LassoIMP = paste0(round(LassoIMP, 4), " (", 
                              round((abs(LassoIMP/FullIMPClean))*100, 2), ")"),
            ElnetIMP = paste0(round(ElnetIMP, 4), " (", 
                              round((abs(ElnetIMP/FullIMPClean))*100, 2), ")")) %>% 
  map_dfc(~ifelse(.x == "0 (0)", " - ", .x)) %>% 
  filter_at(vars(-Names), any_vars(. != " - "))

# Bootstrap coefficients table
plot_boot_ci %>% 
  mutate(CI = paste("[", round(LInf, 4) ,";", round(LSup, 4), "]"),
         Temp = paste(Estimate, CI, sep = "&")) %>% 
  select(Coefficient, Model, Temp) %>% 
  spread(key = Model, value = Temp) %>% 
  separate("Elastic-Net", sep = "&", into = c("ElnetEst", "ElnetCI")) %>% 
  separate("Lasso", sep = "&", into = c("LassoEst", "LassoCI")) %>% 
  separate("Stepwise", sep = "&", into = c("StepwiseEst", "StepwiseCI"))


########### RELATIVE VARIABLE EFFECTS PLOT

# Mutating the data
full_table %>% 
  select(Names, StepIMPClean, LassoIMP, ElnetIMP) %>% 
  mutate(StepIMPClean = ifelse(full_table$StepIMPClean == " - ", "0", StepIMPClean),
         StepIMPClean = as.numeric(StepIMPClean)) %>% 
  filter(StepIMPClean != 0 | LassoIMP != 0 | ElnetIMP !=0, Names != "(Intercept)") %>%
  modify_if(is.numeric, ~./sum(abs(.))) %>% 
  gather(key = "Modelo", val = "Coefficients", -Names) %>% 
  mutate(Modelo = recode(Modelo, 
                         ElnetIMP = "Elastic-Net", 
                         LassoIMP = "Lasso", 
                         StepIMPClean = "Stepwise"),
         Modelo = factor(Modelo, c("Stepwise", "Lasso", "Elastic-Net"), 
                         ordered = TRUE),
         Names = recode(Names, Potassium = "Potas", Sodium = "Sod")) %>% 
  # Creating the plot
  ggplot(aes(x = Names, y = Coefficients, fill = Modelo)) + 
  geom_col(width = 0.8, position = "dodge") + 
  geom_abline(intercept = 0, slope = 0) + 
  theme(panel.grid.minor.y = element_line(color = "snow")) + 
  geom_vline(xintercept = 0.5 + 1:7, color = "grey") + 
  labs(fill = NULL, x = NULL, y = NULL) + 
  scale_fill_viridis_d() + coord_flip() + ggpubr::theme_pubr() +
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  theme(text = element_text(size = 16))


########### RELATIVE VARIABLE EFFECTS PLOT

# Creating object with all ROC info
all_rocs <- tibble(Sensitivity = c(stepRocImp$sensitivities, 
                                   lassoRocImp$sensitivities, 
                                   elnetRocImp$sensitivities),
                   Specificity = c(stepRocImp$specificities, 
                                   lassoRocImp$specificities, 
                                   elnetRocImp$specificities),
                   Model = c(rep("Stepwise", length(stepRocImp$sensitivities)), 
                             rep("Lasso", length(lassoRocImp$sensitivities)),
                             rep("Elastic-Net", length(elnetRocImp$sensitivities)))) %>% 
  mutate(Model = factor(Model, levels = c("Stepwise", "Lasso", "Elastic-Net")))

# Creating ROC curve plot
g1 <- ggplot(all_rocs, aes(x = 1-Specificity, y = Sensitivity, color = Model)) + 
  geom_line(size = 1) + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), linetype = 2, size = 1, color = "black") +
  scale_color_viridis_d() + ggpubr::theme_pubr() + theme(text = element_text(size = 16))

# mutating the +632 bootstrap results
plot_boot <- mutate(bootStepImp, Model = "Stepwise") %>% 
  bind_rows(mutate(bootLassoImp, Model = "Lasso")) %>% 
  bind_rows(mutate(bootElnetImp, Model = "Elastic-Net")) %>% 
  select(Model, ACC, SNS, SPC) %>% 
  gather(-Model, key = "Measure", value = "Estimate") %>% 
  mutate(Model = factor(Model, levels = c("Stepwise", "Lasso", "Elastic-Net")),
         Measure = factor(Measure, levels = c("ACC", "SNS", "SPC")))

# Creating the confidence intervals for +632 bootstrap results
plot_boot_ci <- plot_boot %>% group_by(Measure, Model) %>% 
  summarise(LInf = quantile(Estimate, 0.05, na.rm = T),
            LSup = quantile(Estimate, 0.95, na.rm = T),
            Estimate = mean(Estimate, na.rm = T))

# Creating the bootstrap performance plot
g2 <- ggplot(plot_boot, aes(x = Measure, y = Estimate)) + 
  geom_violin(aes(fill = Model), color = FALSE, alpha = 0.5) + 
  geom_pointrange(aes(x = Measure, y = Estimate, ymin = LInf, ymax = LSup, group = Model), 
                  data = plot_boot_ci, position = position_dodge(width = 0.9), size = 0.30) +
  scale_fill_viridis_d() + scale_color_viridis_d() + 
  geom_hline(yintercept = 0.50) + labs(x = " ", y = NULL) +
  ggpubr::theme_pubr() + theme(text = element_text(size = 16))

# Combining both plots
cowplot::plot_grid(g1, g2, labels = c("(a)", "(b)"))
