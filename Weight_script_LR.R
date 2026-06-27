#weight throughout recovery - data analysis
library(lme4)
library(car)
library(ggplot2)


#food ####
weight <- read.csv("weight.csv")
weight$estimated_dry_weight <- weight$weight_per_larva*0.27
str(weight)
#hist(weight$weight)
#hist(weight$weight_per_larva)
hist(weight$estimated_dry_weight)
weight$recovery_time <- as.factor(weight$recovery_time)
#weight$molecule <- as.factor(weight$molecule)


#plot
plot_weight <- ggplot(weight, aes(x = recovery_time, 
                                  y = estimated_dry_weight, 
                                  color = treatment, 
                                 # shape = treatment,
                                  linetype = treatment,
                                  group = treatment)) +
  # geom_point(size = .5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.15, 
               alpha = 1) +
  scale_color_manual(values = c(Control = "grey40", Frozen = "#54B9E3")) +
  scale_linetype_manual(values = c(Control = "dashed", Frozen = "solid")) +
  #ylim(.3, 1.1) +
  labs(x = "Recovery time (days)", y = "Estimated dry weight per larva (mg)") + 
  theme_classic(base_size = 12) +
  theme(legend.position = c(.15, .15)) + theme(legend.title=element_blank())
plot_weight
#plot_protein_food


#model
model_food <- lm(estimated_dry_weight ~ recovery_time * treatment, 
                 data = weight)
summary(model_food)

model_food_anova <- Anova(model_food, type=3)
model_food_anova


#tukey post hoc
library(emmeans)
emm_model <- emmeans(model_food, specs = pairwise ~ treatment | recovery_time)
emm_model$contrasts
 

###### individual macronutrients #####
#subset the data to test individual molecules
protein_data <- subset(weight, molecule == "protein")
str(protein_data)
#hist(protein_data$weight_per_larva)
hist(protein_data$estimated_dry_weight)


#plot protein averages
plot_protein <- ggplot(protein_data, aes(x = recovery_time, 
                                         y = estimated_dry_weight, 
                                         color = treatment, 
                                         # shape = treatment,
                                         linetype = treatment,
                                         group = treatment)) +
  # geom_point(size = .5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.15, 
               alpha = 1) +
  scale_color_manual(values = c(Control = "grey40", Frozen = "#54B9E3")) +
  scale_linetype_manual(values = c(Control = "dashed", Frozen = "solid")) +
  #ylim(.3, 1.1) +
  labs(x = "Recovery time (days)", y = "Estimated dry weight per larva (mg)") + 
  theme_classic(base_size = 12) +
  theme(legend.position = c(.15, .15)) + theme(legend.title=element_blank())
plot_protein


#fit the model only for proteins
model_protein <- lm(estimated_dry_weight ~ recovery_time * treatment,
                     data = protein_data)
plot(model_protein)
summary(model_protein)

#anova
Anova(model_protein, type=3)


#tukey post hoc
emm_model <- emmeans(model_protein, specs = pairwise ~ treatment | recovery_time)
emm_model$contrasts




##

#lipids
lipid_data <- subset(weight, molecule == "lipid")
str(lipid_data)
#hist(lipid_data$weight_per_larva)
hist(lipid_data$estimated_dry_weight)
dim(lipid_data) 
Q1 <- quantile(lipid_data$estimated_dry_weight, .25)
Q3 <- quantile(lipid_data$estimated_dry_weight, .75)
IQR <- IQR(lipid_data$estimated_dry_weight)
lipid_NF_cleaned<- subset(lipid_data, lipid_data$estimated_dry_weight > (Q1 - 1.5*IQR) & lipid_data$estimated_dry_weight < (Q3 + 1.5*IQR))
dim(lipid_NF_cleaned) 
hist(lipid_NF_cleaned$estimated_dry_weight)


#plot protein averages
plot_lipid <- ggplot(lipid_NF_cleaned, aes(x = recovery_time, 
                                     y = estimated_dry_weight, 
                                     color = treatment, 
                                     # shape = treatment,
                                     linetype = treatment,
                                     group = treatment)) +
  # geom_point(size = .5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.15, 
               alpha = 1) +
  scale_color_manual(values = c(Control = "grey40", Frozen = "#54B9E3")) +
  scale_linetype_manual(values = c(Control = "dashed", Frozen = "solid")) +
  #ylim(.3, 1.1) +
  labs(x = "Recovery time (days)", y = "Estimated dry weight per larva (mg)") + 
  theme_classic(base_size = 12) +
  theme(legend.position = c(.15, .85)) + theme(legend.title=element_blank())
plot_lipid


#fit the model only for lipids
model_lipid <- lm(estimated_dry_weight ~ recovery_time * treatment,
                     data = lipid_NF_cleaned)
plot(model_lipid)
summary(model_lipid)


#glm with gamma distribution
#model_lipid <- glm(estimated_dry_weight ~ recovery_time * treatment,
#                   family = Gamma(link = "log"), data = lipid_data)
#plot(model_lipid)
#summary(model_lipid)

#anova
Anova(model_lipid, type=3)

#tukey post hoc
#emm_model <- emmeans(model_lipid, specs = pairwise ~ treatment | recovery_time)
#emm_model$contrasts

##



#carbs
carb_data <- subset(weight, molecule == "carb")
str(carb_data)
hist(carb_data$weight_per_larva)
hist(carb_data$estimated_dry_weight)


#plot protein averages
plot_carb <- ggplot(carb_data, aes(x = recovery_time, 
                                   y = estimated_dry_weight, 
                                   color = treatment, 
                                   # shape = treatment,
                                   linetype = treatment,
                                   group = treatment)) +
  # geom_point(size = .5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.15, 
               alpha = 1) +
  scale_color_manual(values = c(Control = "grey40", Frozen = "#54B9E3")) +
  scale_linetype_manual(values = c(Control = "dashed", Frozen = "solid")) +
  #ylim(.3, 1.1) +
  labs(x = "Recovery time (days)", y = "Estimated dry weight per larva (mg)") + 
  theme_classic(base_size = 12) +
  theme(legend.position = c(.15, .85)) + theme(legend.title=element_blank())
plot_carb



#fit the model only for carbs
model_carb <- lm(estimated_dry_weight ~ recovery_time * treatment,
                     data = carb_data)
plot(model_carb)
summary(model_carb)

#anova
Anova(model_carb, type=3)

#tukey post hoc
#emm_model <- emmeans(model_carb, specs = pairwise ~ treatment | recovery_time)
#emm_model$contrasts




#
#no food ####
weight_NF <- read.csv("weight_NF.csv")
weight_NF$estimated_dry_weight <- weight_NF$weight_per_larva*0.27
str(weight_NF)
weight_NF$recovery_time <- as.factor(weight_NF$recovery_time)
#hist(weight_NF$weight_per_larva)


#plot
plot_weight_NF <- ggplot(weight_NF, aes(x = recovery_time, 
                                        y = estimated_dry_weight, 
                                        color = treatment, 
                                        # shape = treatment,
                                        linetype = treatment,
                                        group = treatment)) +
  # geom_point(size = .5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.15, 
               alpha = 1) +
  scale_color_manual(values = c(Control = "grey40", Frozen = "#54B9E3")) +
  scale_linetype_manual(values = c(Control = "dashed", Frozen = "solid")) +
  #ylim(.3, 1.1) +
#  ylim(.75, 1.5) +
  labs(x = "Recovery time (days)", y = "Estimated dry weight per larva (mg)") + 
  theme_classic(base_size=12) +
  theme(legend.position = c(.15, .95)) + theme(legend.title=element_blank())
plot_weight_NF


#model
model_NF <- lm(estimated_dry_weight ~ recovery_time * treatment,
               data = weight_NF)
plot(model_NF)
summary(model_NF)

#anova
Anova(model_NF, type=3)

#tukey post hoc
emm_model <- emmeans(model_NF, specs = pairwise ~ treatment | recovery_time)
emm_model$contrasts




#proteins
NF_protein_data <- subset(weight_NF, molecule == "protein")
#str(NF_protein_data)
hist(NF_protein_data$weight_per_larva)


#plot protein averages
plot_NF_protein <- ggplot(NF_protein_data, aes(x = recovery_time, 
                                               y = estimated_dry_weight, 
                                               color = treatment, 
                                               # shape = treatment,
                                               linetype = treatment,
                                               group = treatment)) +
  # geom_point(size = .5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.15, 
               alpha = 1) +
  scale_color_manual(values = c(Control = "grey40", Frozen = "#54B9E3")) +
  scale_linetype_manual(values = c(Control = "dashed", Frozen = "solid")) +
  #ylim(.3, 1.1) +
  labs(x = "Recovery time (days)", y = "Estimated dry weight per larva (mg)") + 
  theme_classic(base_size = 12) +
  theme(legend.position = c(.15, .15)) + theme(legend.title=element_blank())
plot_NF_protein


#fit the model only for proteins
model_NF_protein <- lm(estimated_dry_weight ~ recovery_time * treatment,
                    data = NF_protein_data)
plot(model_NF_protein)
summary(model_NF_protein)

#anova
Anova(model_NF_protein, type=3)


#tukey post hoc
#emm_model <- emmeans(model_NF_protein, specs = pairwise ~ treatment | recovery_time)
#emm_model$contrasts



#lipids
NF_lipid_data <- subset(weight_NF, molecule == "lipid")
#str(NF_lipid_data)
hist(NF_lipid_data$weight_per_larva)


#plot protein averages
plot_NF_lipid <- ggplot(NF_lipid_data, aes(x = recovery_time, 
                                           y = estimated_dry_weight,
                                           color = treatment, 
                                           # shape = treatment,
                                           linetype = treatment,
                                           group = treatment)) +
  # geom_point(size = .5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.15, 
               alpha = 1) +
  scale_color_manual(values = c(Control = "grey40", Frozen = "#54B9E3")) +
  scale_linetype_manual(values = c(Control = "dashed", Frozen = "solid")) +
  #ylim(.3, 1.1) +
  labs(x = "Recovery time (days)", y = "Estimated dry weight per larva (mg)") + 
  theme_classic(base_size = 12) +
  theme(legend.position = c(.15, .15)) + theme(legend.title=element_blank())
plot_NF_lipid


#fit the model only for lipids
model_NF_lipid <- lm(estimated_dry_weight ~ recovery_time * treatment,
                   data = NF_lipid_data)
plot(model_NF_lipid)
summary(model_NF_lipid)

#anova
Anova(model_NF_lipid, type=3)

#tukey post hoc
#emm_model <- emmeans(model_NF_lipid, specs = pairwise ~ treatment | recovery_time)
#emm_model$contrasts


##


#carbs
NF_carb_data <- subset(weight_NF, molecule == "carb")
str(NF_carb_data)
hist(NF_carb_data$weight_per_larva)
##not normally distributed


#plot protein averages
plot_NF_carb <- ggplot(NF_carb_data, aes(x = recovery_time, 
                                         y = estimated_dry_weight,
                                         color = treatment, 
                                         # shape = treatment,
                                         linetype = treatment,
                                         group = treatment)) +
  # geom_point(size = .5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.15, 
               alpha = 1) +
  scale_color_manual(values = c(Control = "grey40", Frozen = "#54B9E3")) +
  scale_linetype_manual(values = c(Control = "dashed", Frozen = "solid")) +
  #ylim(.3, 1.1) +
  labs(x = "Recovery time (days)", y = "Estimated dry weight per larva (mg)") + 
  theme_classic(base_size = 12) +
  theme(legend.position = c(.15, .85)) + theme(legend.title=element_blank())
plot_NF_carb


#fit the model only for carbs
model_NF_carb <- lm(estimated_dry_weight ~ recovery_time * treatment,
                  data = NF_carb_data, na.action = na.omit)
plot(model_NF_carb)
summary(model_NF_carb)

#anova
Anova(model_NF_carb, type=3)

#tukey post hoc
#emm_model <- emmeans(model_NF_carb, specs = pairwise ~ treatment | recovery_time)
#emm_model$contrasts

#





# panel food 
plot_list_all <- list(
  plot_carb + theme(axis.title.x = element_blank()) + theme(legend.position = "none"),
  plot_lipid + theme(axis.title.x = element_blank()) + theme(legend.position = "none"),
  plot_protein + theme(legend.position = "none")
)

#plot panel
plot_all <- ggarrange(plotlist = plot_list_all, 
                      labels = c("A) Carbohydrate", "B) Lipid", "C) Protein"), 
                      hjust = 0.15, ncol = 1, nrow = 3,
                      #common.legend = FALSE, legend = "top", 
                      font.label = list(size = 13, face = "bold", family = "Arial")) + theme(text = element_text(family = "Arial")) 
plot_all


# panel NF
plot_list_all <- list(
  plot_NF_carb + theme(axis.title.x = element_blank()) + theme(legend.position = "none"),
  plot_NF_lipid + theme(axis.title.x = element_blank()) + theme(legend.position = "none"),
  plot_NF_protein + theme(legend.position = "none")
)

#plot panel
plot_all <- ggarrange(plotlist = plot_list_all, 
                      labels = c("A) Carbohydrate", "B) Lipid", "C) Protein"), 
                      hjust = 0.15, ncol = 1, nrow = 3,
                      #  common.legend = TRUE, legend = "top", 
                      font.label = list(size = 13, face = "bold", family = "Arial")) + theme(text = element_text(family = "Arial")) 
plot_all

#
sessionInfo()
#
