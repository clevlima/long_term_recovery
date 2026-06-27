#Project "Molecular mechanisms and energetic costs of recovery from freezing in the polyextremophile midge Belgica antarctica"

#Compiled and annotated codes for part 2: energy stores

#expected files:
#carb.csv
#carb_NF.csv
#lipid.csv
#lipid_NF.csv
#protein.csv
#protein_NF.csv

#if you're getting this code from the paper online page, make sure you separate each dataset into different .csv files with the names listed above.

#set directory
#setwd("_")

#load packages
library(car) #anova
library(ggplot2) #graph
library(ggpubr) #panel

#glm carb with food ####
#load data
carb <- read.csv("carb.csv")
carb$recovery_time <- as.factor(carb$recovery_time)
carb$total_carb_per_larva <- carb$total.carb..ug/carb$N_larvae #total carb per larva, currently being used for the data analysis
carb$estimated_dry_weight_larva <- carb$estimated_dry_weight/carb$N_larvae
str(carb)
hist(carb$total_carb_per_larva)

#summary statistics
carb_summary <- carb %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(total_carb_per_larva),
    sd   = sd(total_carb_per_larva),
    se   = plotrix::std.error(total_carb_per_larva),
    .groups = "drop"
  )
carb_summary



# model carb
#mass_corrected_carb_model <- lm(mass_corrected ~ recovery_time * treatment, data = carb)
carb_model <- lm(total_carb_per_larva ~ recovery_time * treatment + estimated_dry_weight_larva, data = carb)

#verifying assumptions
plot(carb_model)

# model summary
summary(carb_model)

#anova on carb model
carb_anova <- Anova(carb_model, type = 3)
carb_anova


#graph
# it’s fine to use content per larva, but then the text needs to reflect that. You still have “micrograms per mg” in the text. If you are using the content per larva values, it should say “14.1±1.02 ug per larva,” since you don’t do any dividing by weight. Also, if you do end up doing any weight corrections, if you do mg-1, you don’t also need the slash; mg-1 means divided by mg. 


plot_carb_food <- ggplot(carb, aes(x = recovery_time, 
                                   y = total_carb_per_larva, 
                                   color = treatment, 
                                   shape = treatment,
                                   linetype = treatment,
                                   group = treatment)) +
  #geom_point(size = .5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.15, alpha = 1) +
  scale_color_manual(values = c(Control = "grey40", Frozen = "#54B9E3")) +
  scale_linetype_manual(values = c(Control = "dashed", Frozen = "solid")) +
  scale_shape_manual(values = c(1, 16)) +
  #ylim(100, 400) +
  labs(x = "Recovery time (days)", y = "Carbohydrate content (µg per larva)") +
  theme_classic() +
  theme_classic(base_size = 12) +
  theme(legend.position = c(0.15, 0.80),
        legend.title = element_blank())+
  theme(text = element_text(family = "Arial"))
plot_carb_food




#glm carb No Food ####
#load data
carb_NF <- read.csv("carb_NF.csv")
carb_NF$total_carb_per_larva <- carb_NF$total_carb..ug./carb_NF$N_larvae
carb_NF$estimated_dry_weight_larva <- carb_NF$estimated_dry_weight/carb_NF$N_larvae
str(carb_NF)
#plot(carb_NF$total_carb_per_larva~carb_NF$recovery_time)
hist(carb_NF$total_carb_per_larva)


#summary statistics
carb_NF_summary <- carb_NF %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(total_carb_per_larva),
    sd   = sd(total_carb_per_larva),
    se   = plotrix::std.error(total_carb_per_larva),
    .groups = "drop"
  )
carb_NF_summary


#model
#mass_corrected_carb_NF_model <- lm(mass_corrected ~ recovery_time * treatment, data = carb_NF)
carb_NF_model <- lm(total_carb_per_larva ~ recovery_time * treatment+estimated_dry_weight_larva, data = carb_NF, na.action = na.exclude)

#verifying assumptions
plot(carb_NF_model)

#data seem uniformly distributed but the model didn't perform too bad, so i kept the regular linear model

# model summary
summary(carb_NF_model)

#anova
carb_NF_anova <- Anova(carb_NF_model, type = 3)
carb_NF_anova


#graph
carb_NF$recovery_time <- as.factor(carb_NF$recovery_time)
plot_carb_NF <- ggplot(carb_NF, aes(x = recovery_time, 
                                    y = total_carb_per_larva, 
                                    color = treatment, 
                                    shape = treatment,
                                    linetype = treatment,
                                    group = treatment)) +
  #geom_point(size = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.15, alpha = 1) +
  scale_color_manual(values = c(Control = "grey40", Frozen = "#54B9E3")) +
  scale_linetype_manual(values = c(Control = "dashed", Frozen = "solid")) +
  scale_shape_manual(values = c(1, 16)) +
  #ylim(0, 400) +
  labs(x = "Recovery time (days)", y = "Carbohydrate content (µg per larva)") +
  theme_classic(base_size = 12) +
  theme(legend.position = c(0.85, 0.90),
        legend.title = element_blank())+
  theme(text = element_text(family = "Arial"))

plot_carb_NF



#model lipid with food ####
#load data
lipid <- read.csv("lipid.csv")
lipid$recovery_time <- as.factor(lipid$recovery_time)
lipid$total_lipid_per_larva <- lipid$total_lipid/lipid$N_larvae
lipid$estimated_dry_weight_larva <- lipid$estimated_dry_weight/lipid$N_larvae


str(lipid)
#plot(lipid$total_lipid_per_larva~lipid$recovery_time)
hist(lipid$total_lipid_per_larva)

#summary statistics
lipid_summary <- lipid %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(total_lipid_per_larva),
    sd   = sd(total_lipid_per_larva),
    se   = plotrix::std.error(total_lipid_per_larva),
    .groups = "drop"
  )
lipid_summary


#model
#lipid_model <- lm(total_lipid_per_larva ~ recovery_time * treatment,
#                   data = lipid, na.action = na.exclude)

#verifying assumptions
#plot(lipid_model)

#model violated regular lm assumptions, when i run hist(), distribution seem skewed to the left, so i modeled a glm with gamma distribution instead!!

lipid_model_gamma <- glm(total_lipid_per_larva ~ recovery_time * treatment + estimated_dry_weight_larva, family = Gamma(link = "log"), data = lipid, na.action = na.exclude)

#verifying assumptions
plot(lipid_model_gamma)

#much better, so i will go ahead and check the results

# model summary
summary(lipid_model_gamma)

#anova
lipid_anova <- Anova(lipid_model_gamma, type = 3)
lipid_anova

#graph
plot_lipid_food <- ggplot(lipid, aes(x = recovery_time, 
                                     y = total_lipid_per_larva, 
                                     color = treatment, 
                                     shape = treatment,
                                     linetype = treatment,
                                     group = treatment)) +
  #geom_point(size = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.15, alpha = 1) +
  scale_color_manual(values = c(Control = "grey40", Frozen = "#54B9E3")) +
  scale_linetype_manual(values = c(Control = "dashed", Frozen = "solid")) +
  scale_shape_manual(values = c(1, 16)) +
 # ylim(10, 60) +
  labs(x = "Recovery time (days)", y = "Lipid content (µg per larva)") +
  theme_classic() +
  theme_classic(base_size = 12) +
  #theme(legend.position = c(0.85, 0.85),
  #      legend.title = element_blank())+
  theme(legend.position = "none") +
  theme(text = element_text(family = "Arial"))
plot_lipid_food






#model lipid NF ####
#load data
lipid_NF <- read.csv("lipid_NF.csv")
lipid_NF$total_lipid_per_larva <- lipid_NF$total_lipid/lipid_NF$N_larvae
lipid_NF$estimated_dry_weight_larva <- lipid_NF$estimated_dry_weight/lipid_NF$N_larvae

str(lipid_NF)
#plot(lipid_NF$total_lipid_per_larva~lipid_NF$recovery_time)
hist(lipid_NF$total_lipid_per_larva)

#summary statistics
lipid_NF_summary <- lipid_NF %>%
  group_by(recovery_time) %>%
  summarise(
    mean = mean(total_lipid_per_larva),
    sd   = sd(total_lipid_per_larva),
    se   = plotrix::std.error(total_lipid_per_larva),
    .groups = "drop"
  )
lipid_NF_summary


#model
#mass_corrected_lipid_NF_model <- lm(mass_corrected ~ recovery_time * treatment, data = lipid_NF)
lipid_NF_model <- lm(total_lipid_per_larva ~ recovery_time*treatment + estimated_dry_weight_larva, data = lipid_NF, na.action = na.exclude)
#verifying assumptions
plot(lipid_NF_model)

#model summary
summary(lipid_NF_model)

#anova
lipid_NF_anova <- Anova(lipid_NF_model, type = 3)
lipid_NF_anova


#note that the effect of recovery_time on lipid content did not change (still significant) when using the total lipid per larva as response variable

#graph
lipid_NF$recovery_time <- as.factor(lipid_NF$recovery_time)
plot_lipid_NF <- ggplot(lipid_NF, aes(x = recovery_time, 
                                      y = total_lipid_per_larva, 
                                      color = treatment, 
                                      shape = treatment,
                                      linetype = treatment,
                                      group = treatment)) +
  #  geom_point(size = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.15, alpha = 1) +
  scale_color_manual(values = c(Control = "grey40", Frozen = "#54B9E3")) +
  scale_linetype_manual(values = c(Control = "dashed", Frozen = "solid")) +
  scale_shape_manual(values = c(1, 16)) +
  #  ylim(20, 100) +
  labs(x = "Recovery time (days)", y = "Lipid content (µg per larva)") +
  theme_classic(base_size = 12) +
  #theme(legend.position = c(0.15, 0.10),
  #      legend.title = element_blank())+
  theme(legend.position = "none") +
  theme(text = element_text(family = "Arial"))
plot_lipid_NF






#model protein with food ####
#load data
protein <- read.csv("protein.csv")
protein$recovery_time <- as.factor(protein$recovery_time)
protein$total_protein_per_larva <- protein$total_protein/protein$N_larvae
protein$estimated_dry_weight_larva <- protein$estimated_dry_weight/protein$N_larvae
str(protein)
hist(protein$total_protein_per_larva)


#summary statistics
protein_summary <- protein %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(total_protein_per_larva),
    sd   = sd(total_protein_per_larva),
    se   = plotrix::std.error(total_protein_per_larva),
    .groups = "drop"
  )
protein_summary


#model
#mass_corrected_protein_model <- lm(mass_corrected ~ recovery_time * treatment, data = protein)
protein_model <- lm(total_protein_per_larva ~ recovery_time * treatment + estimated_dry_weight_larva, data = protein)

#verifying assumptions
plot(protein_model)

#model summary
summary(protein_model)

#anova
protein_anova <- Anova(protein_model, type = 3)
protein_anova

#the anova reports a marginal but not significant interaction between time and treatment. However i will run a post hoc to be sure

#tukey post hoc
#emm_model <- emmeans(protein_model, specs = pairwise ~ treatment | recovery_time)
#emm_model$contrasts

#looks like there is a significant difference between frozen and control at days 7 and 15 (at least when dry weight is not included in the model), both of which larvae that had been frozen showed a significantly higher protein content than untreated larvae
#Authors discussed and agreed that since 1) protein is higher in frozen larvae and 2) there was still no statistical difference between frozen and control for carbs and lipids, the conclusions are the same: no energy depletion following freezing stress. 
#Additionally, because the results did not change when using mass_corrected vs total macronutrient per larva as response variables, there is a possibility that the difference in weight observed in frozen larvae from which protein was extracted is likely to be an experimental artifact 
  ##we might have placed disproportionally larger larvae in the frozen groups vs the control groups for this macronutrient

#graph
plot_protein_food <- ggplot(protein, aes(x = recovery_time, 
                                         y = total_protein_per_larva, 
                                         color = treatment, 
                                         shape = treatment,
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
  scale_shape_manual(values = c(1, 16)) +
  #ylim(70, 120) +
  labs(x = "Recovery time (days)", y = "Protein content (µg per larva)") +
  theme_classic() +
  theme_classic(base_size = 12) +
  #theme(legend.position = c(0.85, 0.85),
  #      legend.title = element_blank())+
  theme(legend.position = "none") +
  theme(text = element_text(family = "Arial"))
plot_protein_food



#model protein NF ####
#load data
protein_NF <- read.csv("protein_NF.csv")
protein_NF$estimated_dry_weight_larva <- protein_NF$estimated_dry_weight/protein_NF$N_larvae
protein_NF$total_protein_per_larva <- protein_NF$total_protein/protein_NF$N_larvae

#check outlier - protein_NF: C7_4??
str(protein_NF)
hist(protein_NF$total_protein_per_larva)

#it really does seem like an outlier, so i will remove it. 
#the criteria was the interquartile range between 0.25 and 0.75, where anything out of this range is considered an outlier and removed from the original dataset
dim(protein_NF) 
Q1 <- quantile(protein_NF$total_protein_per_larva, .25)
Q3 <- quantile(protein_NF$total_protein_per_larva, .75)
IQR <- IQR(protein_NF$total_protein_per_larva)
protein_NF_cleaned<- subset(protein_NF, protein_NF$total_protein_per_larva > (Q1 - 1.5*IQR) & protein_NF$total_protein_per_larva < (Q3 + 1.5*IQR))
dim(protein_NF_cleaned) 
str(protein_NF_cleaned)
hist(protein_NF_cleaned$total_protein_per_larva)


#summary statistics
protein_NF_summary <- protein_NF_cleaned %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(total_protein_per_larva),
    sd   = sd(total_protein_per_larva),
    se   = plotrix::std.error(total_protein_per_larva),
    .groups = "drop"
  )
protein_NF_summary



#model
#mass_corrected_protein_NF_model <- lm(mass_corrected ~ recovery_time * treatment, data = protein_NF)
protein_NF_model <- lm(total_protein_per_larva ~ recovery_time * treatment + estimated_dry_weight_larva, data = protein_NF_cleaned, na.action = na.exclude)
#verifying assumptions
plot(protein_model)

#model summary
summary(protein_NF_model)

#anova
protein_NF_anova <- Anova(protein_NF_model, type = 3)
protein_NF_anova

#graph
protein_NF_cleaned$recovery_time <- as.factor(protein_NF_cleaned$recovery_time)
plot_protein_NF <- ggplot(protein_NF_cleaned, aes(x = recovery_time, 
                                          y = total_protein_per_larva, 
                                          color = treatment, 
                                          shape = treatment,
                                          linetype = treatment,
                                          group = treatment)) +
  #  geom_point(size = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.15, alpha = 1) +
  scale_color_manual(values = c(Control = "grey40", Frozen = "#54B9E3")) +
  scale_linetype_manual(values = c(Control = "dashed", Frozen = "solid")) +
  scale_shape_manual(values = c(1, 16)) +
    #ylim(27, 38) +
    labs(x = "Recovery time (days)", y = "Protein content content (µg per larva)") +
  theme_classic(base_size = 12) +
  # theme(legend.position = c(0.15, 0.10),
  #      legend.title = element_blank()) +
  theme(legend.position = "none") +
  theme(text = element_text(family = "Arial"))
plot_protein_NF



# panel food 
#make a list with plots for panel (need to run survival script first)
plot_list_all <- list(
  plot_carb_food + theme(axis.title.x = element_blank()) + theme(legend.position = "none"),
  plot_lipid_food + theme(axis.title.x = element_blank()) + theme(legend.position = "none"),
  plot_protein_food + theme(legend.position = "none")
)

#plot panel
plot_all <- ggarrange(plotlist = plot_list_all, 
                      labels = c("A", "B", "C"), 
                      hjust = 0.15, ncol = 1, nrow = 3,
                      #common.legend = FALSE, legend = "top", 
                      font.label = list(size = 12, face = "bold", family = "Arial")) + theme(text = element_text(family = "Arial")) 
plot_all



# panel NF
plot_list_all <- list(
  plot_carb_NF + theme(axis.title.x = element_blank()) + theme(legend.position = "none"),
  plot_lipid_NF + theme(axis.title.x = element_blank()) + theme(legend.position = "none"),
  plot_protein_NF + theme(legend.position = "none")
)

#plot panel
plot_all <- ggarrange(plotlist = plot_list_all, 
                      labels = c("A", "B", "C"), 
                      hjust = 0.15, ncol = 1, nrow = 3,
                      #  common.legend = TRUE, legend = "top", 
                      font.label = list(size = 13, face = "bold", family = "Arial")) + theme(text = element_text(family = "Arial")) 
plot_all

#
sessionInfo()
#