source("https://raw.githubusercontent.com/Khatri-Lab/SUBSPACE/main/Code/subspace_load_data.R")

savemore_baseline <- savemore%>%
  filter(timepoint == "baseline")%>%
  mutate(high_lymphoid_score = ifelse(lymphoid_z_score >= median(filter(scores_w_pheno, severity == "non-severe",infection_type == "Viral", timepoint == "baseline", site != "cchmc", condition == "infected")$lymphoid_z_score, na.rm = T), "High", "Low"),
         high_myeloid_score = ifelse(myeloid_z_score >= median(filter(scores_w_pheno, severity == "non-severe", infection_type == "Viral", timepoint == "baseline", site != "cchmc", condition == "infected")$myeloid_z_score, na.rm = T), "High","Low"),
         subgroup = ifelse(high_lymphoid_score == "High" & high_myeloid_score == "High", "system-wide",
                           ifelse(high_lymphoid_score == "High" & high_myeloid_score == "Low", "lymphoid",
                                  ifelse(high_myeloid_score == "High" & high_lymphoid_score == "Low", "myeloid",
                                         ifelse(high_myeloid_score == "Low" & high_lymphoid_score == "Low", "balanced",NA)))))

savemore_baseline$subgroup = factor(savemore_baseline$subgroup, levels = c("balanced","lymphoid", "myeloid","system-wide"))
savemore_baseline$treatment= factor(savemore_baseline$treatment, levels = c("placebo","anakinra"))
savemore_baseline$high_lymphoid_score = factor(savemore_baseline$high_lymphoid_score, levels = c("Low","High"))
savemore_baseline$high_myeloid_score = factor(savemore_baseline$high_myeloid_score, levels = c("Low","High"))

fisher_low_savemore = fisher.test(filter(savemore_baseline, high_lymphoid_score == "Low")$treatment, filter(savemore_baseline, high_lymphoid_score == "Low")$d30_mort)
fisher_high_savemore = fisher.test(filter(savemore_baseline, high_lymphoid_score == "High")$treatment, filter(savemore_baseline, high_lymphoid_score == "High")$d30_mort)

prop_savemore <- savemore_baseline%>%
  group_by(high_lymphoid_score, treatment)%>%
  dplyr::summarize(mortality_prop = mean(d30_mort == 1, na.rm = TRUE)*100)

table(savemore_baseline$treatment, savemore_baseline$subgroup)

table(savemore_baseline$treatment, savemore_baseline$high_lymphoid_score, savemore_baseline$d30_mort)

savemore_barplot <- ggplot(prop_savemore, aes(x = high_lymphoid_score, y = mortality_prop, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  #adjust labels of x axis to capitalize words
  scale_x_discrete(labels = function(x) stringr::str_to_title(x)) +
  scale_fill_manual(values = c("placebo" = "#999999", "anakinra" = "goldenrod3"))+
  theme(legend.position = "none")+
  xlab("Lymphoid Dysregulation Score") +
  ylab("28-day Mortality (%)")+
  ylim(0,25)+
  annotate("text",x=1, y = 9, label = paste("p =",trail.sigfig(fisher_low_savemore$p.value,4)))+
  annotate("text",x=2, y = 23, label = paste("p =",trail.sigfig(fisher_high_savemore$p.value,4)))+
  theme(panel.grid.major.x = element_blank())

print(savemore_barplot)

km_savemore <- survfit2(Surv(survdays,d30_mort==1) ~ treatment, data = filter(savemore_baseline,high_lymphoid_score == "High"), conf.type = "log-log")

cox_savemore <-coxph(Surv(survdays,d30_mort==1) ~ treatment + age + sex, data = filter(savemore_baseline,high_lymphoid_score == "High"))
summary(cox_savemore)

savemore_survplot<- ggsurvfit(km_savemore, size = 1.5) +
  ylim(.4,1) +
  labs(
    x = "Days",
    y = "Survival"
  ) +
  theme_classic() +
  theme(legend.position =  "none") +
  scale_color_manual(name = "Treatment",values = c("#999999","goldenrod3")) +
  annotate("text",x = 4.4, y = 0.7, hjust = 0, label =  "HR: 0.06 (0.008-0.53)")+
  annotate("text",x = 4.4, y = 0.65, hjust = 0, label =  "p = 0.01")+
  theme(panel.grid.major.y = element_line(color = "grey", size = 0.5))

savemore_survplot

prop_savemore_subgroup <- savemore_baseline%>%
  group_by(subgroup, treatment)%>%
  summarize(mortality_prop = mean(d30_mort == 1, na.rm = TRUE)*100)

savemore_subgroup_barplot <- ggplot(prop_savemore_subgroup, aes(x = subgroup, y = mortality_prop, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  #adjust labels of x axis to capitalize words
  scale_x_discrete(labels = function(x) stringr::str_to_title(x)) +
  scale_fill_manual(values = c("placebo" = "#999999", "anakinra" = "goldenrod3"))+
  theme(legend.position = "none")+
  xlab("Subgroup") +
  ylab("28-day Mortality (%)")+
  theme(panel.grid.major.x = element_blank())

print(savemore_subgroup_barplot)

##########VICTAS##########
victas_treatment_eval <- victas %>%
  mutate(high_lymphoid_score = ifelse(lymphoid_z_score >= median(filter(scores_w_pheno, severity == "severe" & site != "cchmc",  timepoint == "baseline", condition == "infected")$lymphoid_z_score, na.rm = T), "High", "Low"),
         high_myeloid_score = ifelse(myeloid_z_score >= median(filter(scores_w_pheno, severity == "severe" & site != "cchmc", timepoint == "baseline", condition == "infected")$myeloid_z_score, na.rm = T), "High","Low"),
         subgroup = ifelse(high_lymphoid_score == "High" & high_myeloid_score == "High", "system-wide",
                           ifelse(high_lymphoid_score == "High" & high_myeloid_score == "Low", "lymphoid",
                                  ifelse(high_myeloid_score == "High" & high_lymphoid_score == "Low", "myeloid",
                                         ifelse(high_myeloid_score == "Low" & high_lymphoid_score == "Low", "balanced",NA)))))

#Filters out patients who received open-label steroids
victas_no_open <- victas_treatment_eval %>%
  filter(!grepl("Open", treatment), !is.na(treatment))

victas_no_open$subgroup = factor(victas_no_open$subgroup, levels = c("balanced","lymphoid", "myeloid","system-wide"))

table(victas_no_open$treatment, victas_no_open$subgroup)

table(victas_no_open$treatment, victas_no_open$high_lymphoid_score, victas_no_open$d30_mort)

victas_no_open$high_lymphoid_score = factor(victas_no_open$high_lymphoid_score, levels = c("Low","High"))

victas_no_open$treatment = factor(victas_no_open$treatment, levels = c("Placebo","HAT"))

victas_prop <- victas_no_open%>%
  group_by(high_lymphoid_score, treatment)%>%
  dplyr::summarize(mort_prop = mean(d30_mort == 1, na.rm = TRUE)*100)


fisher_low_victas = fisher.test(filter(victas_no_open, high_lymphoid_score == "Low")$treatment, filter(victas_no_open, high_lymphoid_score == "Low")$d30_mort)
fisher_high_victas = fisher.test(filter(victas_no_open, high_lymphoid_score == "High")$treatment, filter(victas_no_open, high_lymphoid_score == "High")$d30_mort)

victas_barplot <- ggplot(victas_prop, aes(x = high_lymphoid_score, y = mort_prop, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge", reverse = T) +
  theme_minimal() +
  scale_fill_manual(values = c("Placebo" = "#999999", "HAT" = "#DC3023"))+
  #adjust labels of x axis to capitalize words
  scale_x_discrete(labels = function(x) stringr::str_to_title(x)) +
  theme(legend.position = "none")+
  xlab("Lymphoid Dysregulation Score") +
  ylab("30-day Mortality (%)")+
  #  ylim(0,50)+
  annotate("text",x=1, y = 20, label = paste("p =",format.pval(fisher_low_victas$p.value, 3)))+
  annotate("text",x=2, y = 48, label = paste("p =",format.pval(fisher_high_victas$p.value,2)))+
  #remove vertical gridlines
  theme(panel.grid.major.x = element_blank())

print(victas_barplot)

log_model_victas <- glm(d30_mort ~ treatment + age + sex ,  data = filter(victas_no_open,high_lymphoid_score == "High"), family = "binomial")
summary(log_model_victas)

km_victas <- survfit2(Surv(survdays,d30_mort==1) ~ treatment, data = filter(victas_no_open,high_lymphoid_score == "High"), conf.type = "log-log")

cox_victas <-coxph(Surv(survdays,d30_mort==1) ~ treatment + sex +age  , filter(victas_no_open,high_lymphoid_score == "High"))
summary(cox_victas)

victas_survplot<- ggsurvfit(km_victas, size = 1.5) +
  ylim(.4,1) +
  labs(
    x = "Days",
    y = "Survival"
  ) +
  theme_classic() +
  theme(legend.position =  "none") +
  scale_color_manual(name = "Treatment",values = c("#999999","#DC3023"))+
  annotate("text",x = 4.4, y = 0.6, hjust = 0, label =  "HR: 0.22 (0.06-0.85)\np = 0.03")+
  theme(panel.grid.major.y = element_line(color = "grey", size = 0.5))

victas_survplot

victas_prop_subgroup <- victas_no_open%>%
  group_by(subgroup, treatment)%>%
  summarize(mort_prop = mean(d30_mort == 1, na.rm = TRUE)*100)

victas_subgroup_barplot <- ggplot(victas_prop_subgroup, aes(x = subgroup, y = mort_prop, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  #adjust labels of x axis to capitalize words
  scale_x_discrete(labels = function(x) stringr::str_to_title(x)) +
  scale_fill_manual(values = c("Placebo" = "#999999", "HAT" = "#DC3023"))+
  theme(legend.position = "none")+
  xlab("Subgroup") +
  ylab("30-day Mortality (%)")+
  theme(panel.grid.major.x = element_blank())

print(victas_subgroup_barplot)

#########VANISH##########
vanish_prop <- vanish%>%
  filter(Characteristics.drug2.per.protocol. != "none")%>%
  group_by(high_lymphoid, treatment)%>%
  dplyr::summarize(mortality_prop = mean(d28_death_yn == 1, na.rm = TRUE)*100)

vanish_prop$treatment = factor(vanish_prop$treatment, levels = c("no hydrocortisone","hydrocortisone"))

table(filter(vanish,Characteristics.drug2.per.protocol. != "none")$treatment, filter(vanish,Characteristics.drug2.per.protocol. != "none")$high_lymphoid, filter(vanish,Characteristics.drug2.per.protocol. != "none")$d28_death_yn)

fisher_low_vanish = fisher.test(filter(vanish, high_lymphoid == "low")$treatment, filter(vanish, high_lymphoid == "low")$d28_death_yn)
fisher_high_vanish = fisher.test(filter(vanish, high_lymphoid == "high")$treatment, filter(vanish, high_lymphoid == "high")$d28_death_yn)

vanish_barplot <- ggplot(vanish_prop, aes(x = high_lymphoid, y = mortality_prop, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  scale_fill_manual(values = c("no hydrocortisone" = "#999999", "hydrocortisone" = "#81BFE9"))+
  #adjust labels of x axis to capitalize words
  scale_x_discrete(labels = function(x) stringr::str_to_title(x)) +
  theme(legend.position = "none")+
  xlab("Lymphoid Dysregulation Score") +
  ylab("28-day Mortality (%)")+
  ylim(0,50)+
  annotate("text",x=1, y = 45, label = paste("p =",trail.sigfig(fisher_low_vanish$p.value,2)))+
  annotate("text",x=2, y = 38, label = paste("p =",trail.sigfig(fisher_high_vanish$p.value,2)))+
  theme(panel.grid.major.x = element_blank())

print(vanish_barplot)

vanish_prop_subgroup<- vanish%>%
  filter(Characteristics.drug2.per.protocol. != "none")%>%
  group_by(subgroup, treatment)%>%
  summarize(mortality_prop = mean(d28_death_yn == 1, na.rm = TRUE)*100)

vanish_subgroup_barplot <- ggplot(vanish_prop_subgroup, aes(x = subgroup, y = mortality_prop, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  #adjust labels of x axis to capitalize words
  scale_x_discrete(labels = function(x) stringr::str_to_title(x)) +
  scale_fill_manual(values = c("no hydrocortisone" = "#999999", "hydrocortisone" = "#81BFE9"))+
  theme(legend.position = "none")+
  xlab("Subgroup") +
  ylab("28-day Mortality (%)")+
  theme(panel.grid.major.x = element_blank())

print(vanish_subgroup_barplot)

