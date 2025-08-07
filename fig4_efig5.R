source("https://raw.githubusercontent.com/Khatri-Lab/SUBSPACE/main/subspace_load_data.R")



#####Public violin plots#####
#Lymphoid
rawr::jt.test(public_score_table$lymphoid_score,as.factor(public_score_table$severity_grades))

public_severity_lymphoid_score_plot <- ggplot(filter(public_score_table, !is.na(severity_grades)), aes(x = as.factor(severity_grades), y = lymphoid_z_score, fill = as.factor(severity_grades), color = as.factor(severity_grades))) +
  ggbeeswarm::geom_quasirandom(size=0.7, width=0.2,show.legend=F, alpha=0.6,varwidth=T,shape=16)+
  geom_boxplot(width=0.12,size=0.2,outlier.shape=NA,fill=NA, color = "black", alpha = 0.7)+
  geom_violin(size=0.2,trim = T,alpha=0.2)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("Healthy" = "#466d2d", "Non-Severe" = "#72bcc5", "Severe" = "#e1b941", "Fatal" = "#6d2d46"))+
  scale_color_manual(values = c("Healthy" = "#466d2d", "Non-Severe" = "#72bcc5", "Severe" = "#e1b941", "Fatal" = "#6d2d46"))+
  annotate("text", x = 1, y = 10, label = "JT p < 2.2e-16")+
  xlab(NULL)+
  ylab("Score")+
  scale_y_continuous(limits = c(-5, 10), breaks = c(-5, 0, 5, 10)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Lymphoid Dysregulation\nScore")

print(public_severity_lymphoid_score_plot)

#myeloid
rawr::jt.test(public_score_table$myeloid_score,as.factor(public_score_table$severity_grades))

public_severity_myeloid_score_plot <- ggplot(filter(public_score_table, !is.na(severity_grades)), aes(x = as.factor(severity_grades), y = myeloid_z_score, fill = as.factor(severity_grades), color = as.factor(severity_grades))) +
  ggbeeswarm::geom_quasirandom(size=0.7, width=0.2,show.legend=F, alpha=0.6,varwidth=T,shape=16)+
  geom_boxplot(width=0.12,size=0.2,outlier.shape=NA,fill=NA, color = "black", alpha = 0.7)+
  geom_violin(size=0.2,trim = T,alpha=0.2)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("Healthy" = "#466d2d", "Non-Severe" = "#72bcc5", "Severe" = "#e1b941", "Fatal" = "#6d2d46"))+
  scale_color_manual(values = c("Healthy" = "#466d2d", "Non-Severe" = "#72bcc5", "Severe" = "#e1b941", "Fatal" = "#6d2d46"))+
  annotate("text", x = 1, y = 10, label = "JT p < 2.2e-16")+
  xlab(NULL)+
  ylab("Score")+
  scale_y_continuous(limits = c(-5, 10), breaks = c(-5, 0, 5, 10)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Myeloid Dysregulation\nScore")

print(public_severity_myeloid_score_plot)

######Public Framework plots######
public_severity_dotplot <- ggplot(public_score_table, aes(x = myeloid_z_score, y = lymphoid_z_score, color = severity))+
  ylim(-5,10)+
  xlim(-5,12.5)+
  geom_rect(aes(xmin = -5, xmax = 1.65, ymin = -5, ymax = 1.65), fill = "#999999", alpha = 0.4, inherit.aes = F)+
  geom_rect(aes(xmin = -5, xmax = 1.65, ymin = 1.65, ymax = 10), fill = "#B0A4E3", alpha = 0.4, inherit.aes = F)+
  geom_rect(aes(xmin = 1.65, xmax = 12.5, ymin = -5, ymax = 1.65), fill = "lightskyblue", alpha = 0.4, inherit.aes = F)+
  geom_rect(aes(xmin = 1.65, xmax = 12.5, ymin = 1.65, ymax = 10), fill = "#F47B00", alpha = 0.4, inherit.aes = F)+
  geom_point(size = 0.15)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  theme_minimal() +
  scale_color_manual(values = c("Non-Severe" = "bisque", "Severe" = "black"))+
  geom_hline(yintercept = 1.65, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 1.65, linetype = "dashed", color = "black")+
  theme(legend.position = "none",axis.line = element_line(size = 1, arrow = arrow(type = "closed", length = unit(0.25, "inches"))), panel.border = element_blank(), panel.background = element_blank())+ 
  ylab("Lymphoid Dysregulation Score")+
  xlab("Myeloid Dysregulation Score")

public_severity_dotplot <- ggMarginal(public_severity_dotplot, groupFill = TRUE, groupColour = F)

print(public_severity_dotplot)

fisher_dysregulated_v_balanced <- fisher.test(public_score_table$subgroup == "balanced", public_score_table$severity)

prop_public_subgroup_severity <- public_score_table%>%
  group_by(subgroup)%>%
  summarize(prop_severe = mean(severity == "Severe", na.rm = TRUE)*100,
            n = n())

public_subgroup_severity_plot <- ggplot(prop_public_subgroup_severity, aes(x = subgroup, y = prop_severe, fill = subgroup)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_fill_manual(values = c("balanced" = "#999999", "lymphoid dysregulation" = "#B0A4E3", "myeloid dysregulation" = "lightskyblue", "hyperinflammatory" = "#F47B00"))+
  #adjust labels of x axis to capitalize words
  scale_x_discrete(labels = function(x) stringr::str_to_title(x)) +
  theme(axis.text.x = element_blank())+
  theme(legend.position = "none")+
  xlab("Subgroup") +
  ylab("Proportion Severe (%)")+
  ylim(0,60)+
  annotate("text", x = 2, y = 60, label = "OR 5.2 (3.9-7.0)\np < 2.2e-16")+
  geom_segment(aes(x = 1, xend = 2.9, y = 55, yend = 55), size = 0.2)+
  geom_segment(aes(x = 2, xend = 4, y = 53, yend = 53), size = 0.2)+
  geom_segment(aes(x = 1, xend = 1, y = 55, yend = 54), size = 0.2)+
  geom_segment(aes(x = 2.9, xend = 2.9, y = 55, yend = 54), size = 0.2)+
  geom_segment(aes(x = 2, xend = 2, y = 53, yend = 52), size = 0.2)+
  geom_segment(aes(x = 4, xend = 4, y = 53, yend = 52), size = 0.2)+
  geom_segment(aes(x = 3, xend = 3, y = 53, yend = 52), size = 0.2)


print(public_subgroup_severity_plot)

######SUBSPACE violin plots#####
#lymphoid
rawr::jt.test(filter(scores_w_pheno, !is.na(severity_grades),  timepoint == "baseline")$lymphoid_score,as.factor(filter(scores_w_pheno, !is.na(severity_grades),  timepoint == "baseline")$severity_grades))

subspace_severity_lymphoid_score_plot <- ggplot(filter(scores_w_pheno, !is.na(severity_grades),  timepoint == "baseline"), aes(x = as.factor(severity_grades), y = lymphoid_z_score, fill = as.factor(severity_grades), color = as.factor(severity_grades))) +
  ggbeeswarm::geom_quasirandom(size=0.7, width=0.2,show.legend=F, alpha=0.6,varwidth=T,shape=16)+
  geom_boxplot(width=0.12,size=0.2,outlier.shape=NA,fill=NA, color = "black", alpha = 0.7)+
  geom_violin(size=0.2,trim = T,alpha=0.2)+
  theme(legend.position = "none")+
  #renames X axes labels
  scale_fill_manual(values = c("Healthy" = "#466d2d", "Non-Severe" = "#72bcc5", "Severe" = "#e1b941", "Fatal" = "#6d2d46"))+
  scale_color_manual(values = c("Healthy" = "#466d2d", "Non-Severe" = "#72bcc5", "Severe" = "#e1b941", "Fatal" = "#6d2d46"))+
  annotate("text", x = 1, y = 10, label = "JT p < 2.2e-16")+
  xlab(NULL)+
  ylab("Score")+
  theme(axis.ticks.x = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Lymphoid Dysregulation\nScore")

print(subspace_severity_lymphoid_score_plot)

#Myeloid
rawr::jt.test(scores_w_pheno$myeloid_score,as.factor(scores_w_pheno$severity_grades))

subspace_severity_myeloid_score_plot <- ggplot(filter(scores_w_pheno, !is.na(severity_grades),  timepoint == "baseline"), aes(x = as.factor(severity_grades), y = myeloid_z_score, fill = as.factor(severity_grades), color = as.factor(severity_grades))) +
  ggbeeswarm::geom_quasirandom(size=0.7, width=0.2,show.legend=F, alpha=0.6,varwidth=T,shape=16)+
  geom_boxplot(width=0.12,size=0.2,outlier.shape=NA,fill=NA, color = "black", alpha = 0.7)+
  geom_violin(size=0.2,trim = T,alpha=0.2)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("Healthy" = "#466d2d", "Non-Severe" = "#72bcc5", "Severe" = "#e1b941", "Fatal" = "#6d2d46"))+
  scale_color_manual(values = c("Healthy" = "#466d2d", "Non-Severe" = "#72bcc5", "Severe" = "#e1b941", "Fatal" = "#6d2d46"))+
  annotate("text", x = 1, y = 7, label = "JT p < 2.2e-16")+
  xlab(NULL)+
  ylab("Score")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Myeloid Dysregulation\nScore")

print(subspace_severity_myeloid_score_plot)

#####Forest plots for mortality########
#lymphoid
stanford_log_model_lymphoid <- glm(d30_mort ~ lymphoid_score, data = stanford, family = "binomial")
summary(stanford_log_model_lymphoid)

amsterdam_log_model_lymphoid <- glm(d30_mort ~ lymphoid_score, data = filter(amsterdam, timepoint == "baseline"), family = "binomial")
summary(amsterdam_log_model_lymphoid)

victas_log_model_lymphoid <- glm(d30_mort ~ lymphoid_score, data = victas, family = "binomial")
summary(victas_log_model_lymphoid)

savemore_log_model_lymphoid <- glm(d30_mort ~ lymphoid_score, data = filter(savemore, timepoint == "baseline"), family = "binomial")
summary(savemore_log_model_lymphoid)

ufl_log_model_lymphoid <- glm(d30_mort ~ lymphoid_score, data = ufl, family = "binomial")
summary(ufl_log_model_lymphoid)

cchmc_log_model_lymphoid <- glm(d30_mort ~ lymphoid_score, data = filter(cchmc, timepoint == "baseline"), family = "binomial")
summary(cchmc_log_model_lymphoid)

charles_log_model_lymphoid <- glm(d30_mort ~ lymphoid_score, data = charles, family = "binomial")
summary(charles_log_model_lymphoid)

acutelines_log_model_lymphoid <- glm(d30_mort ~ lymphoid_score, data = acutelines, family = "binomial")
summary(acutelines_log_model_lymphoid)

summary_log_model_lymphoid <- glm(d30_mort ~ lymphoid_score, data = filter(scores_w_pheno, timepoint == "baseline"), family = "binomial")
summary(summary_log_model_lymphoid)
exp(coef(summary_log_model_lymphoid))
1/confint(summary_log_model_lymphoid)

log_model_lymphoids <- list(Acutelines = acutelines_log_model_lymphoid, Amsterdam = amsterdam_log_model_lymphoid, CCHMC = cchmc_log_model_lymphoid,Charles =  charles_log_model_lymphoid,Savemore = savemore_log_model_lymphoid, Stanford = stanford_log_model_lymphoid, UFL = ufl_log_model_lymphoid,VICTAS = victas_log_model_lymphoid, Summary = summary_log_model_lymphoid)



coefs <- sapply(log_model_lymphoids, function(mod) summary(mod)$coefficients["lymphoid_score", 1])

cis <- t(sapply(log_model_lymphoids, function(mod) confint(mod)["lymphoid_score", ]))

coefs
cis

forest_plot_lymphoid_score <- forestplot(labeltext = names(log_model_lymphoids),
                                         
                                         mean = coefs,
                                         
                                         lower = cis[, 1],
                                         
                                         upper = cis[, 2],
                                         
                                         is.summary = c(rep(FALSE, length(log_model_lymphoids) - 1), TRUE),
                                         
                                         title = "Lymphoid Dysregulation Score",
                                         
                                         xlab = "Log(Odds Ratio): 30-Day Mortality",
                                         
                                         col = fpColors(box = c("black", "black"), line = c("black", "black"), summary = c("black", "black"))
                                         
)

print(forest_plot_lymphoid_score)


#Myeloid

stanford_log_model_myeloid <- glm(d30_mort ~ myeloid_score, data = stanford, family = "binomial")
summary(stanford_log_model_myeloid)

amsterdam_log_model_myeloid <- glm(d30_mort ~ myeloid_score, data = filter(amsterdam, timepoint == "baseline"), family = "binomial")
summary(amsterdam_log_model_myeloid)

victas_log_model_myeloid <- glm(d30_mort ~ myeloid_score, data = victas, family = "binomial")
summary(victas_log_model_myeloid)

savemore_log_model_myeloid <- glm(d30_mort ~ myeloid_score, data = filter(savemore, timepoint == "baseline"), family = "binomial")
summary(savemore_log_model_myeloid)

ufl_log_model_myeloid <- glm(d30_mort ~ myeloid_score, data = ufl, family = "binomial")
summary(ufl_log_model_myeloid)

cchmc_log_model_myeloid <- glm(d30_mort ~ myeloid_score, data = filter(cchmc, timepoint == "baseline"), family = "binomial")
summary(cchmc_log_model_myeloid)

charles_log_model_myeloid <- glm(d30_mort ~ myeloid_score, data = charles, family = "binomial")
summary(charles_log_model_myeloid)

acutelines_log_model_myeloid <- glm(d30_mort ~ myeloid_score, data = acutelines, family = "binomial")
summary(acutelines_log_model_myeloid)

summary_log_model_myeloid <- glm(d30_mort ~ myeloid_score, data = filter(scores_w_pheno, timepoint == "baseline"), family = "binomial")
summary(summary_log_model_myeloid)
exp(coef(summary_log_model_myeloid))
1/confint(summary_log_model_myeloid)



log_model_myeloids <- list(Acutelines = acutelines_log_model_myeloid, Amsterdam = amsterdam_log_model_myeloid, CCHMC = cchmc_log_model_myeloid, Charles = charles_log_model_myeloid ,Savemore = savemore_log_model_myeloid, Stanford = stanford_log_model_myeloid, UFL = ufl_log_model_myeloid,VICTAS = victas_log_model_myeloid, Summary = summary_log_model_myeloid)



coefs <- sapply(log_model_myeloids, function(mod) summary(mod)$coefficients["myeloid_score", 1])

cis <- t(sapply(log_model_myeloids, function(mod) confint(mod)["myeloid_score", ]))

coefs
cis

forest_plot_myeloid_score <- forestplot(labeltext = names(log_model_myeloids),
                                        
                                        mean = coefs,
                                        
                                        lower = cis[, 1],
                                        
                                        upper = cis[, 2],
                                        
                                        is.summary = c(rep(FALSE, length(log_model_myeloids) - 1), TRUE),
                                        
                                        title = "Myeloid Dysregulation Score",
                                        
                                        xlab = "Log(Odds Ratio): 30-Day Mortality",
                                        
                                        col = fpColors(box = c("black", "black"), line = c("black", "black"), summary = c("black", "black"))
                                        
)

print(forest_plot_myeloid_score)


#####SUBSPACE Framework plots for Mortality#####
all_dotplot_mortality <- ggplot(filter(scores_w_pheno, !is.na(d30_mort),  timepoint == "baseline")[order(scores_w_pheno$d30_mort),], aes(x = myeloid_z_score, y = lymphoid_z_score, color = as.factor(d30_mort))) +
  geom_rect(aes(xmin = -5, xmax = 1.65, ymin = -5, ymax = 1.65), fill = "#999999", alpha = 0.4, inherit.aes = F)+
  geom_rect(aes(xmin = -5, xmax = 1.65, ymin = 1.65, ymax = 12), fill = "#B0A4E3", alpha = 0.4, inherit.aes = F)+
  geom_rect(aes(xmin = 1.65, xmax = 10, ymin = -5, ymax = 1.65), fill = "lightskyblue", alpha = 0.4, inherit.aes = F)+
  geom_rect(aes(xmin = 1.65, xmax = 10, ymin = 1.65, ymax = 12), fill = "#F47B00", alpha = 0.4, inherit.aes = F)+
  geom_point(size = 0.15) +
  theme_minimal() +
  ylim(-5,10)+
  xlim(-5,7.5)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+  
  scale_color_manual(values = c("0" = "bisque", "1" = "black"))+
  geom_vline(xintercept = 1.65, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.65, linetype = "dashed", color = "black") +
  theme(legend.position = "none",axis.line = element_line(size = 1, arrow = arrow(type = "closed", length = unit(0.25, "inches"))), panel.border = element_blank(), panel.background = element_blank())+ 
  xlab("Myeloid Dysregulation Score") +
  ylab("Lymphoid Dysregulation Score")

all_dotplot_mortality <- ggMarginal(all_dotplot_mortality, groupFill = TRUE, groupColour = F)

print(all_dotplot_mortality)

fisher_dysregulated_v_balanced <- fisher.test(filter(scores_w_pheno,  timepoint == "baseline")$subgroup == "balanced", filter(scores_w_pheno,  timepoint == "baseline")$d30_mort)
fisher_dysregulated_v_balanced
1/fisher_dysregulated_v_balanced$estimate
1/fisher_dysregulated_v_balanced$conf.int

mort_prop_all_subgroup <- scores_w_pheno%>%
  filter(timepoint == "baseline",  !is.na(d30_mort))%>%
  group_by(subgroup)%>%
  summarize(mortality_prop = mean(d30_mort == 1, na.rm = TRUE)*100,n = n())

mort_prop_all_subgroup_barplot <- ggplot(mort_prop_all_subgroup, aes(x = subgroup, y = mortality_prop, fill = subgroup)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_fill_manual(values = c("balanced" = "#999999", "lymphoid dysregulation" = "#B0A4E3", "myeloid dysregulation" = "lightskyblue", "hyperinflammatory" = "#F47B00"))+
  scale_x_discrete(labels = function(x) stringr::str_to_title(x)) +
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank())+
  xlab("Subgroup") +
  ylab("30-day Mortality (%)")+
  ylim(0, 28)  +
  annotate("text", x = 2, y = 27, label = "OR 3.4 (2.3-5.4) \n p = 4.1e-11")+
  geom_segment(aes(x = 1, xend = 2.9, y = 23, yend = 23), size = 0.2)+
  geom_segment(aes(x = 2, xend = 4, y = 21, yend = 21), size = 0.2)+
  geom_segment(aes(x = 1, xend = 1, y = 23, yend = 22), size = 0.2)+
  geom_segment(aes(x = 2.9, xend = 2.9, y = 23, yend = 22), size = 0.2)+
  geom_segment(aes(x = 2, xend = 2, y = 21, yend = 20), size = 0.2)+
  geom_segment(aes(x = 4, xend = 4, y = 21, yend = 20), size = 0.2)+
  geom_segment(aes(x = 3, xend = 3, y = 21, yend = 20), size = 0.2)

print(mort_prop_all_subgroup_barplot)

mortality_plots <- plot_grid(all_dotplot_mortality, mort_prop_all_subgroup_barplot, labels = c("A","B"))

print(mortality_plots)

########Subspace severity plots#########
all_dotplot_severity <- ggplot(filter(scores_w_pheno, !is.na(severity), condition != "healthy",  timepoint == "baseline")[order(scores_w_pheno$severity),], aes(x = myeloid_z_score, y = lymphoid_z_score, color = as.factor(severity))) +
  geom_rect(aes(xmin = -5, xmax = 1.65, ymin = -5, ymax = 1.65), fill = "#999999", alpha = 0.4, inherit.aes = F)+
  geom_rect(aes(xmin = -5, xmax = 1.65, ymin = 1.65, ymax = 12), fill = "#B0A4E3", alpha = 0.4, inherit.aes = F)+
  geom_rect(aes(xmin = 1.65, xmax = 10, ymin = -5, ymax = 1.65), fill = "lightskyblue", alpha = 0.4, inherit.aes = F)+
  geom_rect(aes(xmin = 1.65, xmax = 10, ymin = 1.65, ymax = 12), fill = "#F47B00", alpha = 0.4, inherit.aes = F)+
  geom_point(size = 0.15) +
  theme_minimal() +
  ylim(-5,10)+
  xlim(-5,7.5)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+  
  scale_color_manual(values = c("non-severe" = "bisque", "severe" = "black"))+
  geom_vline(xintercept = 1.65, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.65, linetype = "dashed", color = "black") +
  theme(legend.position = "none",axis.line = element_line(size = 1, arrow = arrow(type = "closed", length = unit(0.25, "inches"))), panel.border = element_blank(), panel.background = element_blank())+ 
  xlab("Myeloid Dysregulation Score") +
  ylab("Lymphoid Dysregulation Score")

all_dotplot_severity <- ggMarginal(all_dotplot_severity, groupFill = TRUE, groupColour = F)

print(all_dotplot_severity)

fisher_dysregulated_v_balanced <- fisher.test(filter(scores_w_pheno, !is.na(severity), condition != "healthy",  timepoint == "baseline")$subgroup == "balanced", filter(scores_w_pheno,!is.na(severity), condition != "healthy",  timepoint == "baseline")$severity)
fisher_dysregulated_v_balanced
1/fisher_dysregulated_v_balanced$estimate
1/fisher_dysregulated_v_balanced$conf.int

severity_prop_all_subgroup <- scores_w_pheno%>%
  filter(timepoint == "baseline",  !is.na(severity), condition != "healthy")%>%
  group_by(subgroup)%>%
  summarize(severity_prop = mean(severity == "severe", na.rm = TRUE)*100,
            n = n()
  )


severity_prop_all_subgroup_barplot <- ggplot(severity_prop_all_subgroup, aes(x = subgroup, y = severity_prop, fill = subgroup)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_fill_manual(values = c("balanced" = "#999999", "lymphoid dysregulation" = "#B0A4E3", "myeloid dysregulation" = "lightskyblue", "hyperinflammatory" = "#F47B00"))+
  #adjust labels of x axis to capitalize words
  scale_x_discrete(labels = function(x) stringr::str_to_title(x)) +
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank())+
  xlab("Subgroup") +
  ylab("Proportion Severe (%)")+
  ylim(0, 90)  +
  annotate("text", x = 2, y = 85, label = "p < 2.2e-16")+
  annotate("text", x = 2, y = 88, label = "OR 7.1 (5.6 - 8.9)")+
  geom_segment(aes(x = 1, xend = 2.9, y = 83, yend = 83), size = 0.2)+
  geom_segment(aes(x = 2, xend = 4, y = 81, yend = 81), size = 0.3)+
  geom_segment(aes(x = 1, xend = 1, y = 83, yend = 82), size = 0.2)+
  geom_segment(aes(x = 2.9, xend = 2.9, y = 83, yend = 82), size = 0.2)+
  geom_segment(aes(x = 2, xend = 2, y = 81, yend = 80), size = 0.2)+
  geom_segment(aes(x = 4, xend = 4, y = 81, yend = 80), size = 0.2)+
  geom_segment(aes(x = 3, xend = 3, y = 81, yend = 80), size = 0.2)

print(severity_prop_all_subgroup_barplot)


framework_plot <- ggplot() +
  xlim(0,4)+
  ylim(0,4)+
  xlab("Myeloid Dysregulation Score") +
  ylab("Lymphoid Dysregulation Score") +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.line = element_line(size = 1, arrow = arrow(type = "closed", length = unit(0.25, "inches"))),axis.ticks = element_blank(), panel.border = element_blank(), panel.background = element_blank())+
  geom_rect(aes(xmin = 0, xmax = 2, ymin = 0, ymax = 2), fill = "#999999") +
  geom_rect(aes(xmin = 0, xmax = 2, ymin = 2, ymax = 4), fill = "#B0A4E3") +
  geom_rect(aes(xmin = 2, xmax = 4, ymin = 0, ymax = 2), fill = "lightskyblue") +
  geom_rect(aes(xmin = 2, xmax = 4, ymin = 2, ymax = 4), fill = "#F47B00") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  geom_hline(yintercept = 2, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 2, linetype = "dashed", color = "black") +
  annotate("text", x = 1, y = 1, label = "Balanced") +
  annotate("text", x = 3, y = 1, label = "Myeloid\nDysregulation") +
  annotate("text", x = 1, y = 3, label = "Lymphoid\nDysregulation") +
  annotate("text", x = 3, y = 3, label = "System-Wide\nDysregulation") 


print(framework_plot)


jt_plot <- plot_grid(public_severity_myeloid_score_plot,public_severity_lymphoid_score_plot,
                     nrow = 1, align = "hv", axis = 'tblr', labels = c("A","B"),label_size = 12, label_fontfamily = base_font_family, label_colour = black_text_colour)

print(jt_plot)

public_row_1 <- plot_grid(jt_plot,framework_plot, nrow =1 , ncol = 2, align = "hv", axis = 'tblr', labels = c("","C"), label_size = 12, label_fontfamily = base_font_family, label_colour = black_text_colour)

print(public_row_1)


public_row_2 <- plot_grid(public_severity_dotplot,public_subgroup_severity_plot,
                          nrow = 1, align = "hv", axis = 'tblr', labels = c("D","E"), label_size = 12, label_fontfamily = base_font_family, label_colour = black_text_colour)

print(public_row_2)


subspace_violins <- plot_grid(subspace_severity_myeloid_score_plot,subspace_severity_lymphoid_score_plot,
                              nrow = 1, align = "hv", axis = 'tblr', labels = c("F","G"), label_size = 12, label_fontfamily = base_font_family, label_colour = black_text_colour)

print(subspace_violins)

subspace_forest_plots <- plot_grid(forest_plot_myeloid_score, forest_plot_lymphoid_score,
                                   nrow = 1, align = "hv", axis = 'tblr', labels = c("H","I"), label_size = 12, label_fontfamily = base_font_family, label_colour = black_text_colour)

print(subspace_forest_plots)
#Forest plots can't be joined this way. Need to be added afterwards


subspace_top_row <- plot_grid(subspace_violins, subspace_forest_plots,
                              nrow = 1, align = "hv", axis = 'tblr', labels = c("",""), label_size = 12, label_fontfamily = base_font_family, label_colour = black_text_colour)

subspace_severity_plot <- plot_grid(all_dotplot_severity, severity_prop_all_subgroup_barplot,
                                     nrow = 1, align = "hv", axis = 'tblr', labels = c("J", "K"), label_size = 12, label_fontfamily = base_font_family, label_colour = black_text_colour)

print(subspace_severity_plot)


#combine all the plots
public_title <- ggdraw() + draw_label("Framework Applied to Public Data", size = base_text_size, colour = black_text_colour, fontfamily = base_font_family) +
  theme(plot.background = element_rect(fill = "#DCDCDC", colour = NA))

subspace_title <- ggdraw() + draw_label("Framework Applied to SUBSPACE Data", size = base_text_size, colour = black_text_colour, fontfamily = base_font_family) +
  theme(plot.background = element_rect(fill = "#DCDCDC", colour = NA))

subspace_fig4 <- plot_grid(public_title, public_row_1,public_row_2, subspace_title,subspace_top_row,subspace_severity_plot,
                           nrow = 6, align = "hv", axis = 'tblr', label_size = 12, label_fontfamily = base_font_family, label_colour = black_text_colour, rel_heights = c(.1,1,1, 0.1, 1,1))
#Forest plots can't be joined this way. Need to be added afterwards

print(subspace_fig4)

