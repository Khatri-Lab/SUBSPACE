source("https://raw.githubusercontent.com/Khatri-Lab/SUBSPACE/main/subspace_load_data.R")


########Differences in scores by severity on framework########
subspace_summary_stats <- scores_w_pheno %>%
  filter(timepoint == "baseline", !is.na(severity_grades)) %>%
  group_by(severity_grades) %>%
  summarise(
    myeloid_mean = mean(myeloid_z_score, na.rm = TRUE),
    myeloid_sd = sd(myeloid_z_score, na.rm = TRUE),
    lymphoid_mean = mean(lymphoid_z_score, na.rm = TRUE),
    lymphoid_sd = sd(lymphoid_z_score, na.rm = TRUE)
  )

subspace_threshold_dotplot <- ggplot(filter(scores_w_pheno, timepoint == "baseline", !is.na(severity_grades)), aes(x = myeloid_z_score, y = lymphoid_z_score, color = severity_grades))+
  geom_ellipse(inherit.aes = F,data = filter(subspace_summary_stats, severity_grades == "Healthy"), aes(x0 = myeloid_mean, y0 = lymphoid_mean, a = myeloid_sd, b = lymphoid_sd, angle = 0), color = "#7267AE", fill = "#7267AE", alpha = 0.2) +
  geom_ellipse(inherit.aes = F,data = filter(subspace_summary_stats, severity_grades == "Non-Severe"), aes(x0 = myeloid_mean, y0 = lymphoid_mean, a = myeloid_sd, b = lymphoid_sd, angle = 0), color = "#72bcc5",fill = "#72bcc5", alpha = 0.2) +
  geom_ellipse(inherit.aes = F,data = filter(subspace_summary_stats, severity_grades == "Severe"), aes(x0 = myeloid_mean, y0 = lymphoid_mean, a = myeloid_sd, b = lymphoid_sd, angle = 0), color = "#e1b941", fill = "#e1b941", alpha = 0.2) +
  geom_ellipse(inherit.aes = F,data = filter(subspace_summary_stats, severity_grades == "Fatal"), aes(x0 = myeloid_mean, y0 = lymphoid_mean, a = myeloid_sd, b = lymphoid_sd, angle = 0), color = "#6d2d46", fill = "#6d2d46", alpha = 0.2) +
  geom_point(size = 0.15)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  theme_minimal() +
  scale_color_manual(values = c("Healthy" = "#7267AE", "Non-Severe" = "#72bcc5", "Severe" = "#e1b941", "Fatal" = "#6d2d46"))+
  geom_hline(yintercept = 1.65, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 1.65, linetype = "dashed", color = "black")+
  theme(legend.position = "right",axis.line = element_line(size = 1, arrow = arrow(type = "closed", length = unit(0.25, "inches"))), panel.border = element_blank(), panel.background = element_blank())+ 
  ylab("Lymphoid Dysregulation Score")+
  xlab("Myeloid Dysregulation Score")



subspace_threshold_dotplot


scores_w_pheno_baseline <- scores_w_pheno %>%
  filter(timepoint == "baseline", !is.na(severity_grades))%>%
  mutate(myeloid_quintile = ntile(myeloid_z_score, 5),
         lymphoid_quintile = ntile(lymphoid_z_score, 5)) 

subspace_myeloid_quintile_proportions <- scores_w_pheno_baseline %>%
  group_by(myeloid_quintile, severity_grades) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(myeloid_quintile) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

subspace_myeloid_severity_bar_plot <- ggplot(subspace_myeloid_quintile_proportions, aes(x = factor(myeloid_quintile), y = proportion, fill = severity_grades)) +
  geom_bar(stat = "identity") +
  labs(x = "Quintile of Myeloid Score", y = "Proportion", fill = "Severity Grades") +
  theme_minimal() +
  scale_fill_manual(values = c("Healthy" = "#7267AE", "Non-Severe" = "#72bcc5", "Severe" = "#e1b941", "Fatal" = "#6d2d46")) +
  theme(legend.position = "none")


print(subspace_myeloid_severity_bar_plot)


subspace_lymphoid_quintile_proportions <- scores_w_pheno_baseline %>%
  group_by(lymphoid_quintile, severity_grades) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(lymphoid_quintile) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

subspace_lymphoid_severity_bar_plot <- ggplot(subspace_lymphoid_quintile_proportions, aes(x = factor(lymphoid_quintile), y = proportion, fill = severity_grades)) +
  geom_bar(stat = "identity") +
  labs(x = "Quintile of Lymphoid Score", y = "Proportion", fill = "Severity Grades") +
  theme_minimal() +
  scale_fill_manual(values = c("Healthy" = "#7267AE", "Non-Severe" = "#72bcc5", "Severe" = "#e1b941", "Fatal" = "#6d2d46")) +
  theme(legend.position = "none")

print(subspace_lymphoid_severity_bar_plot)


#######Sensitivity and Specificity for myeloid and Lymphoid########
#Myeloid
scores_w_pheno_baseline <- scores_w_pheno_baseline %>%
  mutate(outcome_mild = ifelse(severity_grades %in% c("Healthy"), 0, 1),
         outcome_severe = ifelse(severity_grades %in% c("Severe", "Fatal"), 1, 0),
         outcome_fatal = ifelse(severity_grades %in% c("Fatal"), 1, 0))

# Create a sequence of cut-off values
cutoffs <- seq(min(scores_w_pheno_baseline$myeloid_z_score, na.rm = TRUE), 
               max(scores_w_pheno_baseline$myeloid_z_score, na.rm = TRUE), 
               length.out = 100)

# Initialize a data frame to store results
results_mild <- data.frame()

# Calculate sensitivity and specificity for each cutoff and severity grade
for (cutoff in cutoffs) {
  predictions <- ifelse(scores_w_pheno_baseline$myeloid_z_score >= cutoff, 1, 0)
  
  # Calculate sensitivity
  true_positives <- sum(predictions[scores_w_pheno_baseline$outcome_mild == 1] == 1)
  total_positives <- sum(scores_w_pheno_baseline$outcome_mild == 1)
  sensitivity_mild <- ifelse(total_positives > 0, true_positives / total_positives, NA)  # Use NA if no positives
  
  # Calculate specificity
  true_negatives <- sum(predictions[scores_w_pheno_baseline$outcome_mild == 0] == 0)
  total_negatives <- sum(scores_w_pheno_baseline$outcome_mild == 0)
  specificity_mild <- ifelse(total_negatives > 0, true_negatives / total_negatives, NA)  # Use NA if no negatives
  
  # Store results
  results_mild <- rbind(results_mild, data.frame(cutoff = cutoff, sensitivity_mild = sensitivity_mild, specificity_mild = specificity_mild))
}


# Initialize a data frame to store results
results_severe <- data.frame()

# Calculate sensitivity and specificity for each cutoff and severity grade
for (cutoff in cutoffs) {
  predictions <- ifelse(scores_w_pheno_baseline$myeloid_z_score >= cutoff, 1, 0)
  
  # Calculate sensitivity
  true_positives <- sum(predictions[scores_w_pheno_baseline$outcome_severe == 1 ] == 1)
  total_positives <- sum(scores_w_pheno_baseline$outcome_severe ==1)
  sensitivity_severe <- ifelse(total_positives > 0, true_positives / total_positives, NA)  # Use NA if no positives
  
  # Calculate specificity
  true_negatives <- sum(predictions[scores_w_pheno_baseline$outcome_severe == 0] == 0)
  total_negatives <- sum(scores_w_pheno_baseline$outcome_severe == 0)
  specificity_severe <- ifelse(total_negatives > 0, true_negatives / total_negatives, NA)  # Use NA if no negatives
  
  # Store results
  results_severe <- rbind(results_severe, data.frame(cutoff = cutoff, sensitivity_severe = sensitivity_severe, specificity_severe = specificity_severe))
}


# Initialize a data frame to store results
results_fatal <- data.frame()

# Calculate sensitivity and specificity for each cutoff and severity grade
for (cutoff in cutoffs) {
  predictions <- ifelse(scores_w_pheno_baseline$myeloid_z_score >= cutoff, 1, 0)
  
  # Calculate sensitivity
  true_positives <- sum(predictions[scores_w_pheno_baseline$outcome_fatal == 1 ] == 1)
  total_positives <- sum(scores_w_pheno_baseline$outcome_fatal ==1)
  sensitivity_fatal <- ifelse(total_positives > 0, true_positives / total_positives, NA)  # Use NA if no positives
  
  # Calculate specificity
  true_negatives <- sum(predictions[scores_w_pheno_baseline$outcome_fatal == 0] == 0)
  total_negatives <- sum(scores_w_pheno_baseline$outcome_fatal == 0)
  specificity_fatal <- ifelse(total_negatives > 0, true_negatives / total_negatives, NA)  # Use NA if no negatives
  
  # Store results
  results_fatal <- rbind(results_fatal, data.frame(cutoff = cutoff, sensitivity_fatal = sensitivity_fatal, specificity_fatal = specificity_fatal))
}


sensitivity_specificity_table <- left_join(results_mild, results_severe, by = c("cutoff")) %>%
  left_join(.,results_fatal, by = c("cutoff"))


myeloid_sensitivity_specificity_plot <- ggplot(sensitivity_specificity_table, aes(x = cutoff)) +
  geom_line(aes(y = sensitivity_mild, color = "Mild")) +
  geom_line(aes(y = sensitivity_severe, color = "Severe")) +
  geom_line(aes(y = sensitivity_fatal, color = "Fatal")) +
  geom_line(aes(y = specificity_mild, color = "Mild"), linetype = "dashed") +
  geom_line(aes(y = specificity_severe, color = "Severe"), linetype = "dashed") +
  geom_line(aes(y = specificity_fatal, color = "Fatal"), linetype = "dashed") +
  labs(x = "Myeloid Dysregulation Score", y = "Sensitivity/Specificity", color = "Outcome") +
  scale_color_manual(values = c("Mild" = "#72bcc5", "Severe" = "#e1b941", "Fatal" = "#6d2d46")) +
  theme_minimal() +
  xlim(-2.5, 6) +
  theme(legend.position = "none")

myeloid_sensitivity_specificity_plot


# Lymphoid
scores_w_pheno_baseline <- scores_w_pheno_baseline %>%
  mutate(outcome_mild = ifelse(severity_grades %in% c("Healthy"), 0, 1),
         outcome_severe = ifelse(severity_grades %in% c("Severe", "Fatal"), 1, 0),
         outcome_fatal = ifelse(severity_grades %in% c("Fatal"), 1, 0))

# Create a sequence of cut-off values
cutoffs <- seq(min(scores_w_pheno_baseline$lymphoid_z_score, na.rm = TRUE), 
               max(scores_w_pheno_baseline$lymphoid_z_score, na.rm = TRUE), 
               length.out = 100)

# Initialize a data frame to store results
results_mild <- data.frame()

# Calculate sensitivity and specificity for each cutoff and severity grade
for (cutoff in cutoffs) {
  predictions <- ifelse(scores_w_pheno_baseline$lymphoid_z_score >= cutoff, 1, 0)
  
  # Calculate sensitivity
  true_positives <- sum(predictions[scores_w_pheno_baseline$outcome_mild == 1] == 1)
  total_positives <- sum(scores_w_pheno_baseline$outcome_mild == 1)
  sensitivity_mild <- ifelse(total_positives > 0, true_positives / total_positives, NA)  # Use NA if no positives
  
  # Calculate specificity
  true_negatives <- sum(predictions[scores_w_pheno_baseline$outcome_mild == 0] == 0)
  total_negatives <- sum(scores_w_pheno_baseline$outcome_mild == 0)
  specificity_mild <- ifelse(total_negatives > 0, true_negatives / total_negatives, NA)  # Use NA if no negatives
  
  # Store results
  results_mild <- rbind(results_mild, data.frame(cutoff = cutoff, sensitivity_mild = sensitivity_mild, specificity_mild = specificity_mild))
}


# Initialize a data frame to store results
results_severe <- data.frame()

# Calculate sensitivity and specificity for each cutoff and severity grade
for (cutoff in cutoffs) {
  predictions <- ifelse(scores_w_pheno_baseline$lymphoid_z_score >= cutoff, 1, 0)
  
  # Calculate sensitivity
  true_positives <- sum(predictions[scores_w_pheno_baseline$outcome_severe == 1 ] == 1)
  total_positives <- sum(scores_w_pheno_baseline$outcome_severe ==1)
  sensitivity_severe <- ifelse(total_positives > 0, true_positives / total_positives, NA)  # Use NA if no positives
  
  # Calculate specificity
  true_negatives <- sum(predictions[scores_w_pheno_baseline$outcome_severe == 0] == 0)
  total_negatives <- sum(scores_w_pheno_baseline$outcome_severe == 0)
  specificity_severe <- ifelse(total_negatives > 0, true_negatives / total_negatives, NA)  # Use NA if no negatives
  
  # Store results
  results_severe <- rbind(results_severe, data.frame(cutoff = cutoff, sensitivity_severe = sensitivity_severe, specificity_severe = specificity_severe))
}


# Initialize a data frame to store results
results_fatal <- data.frame()

# Calculate sensitivity and specificity for each cutoff and severity grade
for (cutoff in cutoffs) {
  predictions <- ifelse(scores_w_pheno_baseline$lymphoid_z_score >= cutoff, 1, 0)
  
  # Calculate sensitivity
  true_positives <- sum(predictions[scores_w_pheno_baseline$outcome_fatal == 1 ] == 1)
  total_positives <- sum(scores_w_pheno_baseline$outcome_fatal ==1)
  sensitivity_fatal <- ifelse(total_positives > 0, true_positives / total_positives, NA)  # Use NA if no positives
  
  # Calculate specificity
  true_negatives <- sum(predictions[scores_w_pheno_baseline$outcome_fatal == 0] == 0)
  total_negatives <- sum(scores_w_pheno_baseline$outcome_fatal == 0)
  specificity_fatal <- ifelse(total_negatives > 0, true_negatives / total_negatives, NA)  # Use NA if no negatives
  
  # Store results
  results_fatal <- rbind(results_fatal, data.frame(cutoff = cutoff, sensitivity_fatal = sensitivity_fatal, specificity_fatal = specificity_fatal))
}


sensitivity_specificity_table <- left_join(results_mild, results_severe, by = c("cutoff")) %>%
  left_join(.,results_fatal, by = c("cutoff"))

lymphoid_sensitivity_specificity_plot <- ggplot(sensitivity_specificity_table, aes(x = cutoff)) +
  geom_line(aes(y = sensitivity_mild, color = "Mild")) +
  geom_line(aes(y = sensitivity_severe, color = "Severe")) +
  geom_line(aes(y = sensitivity_fatal, color = "Fatal")) +
  geom_line(aes(y = specificity_mild, color = "Mild"), linetype = "dashed") +
  geom_line(aes(y = specificity_severe, color = "Severe"), linetype = "dashed") +
  geom_line(aes(y = specificity_fatal, color = "Fatal"), linetype = "dashed") +
  labs(x = "Lymphoid Dysregulation Score", y = "Sensitivity/Specificity", color = "Outcome") +
  scale_color_manual(values = c("Mild" = "#72bcc5", "Severe" = "#e1b941", "Fatal" = "#6d2d46")) +
  theme_minimal() +
  xlim(-2.5, 6) +
  theme(legend.position = "none")

lymphoid_sensitivity_specificity_plot


row1 <- plot_grid("",subspace_threshold_dotplot,"",labels = c("","A",""), rel_widths = c(0.2,1,0.2), ncol = 3)

row2 <- plot_grid(subspace_myeloid_severity_bar_plot,subspace_lymphoid_severity_bar_plot, labels = c("B", "C"), ncol = 2)

row3 <- plot_grid(myeloid_sensitivity_specificity_plot, lymphoid_sensitivity_specificity_plot, labels = c("D", "E"), ncol = 2)

subspace_cutoff_plots <- plot_grid( row1, row2, row3, ncol = 1)
print(subspace_cutoff_plots)



########SUBSPACE results by infection type############
intersect_cols <- intersect(colnames(scores_w_pheno), colnames(glue))
scores_w_pheno_w_glue <- rbind(scores_w_pheno[, intersect_cols], glue[, intersect_cols]) 

subspace_summary_stats_infection <- scores_w_pheno_w_glue %>%
  filter(timepoint == "baseline", severity == "severe", infection_type == "Bacterial" | infection_type == "Viral"| infection_type == "non-infected") %>%
  group_by(infection_type) %>%
  summarise(
    myeloid_mean = mean(myeloid_z_score, na.rm = TRUE),
    myeloid_sd = sd(myeloid_z_score, na.rm = TRUE),
    lymphoid_mean = mean(lymphoid_z_score, na.rm = TRUE),
    lymphoid_sd = sd(lymphoid_z_score, na.rm = TRUE)
  )

subspace_threshold_dotplot_infection <- ggplot(filter(scores_w_pheno_w_glue,timepoint == "baseline", severity == "severe", infection_type == "Bacterial" | infection_type == "Viral" | infection_type == "non-infected"), aes(x = myeloid_z_score, y = lymphoid_z_score, color = infection_type))+
  geom_ellipse(inherit.aes = F,data = filter(subspace_summary_stats_infection, infection_type == "Viral"), aes(x0 = myeloid_mean, y0 = lymphoid_mean, a = myeloid_sd, b = lymphoid_sd, angle = 0), color = "mediumpurple2", fill = "mediumpurple2", alpha = 0.2) +
  geom_ellipse(inherit.aes = F,data = filter(subspace_summary_stats_infection, infection_type == "Bacterial"), aes(x0 = myeloid_mean, y0 = lymphoid_mean, a = myeloid_sd, b = lymphoid_sd, angle = 0), color = "goldenrod2",fill = "goldenrod2", alpha = 0.2) +
  geom_ellipse(inherit.aes = F,data = filter(subspace_summary_stats_infection, infection_type == "non-infected"), aes(x0 = myeloid_mean, y0 = lymphoid_mean, a = myeloid_sd, b = lymphoid_sd, angle = 0), color = "palevioletred",fill = "palevioletred", alpha = 0.2) +
  geom_point(size = 0.15)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  theme_minimal() +
  scale_color_manual(values = c("non-infected" = "palevioletred", "Viral" = "mediumpurple2", "Bacterial" = "goldenrod2"))+
  theme(legend.position = "right",axis.line = element_line(size = 1, arrow = arrow(type = "closed", length = unit(0.25, "inches"))), panel.border = element_blank(), panel.background = element_blank())+ 
  ylab("Lymphoid Dysregulation Score")+
  xlab("Myeloid Dysregulation Score")


subspace_threshold_dotplot_infection



########Sensitivity and specificity plots by infection type##############
#Myeloid
scores_w_pheno_baseline_viral <- scores_w_pheno_baseline %>%
  filter(infection_type == "Viral" | condition == "healthy") %>%
  mutate(outcome_mild = ifelse(severity_grades %in% c("Healthy"), 0, 1),
         outcome_severe = ifelse(severity_grades %in% c("Severe", "Fatal"), 1, 0),
         outcome_fatal = ifelse(severity_grades %in% c("Fatal"), 1, 0))

# Create a sequence of cut-off values
cutoffs <- seq(min(scores_w_pheno_baseline_viral$myeloid_z_score, na.rm = TRUE), 
               max(scores_w_pheno_baseline_viral$myeloid_z_score, na.rm = TRUE), 
               length.out = 100)

# Initialize a data frame to store results
results_mild_viral <- data.frame()

# Calculate sensitivity and specificity for each cutoff and severity grade
for (cutoff in cutoffs) {
  predictions <- ifelse(scores_w_pheno_baseline_viral$myeloid_z_score >= cutoff, 1, 0)
  
  # Calculate sensitivity
  true_positives <- sum(predictions[scores_w_pheno_baseline_viral$outcome_mild == 1] == 1)
  total_positives <- sum(scores_w_pheno_baseline_viral$outcome_mild == 1)
  sensitivity_mild <- ifelse(total_positives > 0, true_positives / total_positives, NA)  # Use NA if no positives
  
  # Calculate specificity
  true_negatives <- sum(predictions[scores_w_pheno_baseline_viral$outcome_mild == 0] == 0)
  total_negatives <- sum(scores_w_pheno_baseline_viral$outcome_mild == 0)
  specificity_mild <- ifelse(total_negatives > 0, true_negatives / total_negatives, NA)  # Use NA if no negatives
  
  # Store results
  results_mild_viral <- rbind(results_mild_viral, data.frame(cutoff = cutoff, sensitivity_mild_viral = sensitivity_mild, specificity_mild_viral = specificity_mild))
}


# Initialize a data frame to store results
results_severe_viral <- data.frame()

# Calculate sensitivity and specificity for each cutoff and severity grade
for (cutoff in cutoffs) {
  predictions <- ifelse(scores_w_pheno_baseline_viral$myeloid_z_score >= cutoff, 1, 0)
  
  # Calculate sensitivity
  true_positives <- sum(predictions[scores_w_pheno_baseline_viral$outcome_severe == 1 ] == 1)
  total_positives <- sum(scores_w_pheno_baseline_viral$outcome_severe ==1)
  sensitivity_severe <- ifelse(total_positives > 0, true_positives / total_positives, NA)  # Use NA if no positives
  
  # Calculate specificity
  true_negatives <- sum(predictions[scores_w_pheno_baseline_viral$outcome_severe == 0] == 0)
  total_negatives <- sum(scores_w_pheno_baseline_viral$outcome_severe == 0)
  specificity_severe <- ifelse(total_negatives > 0, true_negatives / total_negatives, NA)  # Use NA if no negatives
  
  # Store results
  results_severe_viral <- rbind(results_severe_viral, data.frame(cutoff = cutoff, sensitivity_severe_viral = sensitivity_severe, specificity_severe_viral = specificity_severe))
}


# Initialize a data frame to store results
results_fatal_viral <- data.frame()

# Calculate sensitivity and specificity for each cutoff and severity grade
for (cutoff in cutoffs) {
  predictions <- ifelse(scores_w_pheno_baseline_viral$myeloid_z_score >= cutoff, 1, 0)
  
  # Calculate sensitivity
  true_positives <- sum(predictions[scores_w_pheno_baseline_viral$outcome_fatal == 1 ] == 1)
  total_positives <- sum(scores_w_pheno_baseline_viral$outcome_fatal ==1)
  sensitivity_fatal <- ifelse(total_positives > 0, true_positives / total_positives, NA)  # Use NA if no positives
  
  # Calculate specificity
  true_negatives <- sum(predictions[scores_w_pheno_baseline_viral$outcome_fatal == 0] == 0)
  total_negatives <- sum(scores_w_pheno_baseline_viral$outcome_fatal == 0)
  specificity_fatal <- ifelse(total_negatives > 0, true_negatives / total_negatives, NA)  # Use NA if no negatives
  
  # Store results
  results_fatal_viral <- rbind(results_fatal_viral, data.frame(cutoff = cutoff, sensitivity_fatal_viral = sensitivity_fatal, specificity_fatal_viral = specificity_fatal))
}


#Myeloid
scores_w_pheno_baseline_bacterial <- scores_w_pheno_baseline %>%
  filter(infection_type == "Bacterial" | condition == "healthy") %>%
  mutate(outcome_mild = ifelse(severity_grades %in% c("Healthy"), 0, 1),
         outcome_severe = ifelse(severity_grades %in% c("Severe", "Fatal"), 1, 0),
         outcome_fatal = ifelse(severity_grades %in% c("Fatal"), 1, 0))


# Initialize a data frame to store results
results_mild_bacterial <- data.frame()

# Calculate sensitivity and specificity for each cutoff and severity grade
for (cutoff in cutoffs) {
  predictions <- ifelse(scores_w_pheno_baseline_bacterial$myeloid_z_score >= cutoff, 1, 0)
  
  # Calculate sensitivity
  true_positives <- sum(predictions[scores_w_pheno_baseline_bacterial$outcome_mild == 1] == 1)
  total_positives <- sum(scores_w_pheno_baseline_bacterial$outcome_mild == 1)
  sensitivity_mild <- ifelse(total_positives > 0, true_positives / total_positives, NA)  # Use NA if no positives
  
  # Calculate specificity
  true_negatives <- sum(predictions[scores_w_pheno_baseline_bacterial$outcome_mild == 0] == 0)
  total_negatives <- sum(scores_w_pheno_baseline_bacterial$outcome_mild == 0)
  specificity_mild <- ifelse(total_negatives > 0, true_negatives / total_negatives, NA)  # Use NA if no negatives
  
  # Store results
  results_mild_bacterial <- rbind(results_mild_bacterial, data.frame(cutoff = cutoff, sensitivity_mild_bacterial = sensitivity_mild, specificity_mild_bacterial = specificity_mild))
}


# Initialize a data frame to store results
results_severe_bacterial <- data.frame()

# Calculate sensitivity and specificity for each cutoff and severity grade
for (cutoff in cutoffs) {
  predictions <- ifelse(scores_w_pheno_baseline_bacterial$myeloid_z_score >= cutoff, 1, 0)
  
  # Calculate sensitivity
  true_positives <- sum(predictions[scores_w_pheno_baseline_bacterial$outcome_severe == 1 ] == 1)
  total_positives <- sum(scores_w_pheno_baseline_bacterial$outcome_severe ==1)
  sensitivity_severe <- ifelse(total_positives > 0, true_positives / total_positives, NA)  # Use NA if no positives
  
  # Calculate specificity
  true_negatives <- sum(predictions[scores_w_pheno_baseline_bacterial$outcome_severe == 0] == 0)
  total_negatives <- sum(scores_w_pheno_baseline_bacterial$outcome_severe == 0)
  specificity_severe <- ifelse(total_negatives > 0, true_negatives / total_negatives, NA)  # Use NA if no negatives
  
  # Store results
  results_severe_bacterial <- rbind(results_severe_bacterial, data.frame(cutoff = cutoff, sensitivity_severe_bacterial = sensitivity_severe, specificity_severe_bacterial = specificity_severe))
}


# Initialize a data frame to store results
results_fatal_bacterial <- data.frame()

# Calculate sensitivity and specificity for each cutoff and severity grade
for (cutoff in cutoffs) {
  predictions <- ifelse(scores_w_pheno_baseline_bacterial$myeloid_z_score >= cutoff, 1, 0)
  
  # Calculate sensitivity
  true_positives <- sum(predictions[scores_w_pheno_baseline_bacterial$outcome_fatal == 1 ] == 1)
  total_positives <- sum(scores_w_pheno_baseline_bacterial$outcome_fatal ==1)
  sensitivity_fatal <- ifelse(total_positives > 0, true_positives / total_positives, NA)  # Use NA if no positives
  
  # Calculate specificity
  true_negatives <- sum(predictions[scores_w_pheno_baseline_bacterial$outcome_fatal == 0] == 0)
  total_negatives <- sum(scores_w_pheno_baseline_bacterial$outcome_fatal == 0)
  specificity_fatal <- ifelse(total_negatives > 0, true_negatives / total_negatives, NA)  # Use NA if no negatives
  
  # Store results
  results_fatal_bacterial <- rbind(results_fatal_bacterial, data.frame(cutoff = cutoff, sensitivity_fatal_bacterial = sensitivity_fatal, specificity_fatal_bacterial = specificity_fatal))
}


sensitivity_specificity_table <- left_join(results_mild_viral, results_severe_viral, by = c("cutoff")) %>%
  left_join(.,results_fatal_viral, by = c("cutoff"))%>%
  left_join(.,results_mild_bacterial, by = c("cutoff")) %>%
  left_join(.,results_severe_bacterial, by = c("cutoff")) %>%
  left_join(.,results_fatal_bacterial, by = c("cutoff"))


myeloid_sensitivity_specificity_plot <- ggplot(sensitivity_specificity_table, aes(x = cutoff)) +
  geom_line(aes(y = sensitivity_severe_viral, color = "Viral")) +
  geom_line(aes(y = specificity_severe_viral, color = "Viral"), linetype = "dashed") +
  geom_line(aes(y = sensitivity_severe_bacterial, color = "Bacterial")) +
  geom_line(aes(y = specificity_severe_bacterial, color = "Bacterial"), linetype= "dashed") +
  labs(x = "Myeloid Dysregulation Score", y = "Sensitivity/Specificity", color = "Outcome") +
  scale_color_manual(values = c("Viral" = "mediumpurple2", "Bacterial" = "goldenrod2")) +
  theme_minimal() +
  # xlim(-2.5, 6) +
  theme(legend.position = "none")

myeloid_sensitivity_specificity_plot




#Lymphoid
scores_w_pheno_baseline_viral <- scores_w_pheno_baseline %>%
  filter(infection_type == "Viral" | condition == "healthy") %>%
  mutate(outcome_mild = ifelse(severity_grades %in% c("Healthy"), 0, 1),
         outcome_severe = ifelse(severity_grades %in% c("Severe", "Fatal"), 1, 0),
         outcome_fatal = ifelse(severity_grades %in% c("Fatal"), 1, 0))

# Create a sequence of cut-off values
cutoffs <- seq(min(scores_w_pheno_baseline_viral$lymphoid_z_score, na.rm = TRUE), 
               max(scores_w_pheno_baseline_viral$lymphoid_z_score, na.rm = TRUE), 
               length.out = 100)

# Initialize a data frame to store results
results_mild_viral <- data.frame()

# Calculate sensitivity and specificity for each cutoff and severity grade
for (cutoff in cutoffs) {
  predictions <- ifelse(scores_w_pheno_baseline_viral$lymphoid_z_score >= cutoff, 1, 0)
  
  # Calculate sensitivity
  true_positives <- sum(predictions[scores_w_pheno_baseline_viral$outcome_mild == 1] == 1)
  total_positives <- sum(scores_w_pheno_baseline_viral$outcome_mild == 1)
  sensitivity_mild <- ifelse(total_positives > 0, true_positives / total_positives, NA)  # Use NA if no positives
  
  # Calculate specificity
  true_negatives <- sum(predictions[scores_w_pheno_baseline_viral$outcome_mild == 0] == 0)
  total_negatives <- sum(scores_w_pheno_baseline_viral$outcome_mild == 0)
  specificity_mild <- ifelse(total_negatives > 0, true_negatives / total_negatives, NA)  # Use NA if no negatives
  
  # Store results
  results_mild_viral <- rbind(results_mild_viral, data.frame(cutoff = cutoff, sensitivity_mild_viral = sensitivity_mild, specificity_mild_viral = specificity_mild))
}


# Initialize a data frame to store results
results_severe_viral <- data.frame()

# Calculate sensitivity and specificity for each cutoff and severity grade
for (cutoff in cutoffs) {
  predictions <- ifelse(scores_w_pheno_baseline_viral$lymphoid_z_score >= cutoff, 1, 0)
  
  # Calculate sensitivity
  true_positives <- sum(predictions[scores_w_pheno_baseline_viral$outcome_severe == 1 ] == 1)
  total_positives <- sum(scores_w_pheno_baseline_viral$outcome_severe ==1)
  sensitivity_severe <- ifelse(total_positives > 0, true_positives / total_positives, NA)  # Use NA if no positives
  
  # Calculate specificity
  true_negatives <- sum(predictions[scores_w_pheno_baseline_viral$outcome_severe == 0] == 0)
  total_negatives <- sum(scores_w_pheno_baseline_viral$outcome_severe == 0)
  specificity_severe <- ifelse(total_negatives > 0, true_negatives / total_negatives, NA)  # Use NA if no negatives
  
  # Store results
  results_severe_viral <- rbind(results_severe_viral, data.frame(cutoff = cutoff, sensitivity_severe_viral = sensitivity_severe, specificity_severe_viral = specificity_severe))
}


# Initialize a data frame to store results
results_fatal_viral <- data.frame()

# Calculate sensitivity and specificity for each cutoff and severity grade
for (cutoff in cutoffs) {
  predictions <- ifelse(scores_w_pheno_baseline_viral$lymphoid_z_score >= cutoff, 1, 0)
  
  # Calculate sensitivity
  true_positives <- sum(predictions[scores_w_pheno_baseline_viral$outcome_fatal == 1 ] == 1)
  total_positives <- sum(scores_w_pheno_baseline_viral$outcome_fatal ==1)
  sensitivity_fatal <- ifelse(total_positives > 0, true_positives / total_positives, NA)  # Use NA if no positives
  
  # Calculate specificity
  true_negatives <- sum(predictions[scores_w_pheno_baseline_viral$outcome_fatal == 0] == 0)
  total_negatives <- sum(scores_w_pheno_baseline_viral$outcome_fatal == 0)
  specificity_fatal <- ifelse(total_negatives > 0, true_negatives / total_negatives, NA)  # Use NA if no negatives
  
  # Store results
  results_fatal_viral <- rbind(results_fatal_viral, data.frame(cutoff = cutoff, sensitivity_fatal_viral = sensitivity_fatal, specificity_fatal_viral = specificity_fatal))
}


#lymphoid
scores_w_pheno_baseline_bacterial <- scores_w_pheno_baseline %>%
  filter(infection_type == "Bacterial" | condition == "healthy") %>%
  mutate(outcome_mild = ifelse(severity_grades %in% c("Healthy"), 0, 1),
         outcome_severe = ifelse(severity_grades %in% c("Severe", "Fatal"), 1, 0),
         outcome_fatal = ifelse(severity_grades %in% c("Fatal"), 1, 0))


# Initialize a data frame to store results
results_mild_bacterial <- data.frame()

# Calculate sensitivity and specificity for each cutoff and severity grade
for (cutoff in cutoffs) {
  predictions <- ifelse(scores_w_pheno_baseline_bacterial$lymphoid_z_score >= cutoff, 1, 0)
  
  # Calculate sensitivity
  true_positives <- sum(predictions[scores_w_pheno_baseline_bacterial$outcome_mild == 1] == 1)
  total_positives <- sum(scores_w_pheno_baseline_bacterial$outcome_mild == 1)
  sensitivity_mild <- ifelse(total_positives > 0, true_positives / total_positives, NA)  # Use NA if no positives
  
  # Calculate specificity
  true_negatives <- sum(predictions[scores_w_pheno_baseline_bacterial$outcome_mild == 0] == 0)
  total_negatives <- sum(scores_w_pheno_baseline_bacterial$outcome_mild == 0)
  specificity_mild <- ifelse(total_negatives > 0, true_negatives / total_negatives, NA)  # Use NA if no negatives
  
  # Store results
  results_mild_bacterial <- rbind(results_mild_bacterial, data.frame(cutoff = cutoff, sensitivity_mild_bacterial = sensitivity_mild, specificity_mild_bacterial = specificity_mild))
}


# Initialize a data frame to store results
results_severe_bacterial <- data.frame()

# Calculate sensitivity and specificity for each cutoff and severity grade
for (cutoff in cutoffs) {
  predictions <- ifelse(scores_w_pheno_baseline_bacterial$lymphoid_z_score >= cutoff, 1, 0)
  
  # Calculate sensitivity
  true_positives <- sum(predictions[scores_w_pheno_baseline_bacterial$outcome_severe == 1 ] == 1)
  total_positives <- sum(scores_w_pheno_baseline_bacterial$outcome_severe ==1)
  sensitivity_severe <- ifelse(total_positives > 0, true_positives / total_positives, NA)  # Use NA if no positives
  
  # Calculate specificity
  true_negatives <- sum(predictions[scores_w_pheno_baseline_bacterial$outcome_severe == 0] == 0)
  total_negatives <- sum(scores_w_pheno_baseline_bacterial$outcome_severe == 0)
  specificity_severe <- ifelse(total_negatives > 0, true_negatives / total_negatives, NA)  # Use NA if no negatives
  
  # Store results
  results_severe_bacterial <- rbind(results_severe_bacterial, data.frame(cutoff = cutoff, sensitivity_severe_bacterial = sensitivity_severe, specificity_severe_bacterial = specificity_severe))
}


# Initialize a data frame to store results
results_fatal_bacterial <- data.frame()

# Calculate sensitivity and specificity for each cutoff and severity grade
for (cutoff in cutoffs) {
  predictions <- ifelse(scores_w_pheno_baseline_bacterial$lymphoid_z_score >= cutoff, 1, 0)
  
  # Calculate sensitivity
  true_positives <- sum(predictions[scores_w_pheno_baseline_bacterial$outcome_fatal == 1 ] == 1)
  total_positives <- sum(scores_w_pheno_baseline_bacterial$outcome_fatal ==1)
  sensitivity_fatal <- ifelse(total_positives > 0, true_positives / total_positives, NA)  # Use NA if no positives
  
  # Calculate specificity
  true_negatives <- sum(predictions[scores_w_pheno_baseline_bacterial$outcome_fatal == 0] == 0)
  total_negatives <- sum(scores_w_pheno_baseline_bacterial$outcome_fatal == 0)
  specificity_fatal <- ifelse(total_negatives > 0, true_negatives / total_negatives, NA)  # Use NA if no negatives
  
  # Store results
  results_fatal_bacterial <- rbind(results_fatal_bacterial, data.frame(cutoff = cutoff, sensitivity_fatal_bacterial = sensitivity_fatal, specificity_fatal_bacterial = specificity_fatal))
}


sensitivity_specificity_table <- left_join(results_mild_viral, results_severe_viral, by = c("cutoff")) %>%
  left_join(.,results_fatal_viral, by = c("cutoff"))%>%
  left_join(.,results_mild_bacterial, by = c("cutoff")) %>%
  left_join(.,results_severe_bacterial, by = c("cutoff")) %>%
  left_join(.,results_fatal_bacterial, by = c("cutoff"))





lymphoid_sensitivity_specificity_plot <- ggplot(sensitivity_specificity_table, aes(x = cutoff)) +
  geom_line(aes(y = sensitivity_severe_viral, color = "Viral")) +
  geom_line(aes(y = specificity_severe_viral, color = "Viral"), linetype = "dashed") +
  geom_line(aes(y = sensitivity_severe_bacterial, color = "Bacterial")) +
  geom_line(aes(y = specificity_severe_bacterial, color = "Bacterial"), linetype= "dashed") +
  labs(x = "Lymphoid Dysregulation Score", y = "Sensitivity/Specificity", color = "Outcome") +
  scale_color_manual(values = c("Viral" = "mediumpurple2", "Bacterial" = "goldenrod2")) +
  theme_minimal() +
  # xlim(-2.5, 6) +
  theme(legend.position = "none")

lymphoid_sensitivity_specificity_plot

row1 <- plot_grid("",subspace_threshold_dotplot_infection, "", labels = c("","A",""), rel_widths = c(0.25,1,0.25), nrow = 1)

row2 <- plot_grid(myeloid_sensitivity_specificity_plot, lymphoid_sensitivity_specificity_plot, labels = c("B","C"), nrow = 1)

final_plot <- plot_grid(row1, row2, ncol = 1)

print(final_plot)



#############Threshold by predominant enrollment location###########

subspace_summary_stats_groups <- scores_w_pheno %>%
  filter(timepoint == "baseline",  !is.na(site), site != "ufl", site != "trinity", site != "imx", condition != "healthy", site != "cchmc") %>%
  mutate(group = case_when(
    site == "savemore" ~ "Non-ICU",
    site == "stanford" ~ "ICU",
    site == "amsterdam" & grepl("panamo", patient_id) ~ "ICU",
    site == "amsterdam" & grepl("optimact", patient_id) ~ "ED",
    site == "amsterdam" & grepl("EB", patient_id) ~ "Non-ICU",
    site == "trinity" ~ "ICU",
    site == "charles" ~ "ICU",
    site == "victas" ~ "ICU",
    site == "acutelines" ~ "ED"
  ))%>%
  group_by(group) %>%
  summarise(
    myeloid_mean = mean(myeloid_z_score, na.rm = TRUE),
    myeloid_sd = sd(myeloid_z_score, na.rm = TRUE),
    lymphoid_mean = mean(lymphoid_z_score, na.rm = TRUE),
    lymphoid_sd = sd(lymphoid_z_score, na.rm = TRUE)
  )

subspace_groups <- scores_w_pheno %>%
  filter(timepoint == "baseline",  !is.na(site), site != "ufl", site != "trinity", site != "imx", condition != "healthy", site != "cchmc") %>%
  mutate(group = case_when(
    site == "savemore" ~ "Non-ICU",
    site == "stanford" ~ "ICU",
    site == "amsterdam" & grepl("panamo", patient_id) ~ "ICU",
    site == "amsterdam" & grepl("optimact", patient_id) ~ "ED",
    site == "amsterdam" & grepl("EB", patient_id) ~ "Non-ICU",
    site == "trinity" ~ "ICU",
    site == "charles" ~ "ICU",
    site == "victas" ~ "ICU",
    site == "acutelines" ~ "ED"
  ))

subspace_threshold_dotplot_group <- ggplot(subspace_groups, aes(x = myeloid_z_score, y = lymphoid_z_score, color = group))+
  geom_ellipse(inherit.aes = F,data = filter(subspace_summary_stats_groups, group == "Non-ICU"), aes(x0 = myeloid_mean, y0 = lymphoid_mean, a = myeloid_sd, b = lymphoid_sd, angle = 0), color = "dodgerblue", fill = "dodgerblue", alpha = 0.2) +
  geom_ellipse(inherit.aes = F,data = filter(subspace_summary_stats_groups, group == "ED"), aes(x0 = myeloid_mean, y0 = lymphoid_mean, a = myeloid_sd, b = lymphoid_sd, angle = 0), color = "goldenrod2",fill = "goldenrod2", alpha = 0.2) +
  geom_ellipse(inherit.aes = F,data = filter(subspace_summary_stats_groups,  group == "ICU"), aes(x0 = myeloid_mean, y0 = lymphoid_mean, a = myeloid_sd, b = lymphoid_sd, angle = 0), color = "palevioletred",fill = "palevioletred", alpha = 0.2) +
  geom_point(size = 0.15)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  theme_minimal() +
  scale_color_manual(values = c("Non-ICU" = "dodgerblue", "ED" = "goldenrod2", "ICU" = "palevioletred"))+
  theme(legend.position = "right",axis.line = element_line(size = 1, arrow = arrow(type = "closed", length = unit(0.25, "inches"))), panel.border = element_blank(), panel.background = element_blank())+ 
  ylab("Lymphoid Dysregulation Score")+
  xlab("Myeloid Dysregulation Score")


subspace_threshold_dotplot_group
