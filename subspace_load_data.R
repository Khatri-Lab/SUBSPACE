library(tidyverse)
library(data.table)
library(pals)
library(rawr)
library(RColorBrewer)
library(pROC)
library(ggpubr)
library(ggbeeswarm)
library(ggrastr)
library(cowplot)
library(interactions)
library(ggcorrplot)
library(ComplexHeatmap)
library(cluster)
library(factoextra)
library(dendextend)
library(stats)
library(scales)
library(pvclust)
library(clustMixType)
library(FactoMineR)
library(igraph)
library(gtsummary)
library(ggalluvial)
library(survival)
library(ggsurvfit)
library(ggdendro)
library(grid)
library(forestplot)
library(lme4)
library(ggExtra)
library(ggplotify)
library(ggforce)
library(grid)
library(gridExtra)

###############################################################
base_text_size <- 10
base_font_family <- "Helvetica"
text_size_labels <- (base_text_size - 1) / 2.835# in mm, ~ 9 pt
line_size_main <- 0.75
line_size_main_mm <- 0.75 / 2.835# in mm, ~0.75 in pt
line_size_text_repel <- 0.5
line_size_text_repel_mm <- 0.5 / 2.835# in mm, ~0.5 in pt
point_size <- 0.25
cex_val <- 1# For the beeswarm spread
gap_size <- unit(4, "pt")
black_text_colour <- "#252525"
white_text_colour <- "#EFEFEF"
black_line_colour <- "#505050"
grey_line_colour <- "#CCCCCC"
###############################################################


###############################################################
text_setting_sub_small_gg <- element_text(size = base_text_size - 1, family = base_font_family, colour = black_text_colour)
text_setting_small_gg <- element_text(size = base_text_size, family = base_font_family, colour = black_text_colour)
text_setting_small_title_gg <- element_text(size = base_text_size, family = base_font_family, colour = black_text_colour, hjust = 0.5, vjust = 0)
basic_gg_theme <- function()
{
  theme_bw(base_size = base_text_size, base_family = base_font_family, base_line_size = line_size_main_mm, base_rect_size = line_size_main_mm) +
    theme(text = text_setting_small_gg, axis.text = text_setting_sub_small_gg, axis.title = text_setting_small_gg,
          legend.text = text_setting_sub_small_gg, legend.title = text_setting_small_gg, legend.title.align = 0,
          plot.title = text_setting_small_title_gg,
          strip.text = text_setting_small_gg,
          line = element_line(linewidth = line_size_text_repel_mm, colour = black_line_colour),
          panel.border = element_rect(linewidth = line_size_main_mm, colour = black_line_colour, fill = "transparent"),
          axis.ticks = element_line(linewidth = line_size_text_repel_mm, colour = black_line_colour),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), plot.background = element_blank(), strip.background = element_blank(),
          legend.background = element_blank(), legend.box.background = element_blank(), legend.key = element_blank(),
          legend.key.height = unit(base_text_size, "pt"), legend.key.width = unit(base_text_size - 1, "pt"), legend.margin = ggplot2::margin(1, 1, 1, 1, "pt"), legend.spacing.y = unit(1.5, "pt"),
          plot.margin = ggplot2::margin(3, 3, 3, 3, unit = "pt")
    )
}
theme_set(basic_gg_theme())
###############################################################


###############################################################
text_setting_sub_small <- gpar(fontsize = base_text_size - 1, fontfamily = base_font_family, col = black_text_colour)
text_setting_sub_small_white <- gpar(fontsize = base_text_size - 1, fontfamily = base_font_family, col = white_text_colour)
text_setting_small <- gpar(fontsize = base_text_size, fontfamily = base_font_family, col = black_text_colour)
text_setting_small_white <- gpar(fontsize = base_text_size, fontfamily = base_font_family, col = white_text_colour)
dend_lines_gp <- gpar(col = black_line_colour, lwd = line_size_text_repel, lty = 1)
dotted_lines_gp <- gpar(col = black_line_colour, lwd = line_size_text_repel, lty = 2)
ht_opt("heatmap_row_names_gp" = text_setting_small)
ht_opt("heatmap_column_names_gp" = text_setting_small)
ht_opt("heatmap_row_title_gp" = text_setting_small)
ht_opt("heatmap_column_title_gp" = text_setting_small)
ht_opt("TITLE_PADDING" = gap_size)
ht_opt("ROW_ANNO_PADDING" = gap_size)
ht_opt("COLUMN_ANNO_PADDING" = gap_size)
ht_opt("show_parent_dend_line" = FALSE)
ht_opt("legend_border" = gpar(fill = "transparent"))
###############################################################

# Prints extra zeros for sigfigs
trail.sigfig <- function(x, ndig = 2)
{
  loc.dec <- str_locate(as.character(x), "\\.")[1]
  str.len <- str_length(as.character(x))
  ndec <- str.len - loc.dec
  rounded <- round(x, digits = ndig)
  nzero <- (loc.dec + ndig) - str_length(rounded)
  output <- if_else(str_length(rounded) < loc.dec + ndig, paste0(rounded, rep("0", nzero)), as.character(rounded))
  return(output)
}

public_score_table <- read.csv("https://raw.githubusercontent.com/Khatri-Lab/SUBSPACE/main/public_score_table.csv")%>%
  filter(!is.na(Etiology))

public_score_table$subgroup = factor(public_score_table$subgroup, levels = c("balanced","lymphoid dysregulation", "myeloid dysregulation","system-wide"))

public_score_table$severity_grades <- factor(public_score_table$severity_grades, levels = c("Healthy", "Non-Severe", "Severe", "Fatal"))


scores_w_pheno <- read.csv("https://raw.githubusercontent.com/Khatri-Lab/SUBSPACE/refs/heads/main/subspace_score_table.csv")

scores_w_pheno$subgroup = factor(scores_w_pheno$subgroup, levels = c("balanced","lymphoid dysregulation", "myeloid dysregulation","system-wide"))
scores_w_pheno$severity_grades <- factor(scores_w_pheno$severity_grades, levels = c("Healthy", "Non-Severe", "Severe", "Fatal"))

acutelines <- scores_w_pheno%>%
  filter(site == "acutelines")

stanford <- scores_w_pheno%>%
  filter(site == "stanford")

amsterdam <- scores_w_pheno%>%
  filter(site == "amsterdam")

cchmc <- scores_w_pheno%>%
  filter(site == "cchmc")

charles <- scores_w_pheno%>%
  filter(site == "charles")

victas <- scores_w_pheno%>%
  filter(site == "victas")

savemore <- scores_w_pheno%>%
  filter(site == "savemore")

ufl <- scores_w_pheno%>%
  filter(site == "ufl")

glue <- read.csv("https://raw.githubusercontent.com/Khatri-Lab/SUBSPACE/refs/heads/main/glue_score_table.csv")

messi <- read.csv("https://raw.githubusercontent.com/Khatri-Lab/SUBSPACE/refs/heads/main/messi_score_table.csv")%>%
  mutate(subgroup = tolower(subgroup))

messi$subgroup = factor(messi$subgroup, levels = c("balanced","lymphoid dysregulation", "myeloid dysregulation","system-wide"))

vanish <- read.csv("https://raw.githubusercontent.com/Khatri-Lab/SUBSPACE/refs/heads/main/vanish_score_table.csv")%>%
  mutate(high_myeloid = ifelse(myeloid_score >= quantile(myeloid_score, probs = 0.5), "high", "low"),
         high_lymphoid = ifelse(lymphoid_score >= quantile(lymphoid_score, probs = 0.5), "high", "low"),
         subgroup = ifelse(high_myeloid == "high" & high_lymphoid == "high", "system-wide",
                           ifelse(high_myeloid == "high", "myeloid dysregulation",
                                  ifelse(high_lymphoid == "high", "lymphoid dysregulation", "balanced"))),
         d28_death_yn = ifelse(Characteristics.outcome.day.28. == "Dead", 1, ifelse(Characteristics.outcome.day.28. == "Alive", 0,NA)),
         treatment = ifelse(Characteristics.drug2.per.protocol. == "Hydrocortisone", "hydrocortisone", "no hydrocortisone"),
         lymphoid_score = scale(lymphoid_score))

vanish$subgroup = factor(vanish$subgroup, levels = c("immunocompetent","lymphoid dysregulation", "myeloid dysregulation","hyperinflammatory"))
vanish$high_lymphoid = factor(vanish$high_lymphoid, levels = c("low","high"))
vanish$high_myeloid = factor(vanish$high_myeloid, levels = c("low","high"))
