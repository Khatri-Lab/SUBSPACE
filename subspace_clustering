library(tidyverse)
library(ComplexHeatmap)
library(cluster)
library(factoextra)
library(dendextend)
library(scales)
library(pvclust)
library(igraph)



set.seed(777)

###Pulls the data and sets it up for clustering
subspace_score_table <- read.csv("~/subspace_score_table.csv")%>%
  filter(!is.na(Etiology))

clustering_scores = c("mod1_score","mod2_score","mod3_score","mod4_score","inflammopathic_score","adaptive_score","coagulopathic_score","yao_IA_score","yao_IC_score","yao_IN_score","davenport_SRSq","cano_SRSq", 
                      "wong_score",
                      "mars1_score","mars2_score","mars3_score","mars4_score")

score_table_clustering <- subspace_score_table%>%
  dplyr::select(all_of(clustering_scores))

score_table_matrix <- as.matrix(score_table_clustering)

score_table_matrix <- scale(score_table_matrix)

##clustering
scores <- data.frame(t(score_table_matrix))

d <- dist(scores, method = "euclidean")

hc <- hclust(d, method = "ward.D2")
plot(hc, cex = 0.6, hang = -1)

sub_grp <- cutree(hc, k = 4)

scores %>%
  mutate(cluster = sub_grp) %>%
  head

hc_plot <- plot(hc, cex = 0.6)
rect.hclust(hc, k = 4, border = 2:4)

#Creates the heatmap
column_order = c("coagulopathic_score", "yao_IC_score","mars1_score", "mod2_score", "inflammopathic_score", "yao_IN_score", "mod1_score", "mars2_score","wong_score", "mod3_score","mars4_score", "mars3_score", "adaptive_score", "mod4_score", "yao_IA_score")

#reorganizes the columns in the matrix
score_table_matrix <- score_table_matrix[, column_order]

#Makes the unsupervised heatmap
public_heatmap <- Heatmap(score_table_matrix,   show_row_dend = FALSE, column_dend_reorder = F, clustering_method_columns = "ward.D2",column_split = 4, top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#A781BA","#F3766E","#1CBDC2","#7CAF41")))),column_title = NULL, column_labels = c("Coagulopathic", "Yao Coagulopathic","MARS 1", "SoM Mod 2", "Inflammopathic", "Yao Innate", "SoM Mod 1", "MARS 2","Wong Score", "SoM Mod 3","MARS 4", "MARS 3", "Adaptive", "SoM Mod 4", "Yao Adaptive")
)

public_heatmap

pca <- fviz_cluster(list(data = scores, cluster = sub_grp))

rownames(pca$data) <- c("SoM Mod 1", "SoM Mod 2", "SoM Mod 3", "SoM Mod 4", "Inflammopathic", "Adaptive", "Coagulopathic", "Yao Adaptive", "Yao Coagulopathic", "Yao Innate", "Davenport SRSq", "Cano-Gamez SRSq", "Wong Score", "MARS 1", "MARS 2", "MARS 3", "MARS 4")


order <- c("SoM Mod 2","Inflammopathic","Davenport SRSq","Yao Innate","SoM Mod 3","SoM Mod 4","Adaptive","Yao Adaptive","Coagulopathic","Yao Coagulopathic","Wong Score","MARS 1","MARS 3","MARS 4")
#reorder pca$data by rowname
#pca$data <- pca$data[c("SoM Mod 1"),]

pca_ggplot <- pca$data %>%
  ggplot(aes(x = x, y = y, color = factor(cluster))) +
  geom_point() +
  scale_color_manual(values = c("1" = "#F3766E", "2" = "#7CAF41", "3" = "#1CBDC2", "4" = "#A781BA")) +  # Adjust colors as needed
  #add labels to each point
  geom_text(aes(label = rownames(pca$data)), hjust = 0, vjust = 0, color = "black") +
  #add circle around the clusters
  geom_polygon(data = pca$data[order,], aes(x = x, y = y, group = cluster, fill = factor(cluster)), alpha = 0.2) +
  scale_fill_manual(values = c("#F3766E", "#7CAF41", "#1CBDC2", "#A781BA")) +
  theme(legend.position = "none")


# Print the ggplot object
print(pca_ggplot)

##########################
#network analysis
network_matrix <- as.matrix(score_table_clustering)


#creates a correlation matrix for the columns in score_table_matrix
correlation_matrix <- cor(network_matrix, method = "spearman")

#generates edges from correlation matrix
network_matrix <- as.matrix(score_table_clustering)


#creates a correlation matrix for the columns in score_table_matrix
correlation_matrix <- cor(network_matrix, method = "spearman")

quantile(abs(correlation_matrix), probs = c(0.25,.5,.75))

#generates edges from correlation matrix
threshold <- 0.366


# Define edge width based on correlation values

t = which(correlation_matrix > threshold & lower.tri(correlation_matrix),arr.ind=TRUE)
t <- cbind(t, correlation_matrix[which(correlation_matrix > threshold & lower.tri(correlation_matrix),arr.ind=TRUE)])

t <- data.frame(t)%>%
  mutate(row = case_when(
    row == 1 ~ rownames(correlation_matrix)[1],
    row == 2 ~ rownames(correlation_matrix)[2],
    row == 3 ~ rownames(correlation_matrix)[3],
    row == 4 ~ rownames(correlation_matrix)[4],
    row == 5 ~ rownames(correlation_matrix)[5],
    row == 6 ~ rownames(correlation_matrix)[6],
    row == 7 ~ rownames(correlation_matrix)[7],
    row == 8 ~ rownames(correlation_matrix)[8],
    row == 9 ~ rownames(correlation_matrix)[9],
    row == 10 ~ rownames(correlation_matrix)[10],
    row == 11 ~ rownames(correlation_matrix)[11],
    row == 12 ~ rownames(correlation_matrix)[12],
    row == 13 ~ rownames(correlation_matrix)[13],
    row == 14 ~ rownames(correlation_matrix)[14],
    row == 15 ~ rownames(correlation_matrix)[15],
    row == 16 ~ rownames(correlation_matrix)[16],
    row == 17 ~ rownames(correlation_matrix)[17]
  ),
  col= case_when(
    col == 1 ~ rownames(correlation_matrix)[1],
    col == 2 ~ rownames(correlation_matrix)[2],
    col == 3 ~ rownames(correlation_matrix)[3],
    col == 4 ~ rownames(correlation_matrix)[4],
    col == 5 ~ rownames(correlation_matrix)[5],
    col == 6 ~ rownames(correlation_matrix)[6],
    col == 7 ~ rownames(correlation_matrix)[7],
    col == 8 ~ rownames(correlation_matrix)[8],
    col == 9 ~ rownames(correlation_matrix)[9],
    col == 10 ~ rownames(correlation_matrix)[10],
    col == 11 ~ rownames(correlation_matrix)[11],
    col == 12 ~ rownames(correlation_matrix)[12],
    col == 13 ~ rownames(correlation_matrix)[13],
    col == 14 ~ rownames(correlation_matrix)[14],
    col == 15 ~ rownames(correlation_matrix)[15],
    col == 16 ~ rownames(correlation_matrix)[16],
    col == 17 ~ rownames(correlation_matrix)[17]
  )
  )


##this adds the correlation to the graph as an edge attribute "V3"
graph_data=graph.data.frame(t,directed=F)
E(graph_data)$width <- case_when(E(graph_data)$V3 <0.4 ~1, 
                                 E(graph_data)$V3 >=0.4 & E(graph_data)$V3 <0.6 ~2,
                                 E(graph_data)$V3 >=0.6 & E(graph_data)$V3 <0.8 ~4,
                                 E(graph_data)$V3 >=0.8 ~5)


# Assign community colors

cluster <- cluster_fast_greedy(graph_data)

members <- membership(cluster)

V(graph_data)$color <- rainbow(max(members))[members]

# Plot the network

plot(graph_data, 
     
     vertex.size=10,
     
     vertex.label.cex=0.7,
     
     edge.width=E(graph_data)$width,
     
     layout=layout_with_fr)
