#Required Packages
library(tidyverse)
library(data.table)
library(pals)
#devtools::install_github('raredd/rawr', force = TRUE)
library(rawr)
library(pROC)
library(ggpubr)
library(scales)
library(ggbeeswarm)
library(ggrastr)
library(ComplexHeatmap)
library(cluster)
library(factoextra)
library(dendextend)
library(stats)
library(scales)
library(pvclust)
library(readr)

sc_data <- read_csv("https://zenodo.org/records/16759647/files/subspace_single_cell_data.csv")

srs_genes <- c("DYRK2", "CCNB1IP1", "TDRD9", "ZAP70", "ARL14EP", "MDC1", "ADGRE3")
srs_down_genes <- c("DYRK2","CCNB1IP1","ZAP70","ARL14EP","MDC1","ADGRE3")
srs_up_genes <- c("TDRD9")

#SoM modules
mod1 = c("NQO2","SLPI","ORM1","KLHL2","ANXA3","TXN","AQP9","BCL6","DOK3","PFKFB4","TYK2")
mod2 = c("BCL2L11","BCAT1","BTBD7","CEP55","HMMR","PRC1","KIF15","CAMP","CEACAM8","DEFA4","LCN2","CTSG","AZU1")
mod3 = c("MAFB","OASL","UBE2L6","VAMP5","CCL2","NAPA","ATG3","VRK2","TMEM123","CASP7")
mod4 = c("DOK2","HLA-DPB1","BUB3","SMYD2","SIDT1","EXOC2","TRIB2","KLRB1")

#endotyping
inflammopathic.up= c("ARG1","LCN2","LTF","OLFM4")
inflammopathic.down = c("HLA-DMB")
adaptive.up = c("YKT6", "PDE4B", "TWISTNB", "BTN2A2", "ZBTB33", "PSMB9", "CAMK4", "TMEM19", "SLC12A7", 
                "TP53BP1", "PLEKHO1", "SLC25A22", "FRS2")                          
adaptive.down = c("GADD45A", "CD24", "S100A12", "STX1A")
coagulopathic.up = c("KCNMB4", "CRISP2", "HTRA1","PPL")
coagulopathic.down = c("RHBDF2", "ZCCHC4", "YKT6", "DDX6", "SENP5", "RAPGEF1","DTX2","RELB")

#Identifies the MARS genes
mars1_up = "BPGM"
mars1_down = "TAP2"
mars2_up = "GADD45A"
mars2_down = "PCGF5"
mars3_up = "AHNAK"
mars3_down = "PDCD10"
mars4_up = "IFIT5"
mars4_down = "NOP53"


#Yao Score
#Yao Score genes
yao_IA_up <- c("ZNF831", "CD3G", "MME", "BTN3A2"#, "HLA-DPA1"
)
yao_IA_down <- "STOM"
yao_IN_up <- c("HK3","SERPINB1")
yao_IN_down <- c("EPB42", "GSPT1"#, "LAT"
)
yao_IC_up <- c("SLC1A5", "IGF2BP2", "ANXA3")
yao_IC_down <- "GBP2"

#Wong Score
#is a geometric mean of all 100 genes
wong_genes = c("APAF1", "ARPC5","ASAH1","ATP2B2","BCL6","BMPR2","BTK","CAMK2D","CAMK2G","CAMK4","CASP1","CASP2","CASP4","CASP8","CD247","CD3E","CD3G","CD79A",
               "CREB1","CREB5","CSNK1A1","CTNNB1","DAPP1","DBT","EP300","FAS","FCGR2A","FCGR2C","FYN","GK","GNAI3","HDAC4","HLA-DMA","HLA-DOA","ICAM3","IL1A",
               "INPP5D","ITGAM","ITGAV","ITGAX","JAK1","JAK2","KAT2B","LAT2","LYN","MAP2K4","MAP3K1","MAP3K3","MAP3K5","MAP3K7","MAP4K1","MAPK14","MDH1","MKNK1",
               "NCOA2","NCR3","NFATC1","PAK2","PDPR","PIAS1","PIK3C2A","PIK3C3","PIK3CA","PIK3CD","PIK3R1","PLCG1","POU2F2","PPP1R12A","PPP2R2A","PPP2R5C",
               "PRKAR1A","PRKCB","PSMB7","PTEN","PTPRC","RAF1","RHOT1","ROCK1","SEMA4F","SEMA6B","SMAD4","SOS1","SOS2","SP1","TAF11","TBK1","TGFBR1","TLE4",
               "TLR1","TLR2","TLR8","TNFSF10","TRA@","TYROBP","UBE3A","USP48","ZAP70","ZDHHC17")

myeloid_detrimental_genes <- c('ANXA3', 'ARG1', 'AZU1', 'CAMP', 'CEACAM8', 'CEP55', 'CRISP2', 'CTSG', 'DEFA4', 'GADD45A', 'HMMR', 'KIF15', 'LCN2', 'LTF', 'OLFM4', 'ORM1', 'PRC1', 'SLPI', 'STOM', 'AQP9', 'BCL6', 'KLHL2', 'PPL', 'HTRA1', 'TYK2', 'SLC1A5', 'STX1A')


neutrophil_protective_genes <- c('ZDHHC17', 'FAS', 'GK', 'ICAM3', 'MME', 'PDE4B', 'PIK3CD', 'PTEN', 'RAF1', 'TLR1', 'PPP1R12A', 'MAPK14', 'SOS2', 'TXN')

monocyte_protective_genes <- c('ASAH1', 'ATG3', 'BCAT1', 'BCL2L11', 'BTK', 'BTN2A2', 'CASP1', 'CCL2', 'CREB1', 'EP300', 'GNAI3', 'IL1A', 'JAK2', 'MAFB', 'MAP3K1', 'MAP3K3', 'PAK2', 'PLEKHO1', 'POU2F2', 'PRKAR1A', 'PRKCB', 'RHBDF2', 'SEMA6B', 'SP1', 'TLE4', 'BMPR2', 'CTNNB1', 'INPP5D', 'ITGAV', 'SLC12A7', 'TBK1', 'VAMP5', 'VRK2', 'YKT6')

lymphoid_protective_genes <- c('ARL14EP', 'BPGM', 'BTN3A2', 'BUB3', 'CAMK4', 'CASP8', 'CCNB1IP1', 'CD247', 'CD3E', 'CD3G', 'DBT', 'DDX6', 'DYRK2', 'JAK1', 'KLRB1', 'MAP4K1', 'NCR3', 'PIK3R1', 'PLCG1', 'PPP2R5C', 'SEMA4F', 'SIDT1', 'SMAD4', 'SMYD2', 'TP53BP1', 'TRIB2', 'ZAP70', 'ZCCHC4', 'ZNF831')

genes_of_interest = c(srs_genes, mod1, mod2, mod3, mod4, inflammopathic.up, inflammopathic.down, adaptive.up, adaptive.down, coagulopathic.up, coagulopathic.down, mars1_up, mars1_down, mars2_up, mars2_down, mars3_up, mars3_down, mars4_up, mars4_down, yao_IA_up, yao_IA_down, yao_IN_up, yao_IN_down, yao_IC_up, yao_IC_down, wong_genes)

sc_scores <- sc_data%>%
  group_by(celltype)%>%
  summarise(mars1_score = mean(mars1_score, na.rm = TRUE),
            coagulopathic_score = mean(coagulopathic_score, na.rm = TRUE),
            yao_IC_score = mean(yao_IC_score, na.rm = TRUE),
            mars2_score = mean(mars2_score, na.rm = TRUE),
            mod1_score = mean(mod1_score, na.rm = TRUE),
            mod2_score = mean(mod2_score, na.rm = TRUE),
            inflammopathic_score = mean(inflammopathic_score, na.rm = TRUE),
            yao_IN_score = mean(yao_IN_score, na.rm = TRUE),
            srs_up_score = mean(srs_up_score, na.rm = TRUE),
            mars3_score = mean(mars3_score, na.rm = TRUE),
            yao_IA_score = mean(yao_IA_score, na.rm = TRUE),
            adaptive_score = mean(adaptive_score, na.rm = TRUE),
            mod4_score = mean(mod4_score, na.rm = TRUE),
            srs_down_score = mean(srs_down_score, na.rm = TRUE),
            mars4_score = mean(mars4_score, na.rm = TRUE), 
            wong_score = mean(wong_score, na.rm = TRUE), 
            mod3_score = mean(mod3_score, na.rm = TRUE)
  )

mat.expr = as.matrix(sc_scores[,-1])
rownames(mat.expr) = sc_scores$celltype
mat.expr = scale(mat.expr)
mat.expr = t(mat.expr)

col.order <- c("Immat Neutrophil","Neutrophil","CD14 Monocyte","CD16 Monocyte", "cDC","pDC","B","PB","CD4 T", "CD8 T", "Treg","NK","Prolif T/NK","HSPC"
)

mat.expr = mat.expr[,col.order]

sc_cell_prop <- sc_data%>%
  filter(!is.na(celltype))%>%
  summarize(`Immat Neutrophil` = log(mean(celltype == "Immat Neutrophil", na.rm = T)*1000),
            Neutrophil = log(mean(celltype == "Neutrophil", na.rm = T)*1000),
            `CD14 Monocyte` = log(mean(celltype == "CD14 Monocyte", na.rm = T)*1000),
            `CD16 Monocyte` = log(mean(celltype == "CD16 Monocyte", na.rm = T)*1000),
            Mast= log(mean(celltype == "Mast", na.rm = T)*1000),
            pDC= log(mean(celltype == "pDC", na.rm = T)*1000),
            cDC= log(mean(celltype == "cDC", na.rm = T)*1000),
            B = log(mean(celltype == "B", na.rm = T)*1000),
            PB = log(mean(celltype == "PB", na.rm = T)*1000),
            `CD4 T` = log(mean(celltype == "CD4 T", na.rm = T)*1000),
            `CD8 T` = log(mean(celltype == "CD8 T", na.rm = T)*1000),
            Treg = log(mean(celltype == "Treg", na.rm = T)*1000),
            NK = log(mean(celltype == "NK", na.rm = T)*1000),
            `Prolif T/NK` = log(mean(celltype == "Prolif T/NK", na.rm = T)*1000),
            HSPC = log(mean(celltype == "HSPC", na.rm = T)*1000)
  )

sc_matrix_int1= t(as.matrix(sc_cell_prop))
colnames(sc_matrix_int1) = "mars1_score"

sc_matrix_int2 <- data.frame(sc_matrix_int1)%>%
  mutate(mars1_score = ifelse(mars1_score>4.5,4.5,
                              ifelse(mars1_score <1.5,1.5,mars1_score)))%>%
  mutate(coagulopathic_score = mars1_score,
         yao_IC_score = mars1_score,
         mars2_score = mars1_score,
         mod1_score = mars1_score,
         mod2_score = mars1_score,
         inflammopathic_score = mars1_score,
         yao_IN_score = mars1_score,
         srs_up_score = mars1_score,
         mars3_score = mars1_score,
         yao_IA_score = mars1_score,
         adaptive_score = mars1_score,
         mod4_score = mars1_score,
         srs_down_score = mars1_score,
         mars4_score = mars1_score,
         wong_score = mars1_score, 
         mod3_score = mars1_score
  )

mat.N = as.matrix(sc_matrix_int2)
mat.N = t(mat.N) 

row_split = rep("group1", 17)
row_split[4:9] = "group2"
row_split[10:14] = "group3"
row_split[15:17] = "group4"

maxvalue = max(mat.expr, na.rm = TRUE)

textsize = 12

### heatmap ----
circlesize =0.2 # change it to adjust the circle size
col_fun = col_fun = circlize::colorRamp2(seq(-maxvalue, maxvalue, length=13), c(dichromat_pal("DarkRedtoBlue.12")(12)[1:6],"white",dichromat_pal("DarkRedtoBlue.12")(12)[7:12]))
cell_heatmap = Heatmap(mat.expr,
              rect_gp = gpar(type = "none"),
              col = col_fun,
              na_col = "grey93",
              row_split = row_split,
              row_title=NULL,
              # right_annotation = row_ha, # if not including the row annotation, commenting out this line
              cluster_columns = F,show_row_names = T,
              cluster_rows = F,show_row_dend = F,show_column_dend = F,row_dend_side = "left",column_dend_side = "top",
              clustering_method_rows = "ward.D2",clustering_method_columns  = "ward.D2",
              clustering_distance_columns  = "euclidean", clustering_distance_rows  = "euclidean",
              row_names_side = "left",
              column_names_side = "top",column_names_rot = 90,column_names_gp = gpar(fontsize = textsize,fontface="plain"),
              left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c("#A781BA","#F47B00","#40A1FF","#80C684")),
                                                               labels = c("", "", "", ""), 
                                                               labels_gp = gpar(col = "white", fontsize = 10))),
              row_labels = c("MARS 1","Coagulopathic","Yao Coagulopathic","MARS 2","SoM Mod 1","SoM Mod 2","Inflammopathic","Yao Innate","SRS Up","MARS 3","Yao Adaptive","Adaptive","SoM Mod 4","SRS Down"," MARS 4", "Wong Score", "SoM Mod 3"),
              row_names_gp = gpar(fontsize = textsize,fontface="plain"),
              column_title = NULL, column_title_side = "top",column_title_gp = gpar(fontsize = textsize,fontface="plain"),
              show_heatmap_legend = F,
              cell_fun =  function(j, i, x, y, width, height, fill){
                grid.rect(x = x, y = y, width = width, height =  height,gp = gpar(lwd=0.2,col = col_fun(mat.expr[i, j]),fill = col_fun(mat.expr[i, j])))},
              heatmap_legend_param=list(title = "scaled expression", title_gp = gpar( fontsize = textsize),labels_gp = gpar(fontsize = textsize),legend_width=unit(0.1,"line"),legend_height=unit(2,"line"),grid_width=unit(1.5,"line"),grid_height=unit(.05,"line"),at = c(-2,2),labels = c("low","high"),by_row = F,direction = "vertical"))

print(cell_heatmap)


#### Looks for cell-specific genes
sc_genes <- sc_data%>%
  mutate(celltype_modified = ifelse(celltype %in% c("Immat Neutrophil"), "Immature Neutrophils", ifelse(celltype == "Neutrophil", "Neutrophils", ifelse(celltype %in% c( "CD14 Monocyte", "CD16 Monocyte"), "Monocytes", ifelse(celltype %in% c("CD4 T", "CD8 T", "Treg","NK","Prolif T/NK"), "T/NKs",NA
  )))))%>%
  filter(!is.na(celltype_modified))%>%
  dplyr::select(celltype_modified,any_of(genes_of_interest))%>%
  group_by(celltype_modified)%>%
  summarise_all(mean, na.rm = TRUE)

mat.expr = as.matrix(sc_genes[,-1])
rownames(mat.expr) = sc_genes$celltype_modified
mat.expr = scale(mat.expr)

#reorder rows into "Immature/Pre Neutrophils", "Neutrophils", "Monocytes", "T/NKs", "Platelets"
mat.expr = mat.expr[c("Immature Neutrophils", "Neutrophils", "Monocytes", "T/NKs"#, "Platelets"
),]

genes <- data.frame(t(mat.expr))%>%
  #if 1 column is >0, returns the name of that column, if more than 1 is, mixed
  mutate(predominant_cell = case_when(
    Immature.Neutrophils - Neutrophils >1 & Immature.Neutrophils - Monocytes >1 & Immature.Neutrophils - T.NKs >1  ~ "Immature/Pre Neutrophils",
    Neutrophils -Immature.Neutrophils >1 & Neutrophils - Monocytes >1 & Neutrophils - T.NKs >1  ~ "Neutrophils",
    Monocytes - Immature.Neutrophils >1 & Monocytes - Neutrophils >1 & Monocytes - T.NKs >1  ~ "Monocytes",
    T.NKs - Immature.Neutrophils >1 & T.NKs - Neutrophils >1 & T.NKs - Monocytes >1  ~ "T/NKs",
    TRUE ~ "Mixed"
  ))%>%
  arrange(predominant_cell)



#Creates the Cell-Specific Figure
sc_genes <- sc_data%>%
  mutate(celltype_modified = ifelse(celltype %in% c("Immat Neutrophil"), "Immature Neutrophils", ifelse(celltype == "Neutrophil", "Neutrophils", ifelse(celltype %in% c( "CD14 Monocyte", "CD16 Monocyte"), "Monocytes", ifelse(celltype %in% c("CD4 T", "CD8 T", "Treg","NK","Prolif T/NK"), "T/NKs",NA
  )))))%>%
  filter(!is.na(celltype_modified))%>%
  dplyr::select(celltype_modified,any_of(genes_of_interest))%>%
  group_by(celltype_modified)%>%
  summarise_all(mean, na.rm = TRUE)

mat.expr = as.matrix(sc_genes[,-1])
rownames(mat.expr) = sc_genes$celltype_modified
mat.expr = scale(mat.expr)
#remove rownames that aren't in myeloid_detrimental_genes
mat.expr = mat.expr[,colnames(mat.expr) %in% c(myeloid_detrimental_genes, neutrophil_protective_genes,monocyte_protective_genes,lymphoid_protective_genes)]

mat.expr = scale(mat.expr)
#adjust row order to "Immature/Pre Neutrophils", "Neutrophils", "Monocytes", "T/NKs"
mat.expr = mat.expr[c("Immature Neutrophils", "Neutrophils", "Monocytes", "T/NKs"),]


mat.expr = mat.expr[,c( "CEP55", "HMMR", "PRC1", "KIF15", "STOM", "CRISP2", "LTF", "LCN2", "OLFM4", "DEFA4", "CEACAM8", "CTSG","AZU1", "GADD45A",  #detrimental
                        "ZDHHC17",  #protective
                        "CAMP", "ARG1", #detrimental
                        "PPP1R12A", #protective
                        "ANXA3", "SLPI", "ORM1", #detrimental
                        "SOS2", "MAPK14", "TXN", #protective
                        "PPL", #detrimental
                        "MME", "ICAM3", "PDE4B", "FAS", #protective
                        "BCL6", "KLHL2", #detrimental
                        "TLR1", "RAF1", #protective
                        "AQP9", #detrimental
                        "GK", "PTEN","PIK3CD","PRKAR1A", "JAK2", "ASAH1", "GNAI3", "SP1", "ATG3", #protective
                        "TYK2", #detrimental
                        "BTK", "SEMA6B", "EP300", "MAP3K3", "PAK2", "CASP1", "CCL2", "MAFB", #protective
                        "HTRA1", #detrimental
                        "BCAT1", "IL1A", "TLE4", "POU2F2", "PLEKHO1", "BCL2L11", "MAP3K1", "PRKCB", "RHBDF2", #protective
                        "SLC1A5", "STX1A", #detrimental
                        "YKT6","VRK2", "CTNNB1", "ITGAV", "VAMP5", "INPP5D", "TBK1", "SLC12A7", "BMPR2", "CREB1", "BTN2A2", "SMAD4", "ZCCHC4", "DBT", "SMYD2", "ARL14EP", "DYRK2", "MAP4K1", "PPP2R5C", "TP53BP1", "JAK1", "BTN3A2", "CASP8", "TRIB2", "PLCG1", "ZNF831", "ZAP70", "CD247", "SIDT1", "NCR3", "KLRB1", "CAMK4", "CD3G", "CD3E", "DDX6", "SEMA4F",  "BPGM", "PIK3R1", "BUB3","CCNB1IP1" #protective
)]

Category = rep("Detrimental", 104)
Category[1:14] = "Detrimental"
Category[15] = "Protective"
Category[16:17] = "Detrimental"
Category[18] = "Protective"
Category[19:21] = "Detrimental"
Category[22:24] = "Protective"
Category[25] = "Detrimental"
Category[26:29] = "Protective"
Category[30:31] = "Detrimental"
Category[32:33] = "Protective"
Category[34] = "Detrimental"
Category[35:43] = "Protective"
Category[44] = "Detrimental"
Category[45:52] = "Protective"
Category[53] = "Detrimental"
Category[54:62] = "Protective"
Category[63:64] = "Detrimental"
Category[65:104] = "Protective"

Category <- factor(Category, levels = c("Detrimental", "Protective"))

maxvalue = max(mat.expr, na.rm = TRUE)
textsize = 12

mat.expr = t(mat.expr)

col_fun = col_fun = circlize::colorRamp2(seq(-maxvalue, maxvalue, length=13), c(dichromat_pal("DarkRedtoBlue.12")(12)[1:6],"white",dichromat_pal("DarkRedtoBlue.12")(12)[7:12]))
gene_heatmap = Heatmap(mat.expr,
              rect_gp = gpar(type = "none"),
              col = col_fun,
              na_col = "grey93",
              # right_annotation = row_ha, # if not including the row annotation, commenting out this line
              cluster_columns = F,#show_row_names = T,
              show_column_names = T,
              row_names_side = "left",
              left_annotation = rowAnnotation(df = data.frame(Category), col = list(Category = c("Detrimental" = "#a97c50", "Protective" = "#9cd6c3")), show_legend = F, show_annotation_name = F),
              cluster_rows = F,show_row_dend = F,show_column_dend = F,row_dend_side = "left",column_dend_side = "top",
              clustering_method_rows = "ward.D2",clustering_method_columns  = "ward.D2",
              clustering_distance_columns  = "euclidean", clustering_distance_rows  = "euclidean",
              column_names_side = "top",column_names_rot = 45,column_names_gp = gpar(fontsize = textsize,fontface="plain"),
              # row_names_gp = gpar(fontsize = textsize,fontface="italic"),
              column_title_side = "bottom",column_title_gp = gpar(fontsize = textsize,fontface="plain"),
              show_heatmap_legend = F,
              cell_fun =  function(j, i, x, y, width, height, fill){
                grid.rect(x = x, y = y, width = width, height =  height,gp = gpar(lwd=0.2,col = "grey90",fill = col_fun(mat.expr[i, j])))},
              heatmap_legend_param=list(title = "scaled expression", title_gp = gpar( fontsize = textsize),labels_gp = gpar(fontsize = textsize),legend_width=unit(0.1,"line"),legend_height=unit(2,"line"),grid_width=unit(1.5,"line"),grid_height=unit(.05,"line"),at = c(-2,2),labels = c("low","high"),by_row = F,direction = "vertical"))


print(gene_heatmap)
print(pht)