#Required Packages
library(tidyverse)
library(DESeq2)

#Pulls and cleans cchmc expression data
imx_expr <- read.csv("~/imx_expr_tpm.csv")

imx_expr <- imx_expr[,-c(1:2)]%>%
  group_by(HGNC_SYM)%>%
  summarise(across(everything(), list(mean))) %>%
  filter(HGNC_SYM != "")%>%
  #mutate all but first column to numerics
  mutate_at(vars(-("HGNC_SYM")), as.numeric)

imx_expr_matrix <- as.matrix(imx_expr[,-1])

rownames(imx_expr_matrix) <- imx_expr$HGNC_SYM

#creates a matrix called coldata that uses the column names of imx_expr_matrix as the rownames and has a column called "pid" with the value "cchmc" in all rows
coldata_imx <- data.frame(t(imx_expr_matrix))%>%
  mutate(pid = "imx",condition = "healthy")%>%
  dplyr::select(pid,condition)



#Pulls and cleans cchmc expression data
cchmc_expr <- read.csv("~/cchmc_expr_tpm.csv")

cchmc_expr <- cchmc_expr[,-c(1:2)]%>%
  group_by(HGNC_SYM)%>%
  summarise(across(everything(), list(mean))) %>%
  filter(HGNC_SYM != "")%>%
  #mutate all but first column to numerics
  mutate_at(vars(-("HGNC_SYM")), as.numeric)

cchmc_expr_matrix <- as.matrix(cchmc_expr[,-1])

rownames(cchmc_expr_matrix) <- cchmc_expr$HGNC_SYM

#creates a matrix called coldata_cchmc that uses the column names of cchmc_expr_matrix as the rownames and has a column called "pid" with the value "cchmc" in all rows
coldata_cchmc <- data.frame(t(cchmc_expr_matrix))%>%
  mutate(pid = "cchmc",condition = "sepsis")%>%
  dplyr::select(pid,condition)


savemore_expr <- read.csv("~/savemore_expr_tpm.csv")

savemore_expr <- savemore_expr[,-c(1:2)]%>%
  group_by(HGNC_SYM)%>%
  summarise(across(everything(), list(mean))) %>%
  filter(HGNC_SYM != "")

savemore_expr_matrix <- as.matrix(savemore_expr[,-1])

rownames(savemore_expr_matrix) <- savemore_expr$HGNC_SYM

coldata_savemore <- data.frame(t(savemore_expr_matrix))%>%
  mutate(pid = "savemore",condition = "covid-19")%>%
  dplyr::select(pid,condition)

#Trinity Data
trinity_expr <- read.csv("~/trinity_expr_tpm.csv")

trinity_expr <- trinity_expr[,-c(1:2)]%>%
  group_by(HGNC_SYM)%>%
  summarise(across(everything(), list(mean))) %>%
  filter(HGNC_SYM != "")%>%
  #mutate all but first column to numerics
  mutate_at(vars(-("HGNC_SYM")), as.numeric)

trinity_expr_matrix <- as.matrix(trinity_expr[,-1])

rownames(trinity_expr_matrix) <- trinity_expr$HGNC_SYM

#creates a matrix called coldata_cchmc that uses the column names of cchmc_expr_matrix as the rownames and has a column called "pid" with the value "cchmc" in all rows
coldata_trinity <- data.frame(t(trinity_expr_matrix))%>%
  mutate(pid = "trinity",condition = "sepsis")%>%
  dplyr::select(pid,condition)


#Pulls the Charles data
charles_expr <- read.csv("~/charles_expr_tpm.csv")

charles_expr <- charles_expr[,-c(1:2)]%>%
  group_by(HGNC_SYM)%>%
  summarise(across(everything(), list(mean))) %>%
  filter(HGNC_SYM != "")%>%
  #mutate all but first column to numerics
  mutate_at(vars(-("HGNC_SYM")), as.numeric)

charles_expr_matrix <- as.matrix(charles_expr[,-1])

rownames(charles_expr_matrix) <- charles_expr$HGNC_SYM


#creates a matrix called coldata_cchmc that uses the column names of cchmc_expr_matrix as the rownames and has a column called "pid" with the value "cchmc" in all rows
coldata_charles <- data.frame(t(charles_expr_matrix))%>%
  mutate(pid = "charles",condition = "sepsis")%>%
  dplyr::select(pid,condition)

#pulls the amsterdam data
amsterdam_expr <- read.csv("~/amsterdam_expr_tpm.csv")

amsterdam_expr <- amsterdam_expr[,-c(1:2)]%>%
  group_by(HGNC_SYM)%>%
  summarise(across(everything(), list(mean))) %>%
  filter(HGNC_SYM != "")%>%
  #mutate all but first column to numerics
  mutate_at(vars(-("HGNC_SYM")), as.numeric)

amsterdam_expr_matrix <- as.matrix(amsterdam_expr[,-1])

rownames(amsterdam_expr_matrix) <- amsterdam_expr$HGNC_SYM

#creates a the coldata
coldata_amsterdam <- data.frame(t(amsterdam_expr_matrix))%>%
  mutate(pid = "amsterdam",condition = "sepsis")%>%
  dplyr::select(pid,condition)

#pulls the stanford data
stanford_expr <- read.csv("~/stanford_expr_tpm.csv")

stanford_expr <- stanford_expr[,-c(1:2)]%>%
  group_by(HGNC_SYM)%>%
  summarise(across(everything(), list(mean))) %>%
  filter(HGNC_SYM != "")%>%
  #mutate all but first column to numerics
  mutate_at(vars(-("HGNC_SYM")), as.numeric)

stanford_expr_matrix <- as.matrix(stanford_expr[,-1])

rownames(stanford_expr_matrix) <- stanford_expr$HGNC_SYM

#creates a the coldata
coldata_stanford <- data.frame(t(stanford_expr_matrix))%>%
  mutate(pid = "stanford",condition = "sepsis")%>%
  dplyr::select(pid,condition)


#pulls the ufl data
ufl_expr <- read.csv("~/ufl_expr_tpm.csv")

ufl_expr <- ufl_expr[,-c(1:2)]%>%
  group_by(HGNC_SYM)%>%
  summarise(across(everything(), list(mean))) %>%
  filter(HGNC_SYM != "")%>%
  #mutate all but first column to numerics
  mutate_at(vars(-("HGNC_SYM")), as.numeric)

ufl_expr_matrix <- as.matrix(ufl_expr[,-1])

rownames(ufl_expr_matrix) <- ufl_expr$HGNC_SYM

#creates a the coldata
coldata_ufl <- data.frame(t(ufl_expr_matrix))%>%
  mutate(pid = "victas",condition = "sepsis/controls")%>%
  dplyr::select(pid,condition)


#pulls the victas data
victas_expr <- read.csv("~/victas_expr_tpm.csv")

victas_expr <- victas_expr[,-c(1:2)]%>%
  group_by(HGNC_SYM)%>%
  summarise(across(everything(), list(mean))) %>%
  filter(HGNC_SYM != "")%>%
  #mutate all but first column to numerics
  mutate_at(vars(-("HGNC_SYM")), as.numeric)

victas_expr_matrix <- as.matrix(victas_expr[,-1])

rownames(victas_expr_matrix) <- victas_expr$HGNC_SYM

#creates a the coldata
coldata_victas <- data.frame(t(victas_expr_matrix))%>%
  mutate(pid = "victas",condition = "sepsis")%>%
  dplyr::select(pid,condition)


common_rows <- Reduce(intersect, list(rownames(amsterdam_expr_matrix), rownames(cchmc_expr_matrix), rownames(charles_expr_matrix), rownames(imx_expr_matrix) ,rownames(savemore_expr_matrix), rownames(stanford_expr_matrix), rownames(trinity_expr_matrix), rownames(victas_expr_matrix),rownames(ufl_expr_matrix)))
amsterdam_expr_matrix <- amsterdam_expr_matrix[common_rows,]
cchmc_expr_matrix <- cchmc_expr_matrix[common_rows,]
charles_expr_matrix <- charles_expr_matrix[common_rows,]
imx_expr_matrix <- imx_expr_matrix[common_rows,]
savemore_expr_matrix <- savemore_expr_matrix[common_rows,]
stanford_expr_matrix <- stanford_expr_matrix[common_rows,]
trinity_expr_matrix <- trinity_expr_matrix[common_rows,]
ufl_expr_matrix <- ufl_expr_matrix[common_rows,]
victas_expr_matrix <- victas_expr_matrix[common_rows,]

pre_norm_data <- cbind(amsterdam_expr_matrix, cchmc_expr_matrix, charles_expr_matrix, imx_expr_matrix, savemore_expr_matrix, stanford_expr_matrix, trinity_expr_matrix,ufl_expr_matrix,victas_expr_matrix)
pre_norm_data <- matrix(data = as.integer(pre_norm_data), nrow = nrow(pre_norm_data), dimnames = list(rownames(pre_norm_data), colnames(pre_norm_data)))
pre_norm_data <- pre_norm_data + 1

coldata <- as.matrix(rbind(coldata_amsterdam, coldata_cchmc, coldata_charles, coldata_imx, coldata_savemore, coldata_stanford, coldata_trinity,coldata_ufl,coldata_victas))

dds <- DESeqDataSetFromMatrix(countData = pre_norm_data,
                              colData = coldata,
                              design = ~ pid)

saveRDS(dds,file= "~/tpm_deseq_with_controls.rds")

normalized_data <- varianceStabilizingTransformation(dds, blind = FALSE, fitType = "parametric")

saveRDS(normalized_data,file= "~/tpm_normalized_data_with_controls.rds")

#limma batch effect removal the vst transformed data 
limma_normalized_data <- normalized_data
assay(limma_normalized_data) <- limma::removeBatchEffect(assay(limma_normalized_data), limma_normalized_data$pid)
hgnc_expr_limma <- assay(limma_normalized_data)

saveRDS(hgnc_expr_limma, "~/tpm_hgnc_expr_limma_with_controls.RDS")
