#Install Dependencies
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("GEOmetadb")
#install.packages('MetaIntegrator')
library(MetaIntegrator)

#BiocManager::install("batchelor")
#devtools::install_github("jknightlab/SepstratifieR")
library(SepstratifieR)

#install.packages("tidyverse")
library(tidyverse)

#install.packages("data.table")
library(data.table)

#install.packages("nnet")
library(nnet)

##################
#Set-Up
#This will need to be adjusted depending on what directories we will be using. In this case, I have set up a folder subspace and within that I have input and output folders for each individual study

#ID for whichever biobank we're running analysis in. 
sub_id = "study_id"

#subspace directory
dir_main = paste0("~/subspace/")

#directory for input files
dir_in = paste0("~/subspace/",sub_id,"/input/")

#directory for output files
dir_out = paste0("~/subspace/",sub_id,"/output/")

####################
#Reads in data
#requires standardized naming of data. In this case, I used expr, key, pheno
#data should be set up with first three columns being "ENSEMBL_GENE_ID","ENTREZ_GENE_ID","HGNC_SYM", followed by sample_ids
expr <- read.csv(paste0(dir_in,sub_id,"_expr.csv")) 

key <- read.csv(paste0(dir_in,sub_id,"_key.csv"))

pheno<- read.csv(paste0(dir_in,sub_id,"_pheno.csv"))

##################
#QC for expression data
#Checks to see if needs to be transformed
boxplot(expr[,c(4:min(10,dim(expr)[2]))])
boxplot(expr[sample(nrow(expr),1000),-c(1:3)])

#Log2 transforms if needed - double check boxplots above and comment out if not needed
expr[,-c(1:3)] = log2(expr[,-c(1:3)] - min(-c(1:3), na.rm = TRUE) + 1)

#re-checks post transformation if needed
boxplot(expr[,c(4:min(10,dim(expr)[2]))])
boxplot(expr[sample(nrow(expr),1000),-c(1:3)])

dim(expr)
View(expr)
any(is.na(expr[,-c(1:3)]))

##################
#HGNC Gene expression matrix set up
#creates a hgnc dataset and removes blank HGNC symbols
hgnc_expr <- expr[,-c(1:2)]%>%
  group_by(HGNC_SYM)%>%
  summarise(across(everything(), list(mean))) %>%
  filter(HGNC_SYM != "")

#Converts this into a matrix with rownames as the HGNC gene symbols, which is required for getting all gene scores except sepstratifier
hgnc_expr = data.frame(hgnc_expr)
rownames(hgnc_expr) = hgnc_expr[,1]
hgnc_expr = hgnc_expr[,-1]
hgnc_expr = as.matrix(hgnc_expr)

##################
#SEPSRATIFIER SETUP
#Converts to just the Ensembl gene names, converts gene names to the format needed for sepstratifier
ensembl_expr <- expr[,-c(2:3)] %>%
  mutate(ENSMBL_GENE_ID = substr(ENSMBL_GENE_ID, 1,15))  %>%
  group_by(ENSMBL_GENE_ID)%>%
  summarise(across(everything(), list(mean))) %>%
  filter(ENSMBL_GENE_ID != "")

ensembl_expr = data.frame(ensembl_expr)
rownames(ensembl_expr) = ensembl_expr[,1]
ensembl_expr = ensembl_expr[,-1]

#transposes data frame to have the ensembl as column names for input into sepstratifier 
sepstratifier_expr <- t(ensembl_expr)

##################
#Functions
#Geometric Mean Function
geomMean <- function(x, na.rm = FALSE){
  if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("argument is not numeric or logical: returning NA")
    return(as.numeric(NA))
  }
  if (na.rm){
    x <- x[!is.na(x)]
  }
  if (any(x < 0)){
    stop("'x' contains negative value(s)")
  }
  if(any(x == 0)){ #if I don't this, having any 0s will result in the geometric mean being 0
    x = x[-c(which(x == 0))]
    if(length(x) == 0){ #this means all values were 0
      return(0)
    }
  }
  return(exp(sum(log(x))/length(x)))
}

#Generic gene score function. Used for zheng, yao
getGeneScores <- function(geneMtx, pos, neg, makePos = TRUE, out.missing=TRUE){
  if(out.missing){
    missingpos = pos[!(pos %in% rownames(geneMtx))]
    missingneg = neg[!(neg %in% rownames(geneMtx))]
    missing = c(missingpos,missingneg)
    if(length(missing)>0){
      cat("Missing these genes:",missing,"\n")
    }
  }
  pos = pos[pos %in% rownames(geneMtx)]
  neg = neg[neg %in% rownames(geneMtx)]
  
  if(makePos){ #If I need to make everything positive
    if(any(geneMtx < 0, na.rm=T)){
      geneMtx <- geneMtx - min(geneMtx, na.rm=T)
    }
  }
  
  if(length(neg)>=1 && length(pos)>=1){
    scores <- apply(geneMtx[pos, , drop=F], 2, geomMean, na.rm=T) - apply(geneMtx[neg, , drop=F], 2, geomMean, na.rm=T)
  }else if(length(pos)>=1){
    scores <- apply(geneMtx[pos, , drop=F], 2, geomMean, na.rm=T)
  }else if(length(neg)>=1){
    scores <- -1*apply(geneMtx[neg, , drop=F], 2, geomMean, na.rm=T)
  }
  
  if(any(is.nan(scores))){ #this means all upgenes or all downgenes were missing from some samples
    nanIndex = which(is.nan(scores))
    for(i in nanIndex){
      my.score = c(geomMean(geneMtx[pos,i,drop=F],na.rm=T),geomMean(geneMtx[neg,i,drop=F],na.rm=T))
      my.score = my.score[!is.nan(my.score)]
      if(length(my.score) == 1){
        scores[i] = my.score
      }else if(length(my.score) == 2){
        scores[i] = my.score[1] - my.score[2]
      }
    }
  }
  
  return(scores)
}

##################
#Calculates Zheng (SoM) Scores
#SoM modules
mod1 = c("NQO2","SLPI","ORM1","KLHL2","ANXA3","TXN","AQP9","BCL6","DOK3","PFKFB4","TYK2")
mod2 = c("BCL2L11","BCAT1","BTBD7","CEP55","HMMR","PRC1","KIF15","CAMP","CEACAM8","DEFA4","LCN2","CTSG","AZU1")
mod3 = c("MAFB","OASL","UBE2L6","VAMP5","CCL2","NAPA","ATG3","VRK2","TMEM123","CASP7")
mod4 = c("DOK2","HLA-DPB1","BUB3","SMYD2","SIDT1","EXOC2","TRIB2","KLRB1")

#calculate module scores
mod1_score = getGeneScores(hgnc_expr,pos=mod1,neg=c())
mod2_score = getGeneScores(hgnc_expr,pos=mod2,neg=c())
mod3_score = getGeneScores(hgnc_expr,pos=mod3,neg=c())
mod4_score = getGeneScores(hgnc_expr,pos=mod4,neg=c())
#calculate combined score
som_score = (mod1_score+mod2_score)/(mod3_score+mod4_score)
protective_score = (mod3_score+mod4_score)
detrimental_score = (mod1_score+mod2_score)

som_score_table<-data.table(names(som_score),som_score) %>%
  setnames("V1", "accession") 

protective_score_table<-data.table(names(protective_score),protective_score) %>%
  setnames("V1", "accession")

detrimental_score_table<-data.table(names(detrimental_score),detrimental_score) %>%
  setnames("V1", "accession")

mod1_table<-data.table(names(mod1_score),mod1_score) %>%
  setnames("V1", "accession")

mod2_table<-data.table(names(mod2_score),mod2_score)%>%
  setnames("V1", "accession")

mod3_table<-data.table(names(mod3_score),mod3_score) %>%
  setnames("V1", "accession")

mod4_table<-data.table(names(mod4_score),mod4_score) %>%
  setnames( "V1", "accession")

zheng_table = full_join(mod1_table,mod2_table, by = "accession") %>%
  full_join(.,mod3_table, by = "accession")%>%
  full_join(.,mod4_table, by = "accession")%>%
  full_join(.,detrimental_score_table, by = "accession")%>%
  full_join(.,protective_score_table, by = "accession")%>%
  full_join(.,som_score_table, by = "accession")%>%
  mutate(zheng_endotype = ifelse(som_score >= 1, "detrimental","protective"))

##################
#Sweeney Endotypes
#the endotype genes
inflammopathic.up= c("ARG1","LCN2","LTF","OLFM4")
inflammopathic.down = c("HLA-DMB")
adaptive.up = c("YKT6", "PDE4B", "TWISTNB", "BTN2A2", "ZBTB33", "PSMB9", "CAMK4", "TMEM19", "SLC12A7", 
                "TP53BP1", "PLEKHO1", "SLC25A22", "FRS2")                          
adaptive.down = c("GADD45A", "CD24", "S100A12", "STX1A")
coagulopathic.up = c("KCNMB4", "CRISP2", "HTRA1","PPL")
coagulopathic.down = c("RHBDF2", "ZCCHC4", "YKT6", "DDX6", "SENP5", "RAPGEF1","DTX2","RELB")

#The function
getGeneScores_sweeney <- function(geneMtx, pos, neg, makePos = TRUE, out.missing=TRUE){
  if(out.missing){
    missingpos = pos[!(pos %in% rownames(geneMtx))]
    missingneg = neg[!(neg %in% rownames(geneMtx))]
    missing = c(missingpos,missingneg)
    if(length(missing)>0){
      cat("Missing these genes:",missing,"\n")
    }
  }
  pos = pos[pos %in% rownames(geneMtx)]
  neg = neg[neg %in% rownames(geneMtx)]
  
  if(makePos){ #If I need to make everything positive
    if(any(geneMtx < 0, na.rm=T)){
      geneMtx <- geneMtx - min(geneMtx, na.rm=T)
    }
  }
  
  if(length(neg)>=1 && length(pos)>=1){
    scores <- apply(geneMtx[pos, , drop=F], 2, geomMean, na.rm=T) - apply(geneMtx[neg, , drop=F], 2, geomMean, na.rm=T)
  }else if(length(pos)>=1){
    scores <- apply(geneMtx[pos, , drop=F], 2, geomMean, na.rm=T)
  }else if(length(neg)>=1){
    scores <- -1*apply(geneMtx[neg, , drop=F], 2, geomMean, na.rm=T)
  }
  
  if(any(is.nan(scores))){ #this means all upgenes or all downgenes were missing from some samples
    nanIndex = which(is.nan(scores))
    for(i in nanIndex){
      my.score = c(geomMean(geneMtx[pos,i,drop=F],na.rm=T),geomMean(geneMtx[neg,i,drop=F],na.rm=T))
      my.score = my.score[!is.nan(my.score)]
      if(length(my.score) == 1){
        scores[i] = my.score
      }else if(length(my.score) == 2){
        scores[i] = my.score[1] - my.score[2]
      }
    }
  }
  
  return(scale(scores))
}

#Gets the scores
score1 = getGeneScores_sweeney(hgnc_expr,pos=inflammopathic.up,neg=inflammopathic.down)
score2 = getGeneScores_sweeney(hgnc_expr,pos=adaptive.up,neg=adaptive.down)
score3= getGeneScores_sweeney(hgnc_expr,pos=coagulopathic.up,neg=coagulopathic.down)

#Puts them into a dataframes
sweeney_table = data.frame(accession = rownames(score1),score1,score2,score3) 

#pulls the LR model for endotyping
multLR<-readRDS(file=paste0(dir_main,"multLR_model.RDS"))

probs <- data.frame(endoProbs=round(predict(multLR, sweeney_table, type="probs"), 3),
                    endotype=predict(multLR, sweeney_table, type="class"))

probs$accession = rownames(probs)

sweeney_table  <- sweeney_table  %>%
  left_join(.,probs, by = "accession")%>%
  mutate(inflammopathic_score = score1,
         adaptive_score = score2,
         coagulopathic_score = score3,
         inflammopathic_prob = endoProbs.1,
         adaptive_prob = endoProbs.2,
         coagulopathic_prob = endoProbs.3,
         sweeney_endotype = ifelse(endotype == 1, "inflammopathic",
                                   ifelse(endotype == 2, "adaptive",
                                          ifelse(endotype == 3,"coagulopathic",NA))))%>%
  select(-c("score1","score2","score3","endoProbs.1","endoProbs.2","endoProbs.3","endotype"))

##################
#Davenport and Cano-Gamez Scores using Sepstratifier
davenport_scores <- stratifyPatients(sepstratifier_expr,gene_set = "davenport", k = 20)

cano_scores <- stratifyPatients(sepstratifier_expr,gene_set = "extended", k = 20)

davenport_table<-data.frame(accession = names(davenport_scores@SRS),davenport_endotype = davenport_scores@SRS,davenport_SRSq = davenport_scores@SRSq,davenport_SRS1_prob = davenport_scores@SRS_probs$SRS1,davenport_SRS2_prob =davenport_scores@SRS_probs$SRS2,davenport_SRS3_prob =davenport_scores@SRS_probs$SRS3) 

cano_table<-data.frame(accession = names(cano_scores@SRS),cano_endotype = cano_scores@SRS,cano_SRSq = cano_scores@SRSq,cano_SRS1_prob = cano_scores@SRS_probs$SRS1,cano_SRS2_prob =cano_scores@SRS_probs$SRS2,cano_SRS3_prob =cano_scores@SRS_probs$SRS3) 

srs_table = full_join(davenport_table,cano_table, by = "accession")

##################
#Yao Score
#Yao Score genes
yao_IA_up <- c("ZNF831", "CD3G", "MME", "BTN3A2", "HLA-DPA1")
yao_IA_down <- "STOM"
yao_IN_up <- c("HK3","SERPINB1")
yao_IN_down <- c("EPB42", "GSPT1", "LAT")
yao_IC_up <- c("SLC1A5", "IGF2BP2", "ANXA3")
yao_IC_down <- "GBP2"

#calculate yao scores
yao_IA_score = getGeneScores(hgnc_expr,pos=yao_IA_up,neg = yao_IA_down)
yao_IN_score = getGeneScores(hgnc_expr,pos=yao_IN_up,neg = yao_IN_down)
yao_IC_score = getGeneScores(hgnc_expr,pos=yao_IC_up,neg= yao_IC_down)

#Creates gene score tables and combines them
yao_IA_table<-data.table(names(yao_IA_score),yao_IA_score) %>%
  setnames("V1", "accession")

yao_IN_table<-data.table(names(yao_IN_score),yao_IN_score)%>%
  setnames("V1", "accession")

yao_IC_table<-data.table(names(yao_IC_score),yao_IC_score) %>%
  setnames("V1", "accession")

yao_table = full_join(yao_IA_table,yao_IC_table, by = "accession") %>%
  full_join(.,yao_IN_table, by = "accession")

##################
#Wong Score
#is a geometric mean of all 100 genes
wong_genes = c("APAF1", "ARPC5","ASAH1","ATP2B2","BCL6","BMPR2","BTK","CAMK2D","CAMK2G","CAMK4","CASP1","CASP2","CASP4","CASP8","CD247","CD3E","CD3G","CD79A",
               "CREB1","CREB5","CSNK1A1","CTNNB1","DAPP1","DBT","EP300","FAS","FCGR2A","FCGR2C","FYN","GK","GNAI3","HDAC4","HLA-DMA","HLA-DOA","ICAM3","IL1A",
               "INPP5D","ITGAM","ITGAV","ITGAX","JAK1","JAK2","KAT2B","LAT2","LYN","MAP2K4","MAP3K1","MAP3K3","MAP3K5","MAP3K7","MAP4K1","MAP4K4", "MAPK1","MAPK14","MDH1","MKNK1",
               "NCOA2","NCR3","NFATC1","PAK2","PDPR","PIAS1","PIK3C2A","PIK3C3","PIK3CA","PIK3CD","PIK3R1","PLCG1","POU2F2","PPP1R12A","PPP2R2A","PPP2R5C",
               "PRKAR1A","PRKCB","PSMB7","PTEN","PTPRC","RAF1","RHOT1","ROCK1","SEMA4F","SEMA6B","SMAD4","SOS1","SOS2","SP1","TAF11","TBK1","TGFBR1","TLE4",
               "TLR1","TLR2","TLR8","TNFSF10","TRA@","TYROBP","UBE3A","USP48","ZAP70","ZDHHC17")


wong_expr <- 2^hgnc_expr[(row.names(hgnc_expr) %in% wong_genes),]


wong_missing <- function(geneMtx, genes, out.missing=TRUE){
  if(out.missing){
    missinggenes = genes[!(genes %in% rownames(geneMtx))]
    missing = c(missinggenes)
    if(length(missing)>0){
      cat("Missing these genes:",missing,"\n")
    }
  }
}

#tells you which genes are missing 
wong_missing(wong_expr,genes = wong_genes)

wong_int1 <- data.frame(t(wong_expr)) %>%
  mutate(geomean = apply(wong_int1, 1, FUN = geomMean, na.rm = TRUE))

wong_table <- (wong_int1 - wong_int1[,"geomean"])^2 %>%
  mutate(wong_score = rowSums(.)/1000000,
         accession = rownames(wong_int2))%>%
  select("accession","wong_score")


##################
#MARS SCORES
#Identifies the MARS genes
mars1_up = c("BPGM")
mars1_down = c("TAP2")
mars2_up = "GADD45A"
mars2_down = "PCGF5"
mars3_up = "AHNAK"
mars3_down = "PDCD10"
mars4_up = "IFIT5"
mars4_down = "NOP53"

#MARS Score function
getGeneScores_mars <- function(geneMtx, pos, neg, makePos = TRUE, out.missing=TRUE){
  if(out.missing){
    missingpos = pos[!(pos %in% rownames(geneMtx))]
    missingneg = neg[!(neg %in% rownames(geneMtx))]
    missing = c(missingpos,missingneg)
    if(length(missing)>0){
      cat("Missing these genes:",missing,"\n")
    }
  }
  pos = pos[pos %in% rownames(geneMtx)]
  neg = neg[neg %in% rownames(geneMtx)]
  
  if(makePos){ #If I need to make everything positive
    if(any(geneMtx < 0, na.rm=T)){
      geneMtx <- geneMtx - min(geneMtx, na.rm=T)
    }
  }
  
  if(length(neg)>=1 && length(pos)>=1){
    scores <- geneMtx[pos, , drop=F] / geneMtx[neg, , drop=F]
  }else if(length(pos)>=1){
    scores <- NA
  }else if(length(neg)>=1){
    scores <- -NA
  }
  
  return(scores)
}

mars1_score = getGeneScores_mars(hgnc_expr,pos=mars1_up,mars1_down)
mars1_table<-data.frame(accession = colnames(mars1_score), mars1_score = mars1_score[1,])%>%
  mutate(cluster_1 = ifelse(mars1_score >= 1.15, "MARS1",NA))

mars2_score = getGeneScores_mars(hgnc_expr,pos=mars2_up,mars2_down)
mars2_table<-data.frame(accession = colnames(mars2_score), mars2_score = mars2_score[1,]) %>%
  mutate(cluster_2 = ifelse(mars2_score >= 1.05, "MARS2",NA))

mars3_score = getGeneScores_mars(hgnc_expr,pos=mars3_up,mars3_down)
mars3_table<-data.frame(accession = colnames(mars3_score), mars3_score = mars3_score[1,]) %>%
  mutate(cluster_3 = ifelse(mars3_score >= 1.03, "MARS3",NA))

mars4_score = getGeneScores_mars(hgnc_expr,pos=mars4_up,mars4_down)
mars4_table<-data.frame(accession = colnames(mars4_score), mars4_score = mars4_score[1,]) %>%
  mutate(cluster_4 = ifelse(mars4_score >= 0.57, "MARS4",NA))

mars_table<- full_join(mars1_table, mars2_table, by = "accession") %>%
  full_join(., mars3_table, by = "accession") %>%
  full_join(., mars4_table, by = "accession") %>%
  unite(.,col = mars_endotype,cluster_1,cluster_2,cluster_3,cluster_4, sep = "/",remove = TRUE,na.rm = TRUE)


##################
myeloid_detrimental_genes <- c('ANXA3', 'ARG1', 'AZU1', 'CAMP', 'CEACAM8', 'CEP55', 'CRISP2', 'CTSG', 'DEFA4', 'GADD45A', 'HMMR', 'KIF15', 'LCN2', 'LTF', 'OLFM4', 'ORM1', 'PRC1', 'SLPI', 'STOM', 'AQP9', 'BCL6', 'KLHL2', 'PPL', 'HTRA1', 'TYK2', 'SLC1A5', 'STX1A')

neutrophil_protective_genes <- c('ZDHHC17', 'FAS', 'GK', 'ICAM3', 'MME', 'PDE4B', 'PIK3CD', 'PTEN', 'RAF1', 'TLR1', 'PPP1R12A', 'MAPK14', 'SOS2', 'TXN')

monocyte_protective_genes <- c('ASAH1', 'ATG3', 'BCAT1', 'BCL2L11', 'BTK', 'BTN2A2', 'CASP1', 'CCL2', 'CREB1', 'EP300', 'GNAI3', 'IL1A', 'JAK2', 'MAFB', 'MAP3K1', 'MAP3K3', 'PAK2', 'PLEKHO1', 'POU2F2', 'PRKAR1A', 'PRKCB', 'RHBDF2', 'SEMA6B', 'SP1', 'TLE4', 'BMPR2', 'CTNNB1', 'INPP5D', 'ITGAV', 'SLC12A7', 'TBK1', 'VAMP5', 'VRK2', 'YKT6')

lymphoid_protective_genes <- c('ARL14EP', 'BPGM', 'BTN3A2', 'BUB3', 'CAMK4', 'CASP8', 'CCNB1IP1', 'CD247', 'CD3E', 'CD3G', 'DBT', 'DDX6', 'DYRK2', 'JAK1', 'KLRB1', 'MAP4K1', 'NCR3', 'PIK3R1', 'PLCG1', 'PPP2R5C', 'SEMA4F', 'SIDT1', 'SMAD4', 'SMYD2', 'TP53BP1', 'TRIB2', 'ZAP70', 'ZCCHC4', 'ZNF831')


myeloid_score = getGeneScores(hgnc_expr,pos=c(myeloid_detrimental_genes),neg = c(neutrophil_protective_genes, monocyte_protective_genes))
myeloid_table <- data.table(names(myeloid_score),myeloid_score)%>%
  setnames("V1","accession")

myeloid_detrimental_score = getGeneScores(hgnc_expr,pos=c(myeloid_detrimental_genes),neg = c())
myeloid_detrimental_table <- data.table(names(myeloid_detrimental_score),myeloid_detrimental_score)%>%
  setnames("V1","accession")

neutrophil_protective_score = getGeneScores(hgnc_expr,pos=neutrophil_protective_genes,neg = c())
neutrophil_protective_table <- data.table(names(neutrophil_protective_score),neutrophil_protective_score)%>%
  setnames("V1","accession")

monocyte_protective_score = getGeneScores(hgnc_expr,pos=monocyte_protective_genes,neg = c())
monocyte_protective_table <- data.table(names(monocyte_protective_score),monocyte_protective_score)%>%
  setnames("V1","accession")

myeloid_protective_score = getGeneScores(hgnc_expr,pos=c(neutrophil_protective_genes, monocyte_protective_genes),neg = c())
myeloid_protective_table <- data.table(names(myeloid_protective_score),myeloid_protective_score)%>%
  setnames("V1","accession")

lymphoid_score = getGeneScores(hgnc_expr,pos=c(),neg = c(lymphoid_protective_genes))
lymphoid_table <- data.table(names(lymphoid_score),lymphoid_score)%>%
  setnames("V1","accession")

lymphoid_protective_score = getGeneScores(hgnc_expr,pos=lymphoid_protective_genes,neg = c())
lymphoid_protective_table <- data.table(names(lymphoid_protective_score),lymphoid_protective_score)%>%
  setnames("V1","accession")

consensus_scores = full_join(myeloid_detrimental_table,neutrophil_protective_table, by = "accession") %>%
  full_join(.,monocyte_protective_table, by = "accession") %>%
  full_join(.,myeloid_protective_table, by = "accession")%>%
  full_join(.,lymphoid_protective_table, by = "accession") %>%
  full_join(.,lymphoid_table, by = "accession") %>%
  full_join(.,myeloid_table, by = "accession") 

##################
#r pools the score tables
score_table <- full_join(zheng_table, sweeney_table, by = "accession") %>%
  full_join(.,yao_table, by = "accession")%>%
  full_join(.,srs_table, by = "accession")%>%
  full_join(.,wong_table, by = "accession") %>%
  full_join(.,mars_table, by = "accession")%>%
  full_join(.,consensus_scores, by = "accession")

write.csv(score_table, file = paste0(dir_out,sub_id,"_score_table.csv"))
