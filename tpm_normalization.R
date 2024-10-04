# packages and functions
library(dplyr)

# counts to TPM is based on Michael's Love answer at https://support.bioconductor.org/p/91218/

cnt2tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
} 

#ID for whichever data we're running analysis on
sub_id = "sub_id"

#subspace directory
dir_main = paste0("~/subspace/")

#directory for input files
dir_in = paste0("~/subspace/",sub_id,"/input/")

#directory for output files
dir_out = paste0("~/subspace/",sub_id,"/output/")

# Data load
data_counts = read.csv(paste0(dir_in,sub_id,"_expr.csv"))%>%
  dplyr::select(ENSMBL_GENE_ID = X, everything())
data_gene_length = read.csv(paste0(dir_main,"gene_code.csv"))
#Parameter definition
my_missing=c(NA,""," ")

#
head(data_counts[,1:5])
str(data_counts[,1:6])

table(data_counts$ENSMBL_GENE_ID %in% data_gene_length$STARid)
data_counts_full = merge(data_gene_length,data_counts,
                         by.x = "STARid", by.y = "ENSMBL_GENE_ID", all = T)
head(data_counts_full[,1:7])

# Create tables for each ID type separately, such that there are no
# missing or duplicated ids
### STARid
counts_starid = data_counts_full %>%
  dplyr::select (-c(2,4,5,6,7))
#removes any rows with NAs across all columns
counts_starid = counts_starid[complete.cases(counts_starid),]
head(counts_starid[,1:5])
table(duplicated(counts_starid$STARid))
table(counts_starid$STARid %in% my_missing)
counts_starid = counts_starid[!duplicated(counts_starid$STARid),]

# convert to TPM
## STARid
tmp = counts_starid[,3:ncol(counts_starid)]
row.names(tmp) = counts_starid$STARid
counts_starid_tpm = as.data.frame.matrix(cnt2tpm(tmp,counts_starid$ExonLength))%>%
  mutate(STARid = row.names(.))%>%
  left_join(.,data_gene_length, by = "STARid")%>%
  dplyr::select(ENSMBL_GENE_ID = ENSEMBL,ENTREZ_GENE_ID = ENTREZID, HGNC_SYM = SYMBOL, everything(),-c(ExonLength,STARid))
head(counts_starid_tpm[,1:5])

write.csv(counts_starid_tpm,paste0(dir_in,sub_id,"_expr_tpm.csv"),quote = F, row.names=FALSE)
