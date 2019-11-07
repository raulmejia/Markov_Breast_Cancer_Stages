gProfileR::gconvert()
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("gProfileR")) {
  BiocManager::install("gProfileR", ask =FALSE)
  library("gProfileR")
}
#####
# Data given by the user
###
args <- commandArgs(trailingOnly = TRUE)
Path_to_your_Matrix<-args[1] # The path to your matrix
Path_to_your_Matrix <- c("/////")
args[2]
args[3]

#####
# reading the data
#####

gene_list_raw <- read.table(Path_to_your_gene_list, header = FALSE, sep="\t", stringsAsFactors = FALSE) 

xprs_and_mirnas_list<-readRDS("/media/rmejia/ADATA/Transcriptome_Profiling/DataMart_miRNA-RNAseq/Analysis/Results/matrices_by_stage/healthy/xprssn_and_mirnas_list.rds")

str(xprs_and_mirnas_list)

str(xprs_and_mirnas_list[[1]])
xprs_and_mirnas_list[[3]][1:4,1:4]
head(xprs_and_mirnas_list[[1]])
names(xprs_and_mirnas_list)
