if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("genefu")) {
  BiocManager::install("genefu", ask =FALSE)
  library("genefu")
}
library(genefu) ;library(xtable) ; library(rmeta) ; library(Biobase) ; library(caret)
library("breastCancerMAINZ")

data(breastCancerData)
cinfo <- colnames(pData(mainz7g))
data.all <- c("transbig7g"=transbig7g, "unt7g"=unt7g, "upp7g"=upp7g,
              "mainz7g"=mainz7g, "nki7g"=nki7g)
idtoremove.all <- NULL
duplres <- NULL


library
(breastCancerMAINZ)
library
(breastCancerTRANSBIG)
library
(breastCancerUPP)
library
(breastCancerUNT)
library
(breastCancerNKI
