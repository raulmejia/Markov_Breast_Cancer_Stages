if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("genefu")) {
  BiocManager::install("genefu", ask =FALSE)
  library("genefu")
}
if (!require("breastCancerTRANSBIG")) {
  BiocManager::install("breastCancerTRANSBIG", ask =FALSE)
  library("breastCancerTRANSBIG")
}
if (!require("breastCancerUNT")) {
  BiocManager::install("breastCancerUNT", ask =FALSE)
  library("breastCancerUNT")
}
if (!require("breastCancerUPP")) {
  BiocManager::install("breastCancerUPP", ask =FALSE)
  library("breastCancerUPP")
}
if (!require("caret")) {
  BiocManager::install("caret", ask =FALSE)
  library("caret")
}

###############################
##### Matrix given by the user
###############################





library(genefu) ;library(xtable) ; library(rmeta) ; library(Biobase) ; library(caret) ; library("breastCancerMAINZ")
library(breastCancerMAINZ); library(breastCancerTRANSBIG) ; library(breastCancerUPP);library(breastCancerUNT) ; library(breastCancerNKI)

data(breastCancerData)
cinfo <- colnames(pData(mainz7g))
data.all <- c("transbig7g"=transbig7g, "unt7g"=unt7g, "upp7g"=upp7g,
              "mainz7g"=mainz7g, "nki7g"=nki7g)
idtoremove.all <- NULL
duplres <- NULL
## No overlaps in the MainZ and NKI datasets.
## Focus on UNT vs UPP vs TRANSBIG
demo.all <- rbind(pData(transbig7g), pData(unt7g), pData(upp7g))
dn2 <- c("TRANSBIG", "UNT", "UPP")
## Karolinska
## Search for the VDXKIU, KIU, UPPU series
ds2 <- c("VDXKIU", "KIU", "UPPU")
demot <- demo.all[complete.cases(demo.all[ , c("series")]) & is.element(demo.all[ , "series"], ds2), ]
# Find the duplicated patients in that series
duplid <- sort(unique(demot[duplicated(demot[ , "id"]), "id"]))
duplrest <- NULL
for(i in 1:length(duplid)) {
  tt <- NULL
  for(k in 1:length(dn2)) {
    myx <- sort(row.names(demot)[complete.cases(demot[ , c("id", "dataset")]) &
                                   demot[ , "id"] == duplid[i] & demot[ , "dataset"] == dn2[k]])
    if(length(myx) > 0) { tt <- c(tt, myx) }
  }
  duplrest <- c(duplrest, list(tt))
}
names(duplrest) <- duplid
duplres <- c(duplres, duplrest)
## Oxford
## Search for the VVDXOXFU, OXFU series
ds2 <- c("VDXOXFU", "OXFU")
demot <- demo.all[complete.cases(demo.all[ , c("series")]) & is.element(demo.all[ , "series"], ds2), ]
# Find the duplicated patients in that series
duplid <- sort(unique(demot[duplicated(demot[ , "id"]), "id"]))
duplrest <- NULL
for(i in 1:length(duplid)) {
  tt <- NULL
  for(k in 1:length(dn2)) {
    myx <- sort(row.names(demot)[complete.cases(demot[ , c("id", "dataset")]) &
                                   demot[ , "id"] == duplid[i] & demot[ , "dataset"] == dn2[k]])
    if(length(myx) > 0) { tt <- c(tt, myx) }
  }
  duplrest <- c(duplrest, list(tt))
}
names(duplrest) <- duplid
duplres <- c(duplres, duplrest)
## Full set duplicated patients
duPL <- sort(unlist(lapply(duplres, function(x) { return(x[-1]) } )))

dn <- c("transbig", "unt", "upp", "mainz", "nki")
dn.platform <- c("affy", "affy", "affy", "affy", "agilent")
res <- ddemo.all <- ddemo.coln <- NULL
for(i in 1:length(dn)) {
  ## load dataset
  dd <- get(data(list=dn[i]))
  #Remove duplicates identified first
  message("obtained dataset!")
  #Extract expression set, pData, fData for each dataset
  ddata <- t(exprs(dd))
  ddemo <- phenoData(dd)@data
  if(length(intersect(rownames(ddata),duPL))>0)
  {
    ddata<-ddata[-which(rownames(ddata) %in% duPL),]
    ddemo<-ddemo[-which(rownames(ddemo) %in% duPL),]
  }
  dannot <- featureData(dd)@data
  #
  #
  #
  #
  # OBSOLETE FUNCTION CALL - OLDER VERSIONS OF GENEFU
  # SubtypePredictions<-subtype.cluster.predict(sbt.model=scmod2.robust,data=ddata,
  # annot=dannot,do.mapping=TRUE,verbose=TRUE)
# CURRENT FUNCTION CALL - NEWEST VERSION OF GENEFU
SubtypePredictions<-molecular.subtyping(sbt.model = "scmod2",data = ddata,
                                        annot = dannot,do.mapping = TRUE)
#Get sample counts pertaining to each subtype
table(SubtypePredictions$subtype)
#Select samples pertaining to Basal Subtype
Basals<-names(which(SubtypePredictions$subtype == "ER-/HER2-"))
#Select samples pertaining to HER2 Subtype
HER2s<-names(which(SubtypePredictions$subtype == "HER2+"))
#Select samples pertaining to Luminal Subtypes
LuminalB<-names(which(SubtypePredictions$subtype == "ER+/HER2- High Prolif"))
LuminalA<-names(which(SubtypePredictions$subtype == "ER+/HER2- Low Prolif"))
#ASSIGN SUBTYPES TO EVERY SAMPLE, ADD TO THE EXISTING PHENODATA
ddemo$SCMOD2<-SubtypePredictions$subtype
ddemo[LuminalB,]$SCMOD2<-"LumB"
ddemo[LuminalA,]$SCMOD2<-"LumA"
ddemo[Basals,]$SCMOD2<-"Basal"
ddemo[HER2s,]$SCMOD2<-"Her2"

PAM50Preds<-molecular.subtyping(sbt.model = "pam50",data=ddata,
                                annot=dannot,do.mapping=TRUE)
table(PAM50Preds$subtype)
ddemo$PAM50<-PAM50Preds$subtype
LumA<-names(PAM50Preds$subtype)[which(PAM50Preds$subtype == "LumA")]
LumB<-names(PAM50Preds$subtype)[which(PAM50Preds$subtype == "LumB")]
ddemo[LumA,]$PAM50<-"LumA"
ddemo[LumB,]$PAM50<-"LumB"

ddemo.all <- rbind(ddemo, ddemo.all)
}

table(ddemo.all$PAM50)
Normals<-rownames(ddemo.all[which(ddemo.all$PAM50 == "Normal"),])
# Obtain the subtype prediction counts for SCMOD2
table(ddemo.all$SCMOD2)

ddemo.all$PAM50<-as.character(ddemo.all$PAM50)
confusionMatrix(factor(ddemo.all[-which(rownames(ddemo.all) %in% Normals),]$SCMOD2),
                factor(ddemo.all[-which(rownames(ddemo.all) %in% Normals),]$PAM50))

# https://bioconductor.org/packages/release/bioc/vignettes/genefu/inst/doc/genefu.pdf

citation("genefu")
help("citation")
citation("limma")