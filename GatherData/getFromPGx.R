library(Biobase)
library(PharmacoGx)

## to download the datasets from server
#GDSC <- downloadPSet("GDSC")
#CCLE <- downloadPSet("CCLE")

## small dataset for testing & debugging
# data("GDSCsmall")
# data("CCLEsmall")
# GDSC <- GDSCsmall
# CCLE <- CCLEsmall

## use local files
load("PSets/GDSC.RData")
# load("PSets/CCLE.RData")

extractDrugData <- function (drugID, dataPset) {
  ## get AUC for the drug
  label.auc <- summarizeSensitivityProfiles(dataPset,
                                           drugs=c(drugID),
                                           sensitivity.measure="auc_published",
                                           summary.stat="median")
  label.auc <- label.auc[1, !is.na(label.auc)]
  
  ## get IC50 for the drug
  label.IC50 <- summarizeSensitivityProfiles(dataPset,
                                           drugs=c(drugID),
                                           sensitivity.measure="ic50_published",
                                           summary.stat="median")
  label.IC50 <- label.IC50[1, !is.na(label.IC50)]
  
  ## get expression data of cell lines with sensitivity measurement for the drug
  data.expression <- summarizeMolecularProfiles(dataPset,
                                               # cell.lines=cc,
                                               mDataType="rna",
                                               verbose=FALSE)
  ## rename genes to EnsemblGeneId
  gInfo <- featureInfo(dataPset, 'rna')
  rna.all <- exprs(data.expression)
  rownames(rna.all) <- gInfo[rownames(rna.all),]$EnsemblGeneId
  
  ## get list of cell lines
  cl.all <- colnames(rna.all)
  cl.auc <- intersect(names(label.auc), cl.all)
  cl.IC50 <- intersect(names(label.IC50), cl.all)
  
  ## filter for cell lines with both expression and measurement data
  rna.auc <- rna.all[, cl.auc]
  rna.auc <- rna.auc[, colSums(is.na(rna.auc))==0]
  cl.auc <- colnames(rna.auc)
  label.auc <- label.auc[cl.auc]
  
  rna.IC50 <- rna.all[, cl.IC50]
  rna.IC50 <- rna.IC50[, colSums(is.na(rna.IC50))==0]
  cl.IC50 <- colnames(rna.IC50)
  label.IC50 <- label.IC50[cl.IC50]
  
  stopifnot(all(colnames(rna.auc) == colnames(label.auc)))
  stopifnot(all(colnames(rna.IC50) == colnames(label.IC50)))
  
  return(list('rna.auc'=rna.auc, 'rna.IC50'=rna.IC50, 'auc.cont'=label.auc, 'IC50.cont'=label.IC50))
}

### Bortezomib
sampleinfo.gdsc <- cellInfo(GDSC)
temp <- extractDrugData('Bortezomib', GDSC)
bortezomib <- list('gdsc_AUC'=temp$rna.auc, 'gdsc_IC50'=temp$rna.IC50)
bortezomib.labels <- list('AUC.cont'=temp$auc.cont, 'IC50.cont'=temp$IC50.cont)

## binarize the sensitivity measure, the new (as of jan2016) way
source("callSensitivity.R")
bortezomib.labels$AUC <- sapply(callSensitivity(t(data.frame(bortezomib.labels$AUC.cont)))[1, ], as.logical)
bortezomib.labels$IC50 <- sapply(callSensitivity(t(data.frame(bortezomib.labels$IC50.cont)))[1, ], as.logical)

# indices of cell lines to the sampleinfo table
bortezomib.labels$AUC_ind <- which(rownames(sampleinfo.gdsc) %in% colnames(bortezomib$rna.auc))
bortezomib.labels$IC50_ind <- which(rownames(sampleinfo.gdsc) %in% colnames(bortezomib$rna.IC50))

## save the data set to RData
# save(bortezomib, bortezomib.labels, sampleinfo.gdsc, file='bortezomib_new.RData')

# debug out
dim(bortezomib$gdsc_AUC)
length(bortezomib.labels$AUC.cont)
length(bortezomib.labels$AUC)
dim(bortezomib$gdsc_IC50)
length(bortezomib.labels$IC50.cont)
length(bortezomib.labels$IC50)

View(bortezomib$gdsc_AUC)
View(sampleinfo.gdsc)
View(bortezomib.labels$AUC)



