library(Biobase)
library(PharmacoGx)

## to download the datasets from server
#GDSC <- downloadPSet("GDSC")
#CCLE <- downloadPSet("CCLE")

## small dataset for testing & debugging
data("GDSCsmall")
data("CCLEsmall")
GDSC <- GDSCsmall
CCLE <- CCLEsmall

## use local files
#load("PSets/GDSC.RData")
# load("PSets/CCLE.RData")

extractDrugData <- function (drugID, dataPset) {
  ## get AUC for the drug
  drug.auc <- summarizeSensitivityProfiles(dataPset,
                                           drugs=c(drugID),
                                           sensitivity.measure="auc_recomputed",
                                           summary.stat="median")
  drug.auc <- drug.auc[1, !is.na(drug.auc)]
  
  ## get IC50 for the drug
  drug.IC50 <- summarizeSensitivityProfiles(dataPset,
                                           drugs=c(drugID),
                                           sensitivity.measure="ic50_recomputed",
                                           summary.stat="median")
  drug.IC50 <- drug.IC50[1, !is.na(drug.IC50)]
  
  ## get list of cell lines
  cc.auc <- names(drug.auc)
  cc.IC50 <- names(drug.IC50)
  
  ## get expression data of cell lines with sensitivity measurement for the drug
  data.expression <- summarizeMolecularProfiles(dataPset,
                                               # cell.lines=cc,
                                               mDataType="rna",
                                               verbose=FALSE)
  ## rename genes to EnsemblGeneId
  gInfo <- featureInfo(dataPset, 'rna')
  rna <- exprs(data.expression)
  rownames(rna) <- gInfo[rownames(rna),]$EnsemblGeneId
  
  ## TODO: filter for cell lines with both expression and measurement data
  auc.rna <- rna
  IC50.rna <- rna
  
  return(list('auc.rna'=auc.rna, 'auc.cont'=drug.auc, 'IC50.rna'=IC50.rna, 'IC50.cont'=drug.IC50))
}

### Bortezomib
bortezomib <- extractDrugData('Bortezomib', GDSC)

## binarize the sensitivity measure
source("callSensitivity.R")
bortezomib$auc.bin <- callSensitivity(t(data.frame(bortezomib$auc.cont)))


dim(bortezomib$auc.rna)
dim(bortezomib$auc.cont)
dim(bortezomib$auc.bin)
View(bortezomib$auc.rna)
View(bortezomib$auc.cont)
View(bortezomib$auc.bin)





