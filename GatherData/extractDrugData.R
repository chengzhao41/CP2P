require(PharmacoGx)

extractDrugData <- function (drugID, dataPset, mDataType='rna') {
  ## get AUC for the drug
  label.auc <- summarizeSensitivityProfiles(dataPset,
                                            drugs=c(drugID),
                                            sensitivity.measure="auc_recomputed",
                                            summary.stat="median")
  label.auc <- label.auc[1, !is.na(label.auc)]
  
  ## get IC50 for the drug
  label.IC50 <- summarizeSensitivityProfiles(dataPset,
                                             drugs=c(drugID),
                                             sensitivity.measure="ic50_recomputed",
                                             summary.stat="median")
  label.IC50 <- label.IC50[1, !is.na(label.IC50)]
  
  ## get Slope for the drug
  label.slope <- summarizeSensitivityProfiles(dataPset,
                                              drugs=c(drugID),
                                              sensitivity.measure="slope_recomputed",
                                              summary.stat="median")
  label.slope <- label.slope[1, !is.na(label.slope)]
  
  ## get expression data of cell lines with sensitivity measurement for the drug
  data.expression <- summarizeMolecularProfiles(dataPset,
                                                # cell.lines=cc,
                                                mDataType=mDataType,
                                                verbose=FALSE)
  ## rename genes to EnsemblGeneId
  gInfo <- featureInfo(dataPset, mDataType)
  rna.all <- exprs(data.expression)
  rownames(rna.all) <- gInfo[rownames(rna.all),]$EnsemblGeneId
  
  ## get list of cell lines
  cl.all <- colnames(rna.all)
  cl.auc <- intersect(names(label.auc), cl.all)
  cl.IC50 <- intersect(names(label.IC50), cl.all)
  cl.slope <- intersect(names(label.slope), cl.all)
  
  ## filter for cell lines with both expression and measurement data
  rna.auc <- rna.all[, cl.auc]
  rna.auc <- rna.auc[, colSums(is.na(rna.auc))==0]
  cl.auc <- colnames(rna.auc)
  label.auc <- label.auc[cl.auc]
  
  rna.IC50 <- rna.all[, cl.IC50]
  rna.IC50 <- rna.IC50[, colSums(is.na(rna.IC50))==0]
  cl.IC50 <- colnames(rna.IC50)
  label.IC50 <- label.IC50[cl.IC50]
  
  rna.slope <- rna.all[, cl.slope]
  rna.slope <- rna.slope[, colSums(is.na(rna.slope))==0]
  cl.slope <- colnames(rna.slope)
  label.slope <- label.slope[cl.slope]
  
  stopifnot(all(colnames(rna.auc) == colnames(label.auc)))
  stopifnot(all(colnames(rna.IC50) == colnames(label.IC50)))
  stopifnot(all(colnames(rna.slope) == colnames(label.slope)))
  
  return(list('rna.auc'=rna.auc, 'rna.IC50'=rna.IC50, 'rna.slope'=rna.slope, 
              'auc.cont'=label.auc, 'IC50.cont'=label.IC50, 'slope.cont'=label.slope))
}