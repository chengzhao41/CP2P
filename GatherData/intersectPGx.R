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
#load("PSets/CCLE.RData")

## intersect GDSC & CCLE (modified from PharamcoGx vignette)
intersectGDSCnCCLE <- function (GDSC, CCLE) {
  commonGenes <- intersect(rownames(featureInfo(GDSC, "rna")),
                           rownames(featureInfo(CCLE, 'rna')))
  
  common <- intersectPSet(list('CCLE'=CCLE,
                               'GDSC'=GDSC),
                          intersectOn=c("cell.lines", "drugs"))
  
  GDSCexpression <- summarizeMolecularProfiles(common$GDSC,
                                               cellNames(common$GDSC),
                                               mDataType="rna",
                                               features=commonGenes,
                                               verbose=FALSE)
  
  CCLEexpression <- summarizeMolecularProfiles(common$CCLE,
                                               cellNames(common$CCLE),
                                               mDataType="rna",
                                               features=commonGenes,
                                               verbose=FALSE)
  
  GDSC.auc <- summarizeSensitivityProfiles(
    common$GDSC,
    sensitivity.measure='auc_published',
    summary.stat="median")
  CCLE.auc <- summarizeSensitivityProfiles(
    common$CCLE,
    sensitivity.measure='auc_published',
    summary.stat="median")
  
  ## names of common genes and cell lines
  gg <- featureInfo(common[[1]], 'rna')
  cc <- cellNames(common[[1]])
  
  ## rename genes to EntrezGeneId
  gInfo <- featureInfo(GDSC, 'rna')
  gdsc.rna <- exprs(GDSCexpression)
  rownames(gdsc.rna) <- gInfo[rownames(gdsc.rna),]$EntrezGeneId
  
  gInfo <- featureInfo(CCLE, 'rna')
  ccle.rna <- exprs(CCLEexpression)
  rownames(ccle.rna) <- gInfo[rownames(ccle.rna),]$EntrezGeneId
  
  return(list(gdsc.rna, ccle.rna, GDSC.auc, CCLE.auc))
}


## example usage: intersect GDSC, CCLE and then binarize the sensitivity measure
source("callSensitivity.R")
its <- intersectGDSCnCCLE(GDSC, CCLE)
a.bin <- callSensitivity(its[[3]])
b.bin <- callSensitivity(its[[4]])
