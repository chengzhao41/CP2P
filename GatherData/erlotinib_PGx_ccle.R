library(Biobase)
library(PharmacoGx)

setwd("~/Dropbox/CP2P")

source("GatherData/extractDrugData.R")
source("GatherData/computeSlopeMeasure.R")
source("GatherData/dichotomizeSensitivity.R")


## to download the datasets from server
# GDSC <- downloadPSet("GDSC")
# CCLE <- downloadPSet("CCLE")

## small dataset for testing & debugging
# data("GDSCsmall"); GDSC <- GDSCsmall
# data("CCLEsmall"); CCLE <- CCLEsmall

## use local Pset files
# load("GatherData/PSets/GDSC.RData")
load("GatherData/PSets/CCLE.RData")

## cap all the ic50 at maximum concentartion
tt <- which(as.numeric(CCLE@sensitivity$profiles$ic50_recomputed) > as.numeric(CCLE@sensitivity$info$Dose8.uM))
CCLE@sensitivity$profiles[tt, "ic50_recomputed"] <- CCLE@sensitivity$info[tt, "Dose8.uM"]

## compute "slope_recomputed" measure
if (! "slope_recomputed" %in% names(CCLE@sensitivity$profiles)) {
  CCLE <- computeSlopeMeasure(CCLE)
  #length(CCLE@sensitivity$profiles$slope_recomputed)
  #length(CCLE@sensitivity$profiles$auc_recomputed)
}

### Erlotinib
temp <- extractDrugData('Erlotinib', CCLE)
erlotinib <- list('ccle_AUC'=temp$rna.auc, 'ccle_IC50'=temp$rna.IC50, 'ccle_slope'=temp$rna.slope)
erlotinib.labels <- list('AUC.cont'=temp$auc.cont, 'IC50.cont'=temp$IC50.cont, 'slope.cont'=temp$slope.cont)

# dichotomize
erlotinib.labels <- dichotomizeSensitivity(erlotinib.labels)

# indices of cell lines to the sampleinfo table
sampleinfo.ccle <- cellInfo(CCLE)
erlotinib.labels$AUC_ind <- which(rownames(sampleinfo.ccle) %in% colnames(erlotinib$ccle_AUC))
erlotinib.labels$IC50_ind <- which(rownames(sampleinfo.ccle) %in% colnames(erlotinib$ccle_IC50))
erlotinib.labels$slope_ind <- which(rownames(sampleinfo.ccle) %in% colnames(erlotinib$ccle_slope))

## save the data set to RData
# transpose
erlotinib$ccle_AUC <- t(erlotinib$ccle_AUC)
erlotinib$ccle_IC50 <- t(erlotinib$ccle_IC50)
erlotinib$ccle_slope <- t(erlotinib$ccle_slope)
save(erlotinib, erlotinib.labels, sampleinfo.ccle, file='Erlotinib/WS/erlotinib_ccle.RData')

# debug out
table(erlotinib.labels$AUC)
table(erlotinib.labels$IC50)
table(erlotinib.labels$slope)
table(erlotinib.labels$AUCnew)
table(erlotinib.labels$IC50new)
table(erlotinib.labels$slopenew)
sum(erlotinib.labels$AUC == erlotinib.labels$AUCnew) / length(erlotinib.labels$AUC)
sum(erlotinib.labels$IC50 == erlotinib.labels$IC50new) / length(erlotinib.labels$IC50)
sum(erlotinib.labels$slope == erlotinib.labels$slopenew) / length(erlotinib.labels$slope)

dim(erlotinib$ccle_AUC)
length(erlotinib.labels$AUC.cont)
length(erlotinib.labels$AUC)
dim(erlotinib$ccle_IC50)
length(erlotinib.labels$IC50.cont)
length(erlotinib.labels$IC50)
dim(erlotinib$ccle_slope)
length(erlotinib.labels$slope.cont)
length(erlotinib.labels$slope)

#View(erlotinib$ccle_AUC)
#View(sampleinfo.ccle)
View(erlotinib.labels$AUC)
View(erlotinib.labels$slope)


