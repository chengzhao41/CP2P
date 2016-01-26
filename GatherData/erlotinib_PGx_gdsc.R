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
load("GatherData/PSets/GDSC.RData")
# load("GatherData/PSets/CCLE.RData")

## cap all the ic50 at maximum concentartion
tt <- which(as.numeric(GDSC@sensitivity$profiles$ic50_recomputed) > as.numeric(GDSC@sensitivity$info$max.Dose.uM))
GDSC@sensitivity$profiles[tt, "ic50_recomputed"] <- GDSC@sensitivity$info[tt, "max.Dose.uM"]

## compute "slope_recomputed" measure
if (! "slope_recomputed" %in% names(GDSC@sensitivity$profiles)) {
  GDSC <- computeSlopeMeasure(GDSC)
  #length(GDSC@sensitivity$profiles$slope_recomputed)
  #length(GDSC@sensitivity$profiles$auc_recomputed)
}

### Erlotinib
temp <- extractDrugData('Erlotinib', GDSC)
erlotinib <- list('gdsc_AUC'=temp$rna.auc, 'gdsc_IC50'=temp$rna.IC50, 'gdsc_slope'=temp$rna.slope)
erlotinib.labels <- list('AUC.cont'=temp$auc.cont, 'IC50.cont'=temp$IC50.cont, 'slope.cont'=temp$slope.cont)

# dichotomize
erlotinib.labels <- dichotomizeSensitivity(erlotinib.labels)

# indices of cell lines to the sampleinfo table
sampleinfo.gdsc <- cellInfo(GDSC)
erlotinib.labels$AUC_ind <- which(rownames(sampleinfo.gdsc) %in% colnames(erlotinib$gdsc_AUC))
erlotinib.labels$IC50_ind <- which(rownames(sampleinfo.gdsc) %in% colnames(erlotinib$gdsc_IC50))
erlotinib.labels$slope_ind <- which(rownames(sampleinfo.gdsc) %in% colnames(erlotinib$gdsc_slope))

## save the data set to RData
# transpose
erlotinib$gdsc_AUC <- t(erlotinib$gdsc_AUC)
erlotinib$gdsc_IC50 <- t(erlotinib$gdsc_IC50)
erlotinib$gdsc_slope <- t(erlotinib$gdsc_slope)
save(erlotinib, erlotinib.labels, sampleinfo.gdsc, file='Erlotinib/WS/erlotinib_gdsc.RData')

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

dim(erlotinib$gdsc_AUC)
length(erlotinib.labels$AUC.cont)
length(erlotinib.labels$AUC)
dim(erlotinib$gdsc_IC50)
length(erlotinib.labels$IC50.cont)
length(erlotinib.labels$IC50)
dim(erlotinib$gdsc_slope)
length(erlotinib.labels$slope.cont)
length(erlotinib.labels$slope)

#View(erlotinib$gdsc_AUC)
#View(sampleinfo.gdsc)
View(erlotinib.labels$AUC)
View(erlotinib.labels$slope)


