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

### Bortezomib
temp <- extractDrugData('Bortezomib', GDSC)
bortezomib <- list('gdsc_AUC'=temp$rna.auc, 'gdsc_IC50'=temp$rna.IC50, 'gdsc_slope'=temp$rna.slope)
bortezomib.labels <- list('AUC.cont'=temp$auc.cont, 'IC50.cont'=temp$IC50.cont, 'slope.cont'=temp$slope.cont)

# dichotomize
bortezomib.labels <- dichotomizeSensitivity(bortezomib.labels)

# indices of cell lines to the sampleinfo table
sampleinfo.gdsc <- cellInfo(GDSC)
bortezomib.labels$AUC_ind <- which(rownames(sampleinfo.gdsc) %in% colnames(bortezomib$gdsc_AUC))
bortezomib.labels$IC50_ind <- which(rownames(sampleinfo.gdsc) %in% colnames(bortezomib$gdsc_IC50))
bortezomib.labels$slope_ind <- which(rownames(sampleinfo.gdsc) %in% colnames(bortezomib$gdsc_slope))

## save the data set to RData
# transpose
bortezomib.labels$AUC_ind <- t(bortezomib.labels$AUC_ind)
bortezomib.labels$IC50_ind <- t(bortezomib.labels$IC50_ind)
bortezomib.labels$slope_ind <- t(bortezomib.labels$slope_ind)
save(bortezomib, bortezomib.labels, sampleinfo.gdsc, file='Bortezomib/WS/bortezomib_gdsc.RData')

# debug out
table(bortezomib.labels$AUC)
table(bortezomib.labels$IC50)
table(bortezomib.labels$slope)
table(bortezomib.labels$AUCnew)
table(bortezomib.labels$IC50new)
table(bortezomib.labels$slopenew)
sum(bortezomib.labels$AUC == bortezomib.labels$AUCnew) / length(bortezomib.labels$AUC)
sum(bortezomib.labels$IC50 == bortezomib.labels$IC50new) / length(bortezomib.labels$IC50)
sum(bortezomib.labels$slope == bortezomib.labels$slopenew) / length(bortezomib.labels$slope)

dim(bortezomib$gdsc_AUC)
length(bortezomib.labels$AUC.cont)
length(bortezomib.labels$AUC)
dim(bortezomib$gdsc_IC50)
length(bortezomib.labels$IC50.cont)
length(bortezomib.labels$IC50)
dim(bortezomib$gdsc_slope)
length(bortezomib.labels$slope.cont)
length(bortezomib.labels$slope)

#View(bortezomib$gdsc_AUC)
#View(sampleinfo.gdsc)
View(bortezomib.labels$AUC)
View(bortezomib.labels$slope)


