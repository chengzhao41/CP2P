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

### Docetaxel
temp <- extractDrugData('Docetaxel', GDSC)
docetaxel <- list('gdsc_AUC'=temp$rna.auc, 'gdsc_IC50'=temp$rna.IC50, 'gdsc_slope'=temp$rna.slope)
docetaxel.labels <- list('AUC.cont'=temp$auc.cont, 'IC50.cont'=temp$IC50.cont, 'slope.cont'=temp$slope.cont)

# dichotomize
docetaxel.labels <- dichotomizeSensitivity(docetaxel.labels)

# indices of cell lines to the sampleinfo table
sampleinfo.gdsc <- cellInfo(GDSC)
docetaxel.labels$AUC_ind <- which(rownames(sampleinfo.gdsc) %in% colnames(docetaxel$gdsc_AUC))
docetaxel.labels$IC50_ind <- which(rownames(sampleinfo.gdsc) %in% colnames(docetaxel$gdsc_IC50))
docetaxel.labels$slope_ind <- which(rownames(sampleinfo.gdsc) %in% colnames(docetaxel$gdsc_slope))

## save the data set to RData
# transpose
docetaxel$gdsc_AUC <- t(docetaxel$gdsc_AUC)
docetaxel$gdsc_IC50 <- t(docetaxel$gdsc_IC50)
docetaxel$gdsc_slope <- t(docetaxel$gdsc_slope)
save(docetaxel, docetaxel.labels, sampleinfo.gdsc, file='Docetaxel/WS/docetaxel_gdsc.RData')

# debug out
table(docetaxel.labels$AUC)
table(docetaxel.labels$IC50)
table(docetaxel.labels$slope)
table(docetaxel.labels$AUCnew)
table(docetaxel.labels$IC50new)
table(docetaxel.labels$slopenew)
sum(docetaxel.labels$AUC == docetaxel.labels$AUCnew) / length(docetaxel.labels$AUC)
sum(docetaxel.labels$IC50 == docetaxel.labels$IC50new) / length(docetaxel.labels$IC50)
sum(docetaxel.labels$slope == docetaxel.labels$slopenew) / length(docetaxel.labels$slope)

dim(docetaxel$gdsc_AUC)
length(docetaxel.labels$AUC.cont)
length(docetaxel.labels$AUC)
dim(docetaxel$gdsc_IC50)
length(docetaxel.labels$IC50.cont)
length(docetaxel.labels$IC50)
dim(docetaxel$gdsc_slope)
length(docetaxel.labels$slope.cont)
length(docetaxel.labels$slope)

#View(docetaxel$gdsc_AUC)
#View(sampleinfo.gdsc)
View(docetaxel.labels$AUC)
View(docetaxel.labels$slope)


