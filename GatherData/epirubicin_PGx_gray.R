library(Biobase)
library(PharmacoGx)

setwd("~/Dropbox/CP2P")

source("GatherData/extractDrugData.R")
source("GatherData/computeSlopeMeasure.R")
source("GatherData/dichotomizeSensitivity.R")

## use local Pset files
load("GatherData/PSets/GRAY.RData")

## cap all the ic50 at maximum concentartion
tt <- which(as.numeric(GRAY@sensitivity$profiles$ic50_recomputed) > as.numeric(GRAY@sensitivity$info$max.Dose.uM))
GRAY@sensitivity$profiles[tt, "ic50_recomputed"] <- GRAY@sensitivity$info[tt, "max.Dose.uM"]

## compute "slope_recomputed" measure
if (! "slope_recomputed" %in% names(GRAY@sensitivity$profiles)) {
  GRAY <- computeSlopeMeasure(GRAY)
  #length(GRAY@sensitivity$profiles$slope_recomputed)
  #length(GRAY@sensitivity$profiles$auc_recomputed)
}

### Epirubicin
temp <- extractDrugData('Epirubicin', GRAY, mDataType='rnaseq')
epirubicin <- list('gray_AUC'=temp$rna.auc, 'gray_IC50'=temp$rna.IC50, 'gray_slope'=temp$rna.slope)
epirubicin.labels <- list('AUC.cont'=temp$auc.cont, 'IC50.cont'=temp$IC50.cont, 'slope.cont'=temp$slope.cont)

# dichotomize
epirubicin.labels <- dichotomizeSensitivity(epirubicin.labels)

# indices of cell lines to the sampleinfo table
sampleinfo.gray <- cellInfo(GRAY)
epirubicin.labels$AUC_ind <- which(rownames(sampleinfo.gray) %in% colnames(epirubicin$gray_AUC))
epirubicin.labels$IC50_ind <- which(rownames(sampleinfo.gray) %in% colnames(epirubicin$gray_IC50))
epirubicin.labels$slope_ind <- which(rownames(sampleinfo.gray) %in% colnames(epirubicin$gray_slope))

## save the data set to RData
# transpose
epirubicin$gray_AUC <- t(epirubicin$gray_AUC)
epirubicin$gray_IC50 <- t(epirubicin$gray_IC50)
epirubicin$gray_slope <- t(epirubicin$gray_slope)
save(epirubicin, epirubicin.labels, sampleinfo.gray, file='Epirubicin/WS/epirubicin_gray.RData')

# debug out
table(epirubicin.labels$AUC)
table(epirubicin.labels$IC50)
table(epirubicin.labels$slope)
table(epirubicin.labels$AUCnew)
table(epirubicin.labels$IC50new)
table(epirubicin.labels$slopenew)
sum(epirubicin.labels$AUC == epirubicin.labels$AUCnew) / length(epirubicin.labels$AUC)
sum(epirubicin.labels$IC50 == epirubicin.labels$IC50new) / length(epirubicin.labels$IC50)
sum(epirubicin.labels$slope == epirubicin.labels$slopenew) / length(epirubicin.labels$slope)

dim(epirubicin$gray_AUC)
length(epirubicin.labels$AUC.cont)
length(epirubicin.labels$AUC)
dim(epirubicin$gray_IC50)
length(epirubicin.labels$IC50.cont)
length(epirubicin.labels$IC50)
dim(epirubicin$gray_slope)
length(epirubicin.labels$slope.cont)
length(epirubicin.labels$slope)

#View(epirubicin$gray_AUC)
#View(sampleinfo.gray)
View(epirubicin.labels$AUC)
View(epirubicin.labels$slope)


