library(Biobase)
library(PharmacoGx)

setwd("~/Dropbox/CP2P")

source("GatherData/extractDrugData.R")
source("GatherData/computeSlopeMeasure.R")
source("GatherData/dichotomizeSensitivity.R")

## use local Pset files
load("GatherData/PSets/GSE32036.RData")

### GSE32036
GSE32036@sensitivity$profiles$ic50_published <- as.numeric(GSE32036@sensitivity$profiles$ic50_published) 
temp <- extractDrugData_forGSE32036('Erlotinib', GSE32036)
battle <- list('IC50'=temp$rna.IC50)
battle.labels <- list('IC50.cont'=temp$IC50.cont)

## Dichotomize
# old-style waterfall
battle.labels <- c(dichotomizeSensitivityWaterfall(battle.labels), battle.labels)

## save the data set to RData
# transpose
battle$IC50 <- t(battle$IC50)
save(battle, battle.labels, file='Erlotinib/WS/erlotinib_battle.RData')
