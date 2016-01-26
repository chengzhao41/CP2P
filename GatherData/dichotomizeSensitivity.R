dichotomizeSensitivity <- function(labels) {
  ## binarize the sensitivity measure, the new (as of jan2016) way
  source("GatherData/callSensitivity.R")
  labels$AUCnew <- sapply(callSensitivity(t(data.frame(labels$AUC.cont)))[1, ], as.logical)
  labels$IC50new <- !sapply(callSensitivity(t(data.frame(labels$IC50.cont)))[1, ], as.logical)
  labels$slopenew <- sapply(callSensitivity(t(data.frame(labels$slope.cont)))[1, ], as.logical)
  
  ## binarize the sensitivity measure, using the old "water fall" method
  source('Common/drug_cut/callingWaterfall.R')
  source('Common/drug_cut/distancePointLine.R')
  source('Common/drug_cut/distancePointSegment.R')
  labels$AUC <- callingWaterfall(t(data.frame(labels$AUC.cont)), type="AUC")
  labels$AUC <- labels$AUC != "resistant"
  names(labels$AUC) <- names(labels$AUC.cont)
  
  labels$IC50 <- callingWaterfall(t(data.frame(labels$IC50.cont)), type="IC50")
  labels$IC50 <- labels$IC50 != "resistant"
  names(labels$IC50) <- names(labels$IC50.cont)
  
  labels$slope <- callingWaterfall(t(data.frame(labels$slope.cont)), type="AUC") # there is no "type" for slope
  labels$slope <- labels$slope != "resistant"
  names(labels$slope) <- names(labels$slope.cont)
  
#   slope_cutoff <- 0.27
#   labels$slope <- labels$slope.cont
#   labels$slope[which(labels$slope.cont < slope_cutoff)] <- 0
#   labels$slope[which(labels$slope.cont >= slope_cutoff)] <- 1
#   labels$slope <- as.logical(labels$slope)
  
  return(labels)
}