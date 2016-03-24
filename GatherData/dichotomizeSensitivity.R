dichotomizeSensitivityWaterfall <- function(labels) {
  ## binarize the sensitivity measure, using the old "waterfall" method
  source('Common/drug_cut/callingWaterfall.R')
  source('Common/drug_cut/distancePointLine.R')
  source('Common/drug_cut/distancePointSegment.R')
  lls <- NULL
  
  if ("AUC.cont" %in% names(labels)) {
    lls$AUC <- callingWaterfall(t(data.frame(labels$AUC.cont)), type="AUC")
    lls$AUC <- lls$AUC != "resistant"
    names(lls$AUC) <- names(labels$AUC.cont)
  }
  if ("IC50.cont" %in% names(labels)) {
    lls$IC50 <- callingWaterfall(t(data.frame(labels$IC50.cont)), type="IC50")
    lls$IC50 <- lls$IC50 != "resistant"
    names(lls$IC50) <- names(labels$IC50.cont)
  }
  if ("slope.cont" %in% names(labels)) {
    lls$slope <- callingWaterfall(t(data.frame(labels$slope.cont)), type="AUC") # there is no "type" for slope
    lls$slope <- lls$slope != "resistant"
    names(lls$slope) <- names(labels$slope.cont)
    
    ## hard cutoff for slope used by Benjamin's lab before
#     slope_cutoff <- 0.29
#     lls$slope <- labels$slope.cont
#     lls$slope[which(labels$slope.cont < slope_cutoff)] <- 0
#     lls$slope[which(labels$slope.cont >= slope_cutoff)] <- 1
#     lls$slope <- as.logical(lls$slope)
  }
  return(lls)
}

dichotomizeSensitivityNew <- function(labels, pgxPset, drugName, use_all_data=FALSE) {
  ## binarize the sensitivity measure, the new (as of jan2016) way
  ## use_all_data: threshold computed on whole data set (like CCLE or GDSC) with all drugs together (can get better stability)
  source("GatherData/callSensitivity.R")
  
  drugs <- drugNames(pgxPset)
  if (use_all_data == FALSE) {
    drugs <- c(drugName)
  }
  
  lls <- NULL
  
  ## binarize "auc_recomputed" measure
  if ("AUC.cont" %in% names(labels)) {
    all.auc <- summarizeSensitivityProfiles(pgxPset, drugs=drugs, sensitivity.measure="auc_recomputed")
    temp <- callSensitivity(all.auc)
    print(paste(drugName, 'AUC callSensitvity cutoff:', temp$cutoff))

    auc.bin <- temp$sensitivity.call[drugName, names(labels$AUC.cont)]
    lls$AUC.new <- sapply(auc.bin, as.logical)
  }
  ## binarize "ic50_recomputed" measure
  if ("IC50.cont" %in% names(labels)) {
    all.ic50 <- summarizeSensitivityProfiles(pgxPset, drugs=drugs, sensitivity.measure="ic50_recomputed")
    temp <- callSensitivity(all.ic50)
    print(paste(drugName, 'IC50 callSensitvity cutoff:', temp$cutoff))
    
    ic50.bin <- temp$sensitivity.call[drugName, names(labels$IC50.cont)]
    lls$IC50.new <- !sapply(ic50.bin, as.logical) #negate as higher IC50 means less sensitivity
  }
  ## binarize "slope_recomputed" measure
  if ("slope.cont" %in% names(labels)) {
    all.slope <- summarizeSensitivityProfiles(pgxPset, drugs=drugs, sensitivity.measure="slope_recomputed")
    temp <- callSensitivity(all.slope)
    print(paste(drugName, 'Slope callSensitvity cutoff:', temp$cutoff))
    
    slope.bin <- temp$sensitivity.call[drugName, names(labels$slope.cont)]
    lls$slope.new <- sapply(slope.bin, as.logical)
  }
  return(lls)
}