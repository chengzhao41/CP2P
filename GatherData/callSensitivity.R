bimod <- function (data, model=c("E", "V"), do.scale=TRUE, verbose=FALSE) 
{
  require(mclust)
  mystatus <- mystatus.proba <- rep(NA, length(data))
  #names(mystatus) <- names(mystatus.proba) <- dimnames(data)[[1]]
  res <- matrix(NA, nrow=3, ncol=2, dimnames=list(c("mean", "variance", "proportion"), paste("cluster", 1:2, sep=".")))
  
  #Only 2 Gaussians
  rr2 <- mclust::Mclust(data, modelNames=model, G=2)
  res[1, ] <- rr2$parameters$mean
  res[2, ] <- rr2$parameters$variance$sigmasq
  res[3, ] <- rr2$parameters$pro
  
  ## bimodality index (BI)
  smd <- abs(res[1, 2] - res[1, 1]) / sqrt((res[2, 2] + res[2, 1]) / 2)
  bi <- sqrt(res[3, 2] * (1 - res[3, 2])) * smd
  
  #classification
  mystatus <- as.numeric(rr2$classification == 2)
  mystatus.proba <- rr2$z[ , 2, drop=TRUE]
  return(list("status"=mystatus, "status1.proba"=mystatus.proba, "gaussians"=res,  "BI"=bi))
}

## binarize continuous sensitivity measure
callSensitivity <- function (cs) {
  xx <- as.vector(cs)
  xx <- xx[which(!is.na(xx))]
  data.bimod <- bimod(xx, model=c("E", "V"), do.scale=TRUE, verbose=FALSE) 
  data.cutoff <- qnorm(.99, mean = data.bimod$gaussians["mean","cluster.1"], sd = sqrt(data.bimod$gaussians["variance","cluster.1"]))
  
  sensitivity.call <- matrix(NA, ncol=ncol(cs), nrow=nrow(cs), dimnames=list(rownames(cs), colnames((cs))))
  sensitivity.call[which(cs < data.cutoff)] <- 0
  sensitivity.call[which(cs >= data.cutoff)] <- 1
  
  return(list('sensitivity.call'=sensitivity.call, 'cutoff'=data.cutoff))
}

### Example
# library(PharmacoGx)
# data("GDSCsmall")
# #load("GatherData/PSets/GDSC.RData")
# 
# auc <- summarizeSensitivityProfiles(GDSCsmall, sensitivity.measure="auc_recomputed")
# auc.bin <- callSensitivity(auc)
