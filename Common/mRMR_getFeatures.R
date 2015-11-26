mRMR_getFeatures <- function(data, labels, feature_count, solution_count = 1) {
  
  # error checking
  stopifnot(dim(data)[1] >= 10)
  stopifnot(dim(data)[2] > 10)
  stopifnot(dim(data)[1] == length(labels))
  stopifnot(feature_count > 10)
  stopifnot(solution_count >= 1)
  
  require(mRMRe)
  set.thread.count(8)
  
  temp.mRMR.data <- mRMR.data(data = data.frame(labels, data))  
  gc()
  temp.mRMR.result <- mRMR.ensemble(data = temp.mRMR.data, target_indices = c(1), solution_count, feature_count)  
  temp.mRMR.features = unique(as.vector(temp.mRMR.result@filters[[1]])) - 1
  
  return(temp.mRMR.features)
}
