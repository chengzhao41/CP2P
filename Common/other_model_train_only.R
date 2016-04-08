Other_Model_Train_Only <- function(data, ground_truth, partition, selected_features = NULL, NFOLDS = 5, N_CV_REPEATS = 2, run_ind = 0, type_measure = "auc") {
  
  # uses all samples in CV, as opposed to Other_Model_Predict2 which uses only patients
  
  # Checking for input parameters ---------------------------------
  require("glmnet")
  require("foreach")
  require("doParallel")
  require("ROCR")
  require("caret")
  
  stopifnot(abs(mean(data)) < 0.01)
  stopifnot(dim(data)[1] > 10)
  stopifnot(dim(data)[2] > 10)
  #stopifnot(!is.null(partition$patients))
  #stopifnot(!is.null(partition$cell_lines))
  stopifnot(typeof(partition) == "list")  
  #stopifnot(length(ground_truth) == dim(data)[1])  
  
  if (is.null(selected_features)) {
    selected_features = 1:dim(data)[2]
  }
  
  type_measure.glmnet = "auc"
  type_measure.other_models = "ROC"
  type_measure.performance = "auc"   
  
  if (length(partition[[run_ind]]$training_index) < 30) {
    type_measure.glmnet = "class"
  }
  
  if (type_measure == "auc") {
    ground_truth <- as.factor(ground_truth)  
    stopifnot(length(levels(ground_truth)) == 2)
    levels(ground_truth) = c("a", "b")
  } else if (type_measure == "acc") {
    ground_truth <- as.numeric(ground_truth)  
    type_measure.glmnet = "class"
    type_measure.other_models = "ROC"
    type_measure.performance = "acc"    
  }
  
  # 4) Random Forest  
  other_model.y = as.factor(ground_truth[partition[[run_ind]]$training_index])
  levels(other_model.y) <- c("a", "b")
  
  fitControl <- trainControl(
    method = "repeatedcv",
    number = min(10, NFOLDS),
    classProbs = TRUE,     
    repeats = 1,
    allowParallel = TRUE,
    summaryFunction = twoClassSummary)   
  
  rf.model <- train(x = data[partition[[run_ind]]$training_index, selected_features]
                    , y = other_model.y
                    , method = "rf"
                    , trControl = fitControl                   
                    , verbose = FALSE
                    , metric = "ROC")
  
  # 5) SVM with linear kernel  
  svmLinear.model <- train(x = data[partition[[run_ind]]$training_index, selected_features]
                           , y = other_model.y
                           , method = "svmLinear"
                           , trControl = fitControl                   
                           , verbose = FALSE
                           , metric = type_measure.other_models                 
  )
  
  # 6) SVM with radial kernel
  svmRadial.model <- train(x = data[partition[[run_ind]]$training_index, selected_features]
                           , y = other_model.y
                           , method = "svmRadial"
                           , trControl = fitControl
                           , verbose = FALSE
                           , metric = type_measure.other_models
  )
  
  return (list(
    elasticNet.model = elasticNet.model
    , lasso.model = lasso.model
    , ridge.model = ridge.model
    , rf.model = rf.model
    , svmLinear.model = svmLinear.model
    , svmRadial.model = svmRadial.model
  ))
}
