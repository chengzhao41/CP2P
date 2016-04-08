Other_Model_Predict_Only <- function(data, ground_truth, partition, selected_features = NULL, NFOLDS = 5, N_CV_REPEATS = 2, run_ind = 0, type_measure = "auc", input.models) {
  
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
  
  # Setup ---------------------------------    
  # 1) Elastic Net Logistic Regression

  elasticNet.model = input.models$elasticNet.model
  
  if (length(partition[[run_ind]]$test_index) > 1) {
    temp.predictions = predict(elasticNet.model, data[partition[[run_ind]]$test_index, selected_features], type = "response")  
    elasticNet.predictions <- temp.predictions[, temp.min_lambda_index]  
  } else {
    temp.predictions = predict(elasticNet.model, data[, selected_features], type = "response")  
    elasticNet.predictions <- temp.predictions[partition[[run_ind]]$test_index, temp.min_lambda_index]      
  }
  temp.l <- min(length(elastic_net.cv_model$lambda), length(elasticNet.model$lambda)) 
  stopifnot(elastic_net.cv_model$lambda[1:temp.l] == elasticNet.model$lambda[1:temp.l])  
  
  if (type_measure == "auc" && length(table(ground_truth[partition[[run_ind]]$test_index])) == 2) {
    temp.predictions <- prediction(elasticNet.predictions, ground_truth[partition[[run_ind]]$test_index])
    elasticNet.test_auc <- unlist(slot(performance(temp.predictions, type_measure.performance), "y.values"))          
    print(paste("Elastic Net test AUC:", elasticNet.test_auc))  
  } else {
    elasticNet.test_auc = sum(round(elasticNet.predictions) == ground_truth[partition[[run_ind]]$test_index]) / length(elasticNet.predictions)
    print(paste("Elastic Net test acc:", elasticNet.test_auc))  
  }
  
  elasticNet = list(    
    predict = elasticNet.predictions  
    , test_auc = elasticNet.test_auc    
    , alpha = elastic_net.ALPHA[temp.best_alpha_index]        
  )
  
  rm(list = ls(pattern="elasticNet."))
  rm(list = ls(pattern="elastic_net."))
  rm(list = ls(pattern="temp"))
  rm(list = ls(pattern="alpha."))  
  
  # 2) Lasso Logistic Regression  
  lasso.model = input.models$lasso.model
  
  if (length(partition[[run_ind]]$test_index) > 1) {
    temp.predictions = predict(lasso.model, data[partition[[run_ind]]$test_index, selected_features], type = "response")  
    lasso.predictions <- temp.predictions[, temp.min_lambda_index]  
  } else {
    temp.predictions = predict(lasso.model, data[, selected_features], type = "response")  
    lasso.predictions <- temp.predictions[partition[[run_ind]]$test_index, temp.min_lambda_index]      
  }  
  
  if (type_measure == "auc" && length(table(ground_truth[partition[[run_ind]]$test_index])) == 2) {
    temp.predictions <- prediction(lasso.predictions, ground_truth[partition[[run_ind]]$test_index])
    lasso.test_auc <- unlist(slot(performance(temp.predictions, type_measure.performance), "y.values"))          
    print(paste("lasso test AUC:", lasso.test_auc))  
  } else {
    lasso.test_auc = sum(round(lasso.predictions) == ground_truth[partition[[run_ind]]$test_index]) / length(lasso.predictions)
    print(paste("lasso test acc:", lasso.test_auc))  
  }
  
  lasso = list(    
    predict = lasso.predictions
    , test_auc = lasso.test_auc    
  )
  
  rm(list = ls(pattern="lasso."))  
  rm(list = ls(pattern="temp"))
  
  # 3) Ridge Regularized Logistic Regression  
  ridge.model = input.models$ridge.model
  
  if (length(partition[[run_ind]]$test_index) > 1) {
    temp.predictions = predict(ridge.model, data[partition[[run_ind]]$test_index, selected_features], type = "response")  
    ridge.predictions <- temp.predictions[, temp.min_lambda_index]  
  } else {
    temp.predictions = predict(ridge.model, data[, selected_features], type = "response")  
    ridge.predictions <- temp.predictions[partition[[run_ind]]$test_index, temp.min_lambda_index]      
  }  
  
  if (type_measure == "auc" && length(table(ground_truth[partition[[run_ind]]$test_index])) == 2) {
    temp.predictions <- prediction(ridge.predictions, ground_truth[partition[[run_ind]]$test_index])
    ridge.test_auc <- unlist(slot(performance(temp.predictions, type_measure.performance), "y.values"))          
    print(paste("Ridge test AUC:", ridge.test_auc))  
  } else {
    ridge.test_auc = sum(round(ridge.predictions) == ground_truth[partition[[run_ind]]$test_index]) / length(ridge.predictions)
    print(paste("Ridge test acc:", ridge.test_auc))  
  }
  
  ridge = list(    
    predict = ridge.predictions   
    , test_auc = ridge.test_auc
  )
  
  rm(list = ls(pattern="ridge."))  
  rm(list = ls(pattern="temp"))  
  
  # 4) Random Forest  
  rf.model <- input.models$rf.model
  
  if (length(partition[[run_ind]]$test_index) > 1) {
    rf.predict <- predict(rf.model 
                          , newdata = data[partition[[run_ind]]$test_index, selected_features]
                          , type = "prob")
  } else {
    rf.predict <- predict(rf.model 
                          , newdata = data[, selected_features]
                          , type = "prob")
    rf.predict <- rf.predict[partition[[run_ind]]$test_index, ]
  }
  
  if (type_measure == "auc" && length(table(ground_truth[partition[[run_ind]]$test_index])) == 2) { 
    temp.predict <- prediction(rf.predict$b, ground_truth[partition[[run_ind]]$test_index])  
    rf.test_auc <- unlist(slot(performance(temp.predict, "auc"), "y.values"))  
    print(paste("RF test AUC:", rf.test_auc))  
  } else {
    rf.test_auc <- sum(round(rf.predict$b) == ground_truth[partition[[run_ind]]$test_index]) / length(rf.predict$b)
    print(paste("RF test acc:", rf.test_auc))  
  }
  
  rf = list(
    predict = rf.predict$b
    , test_auc = rf.test_auc
  )
  
  rm(list = ls(pattern="rf."))
  rm(list = ls(pattern="temp."))
  
  # 5) SVM with linear kernel  
  svmLinear.model <- input.models$svmLinear.model
    
  if (length(partition[[run_ind]]$test_index) > 1) {
    svmLinear.predict <- predict(svmLinear.model 
                                 , newdata = data[partition[[run_ind]]$test_index, selected_features]
                                 , type = "prob")
  } else {
    svmLinear.predict <- predict(svmLinear.model 
                                 , newdata = data[, selected_features]
                                 , type = "prob")
    svmLinear.predict <- svmLinear.predict[partition[[run_ind]]$test_index, ]
  }
  
  if (type_measure == "auc" && length(table(ground_truth[partition[[run_ind]]$test_index])) == 2) {
    temp.predict <- prediction(svmLinear.predict$b, ground_truth[partition[[run_ind]]$test_index])  
    svmLinear.test_auc <- unlist(slot(performance(temp.predict, "auc"), "y.values"))  
    print(paste("SVM linear test AUC:", svmLinear.test_auc))     
  } else {
    svmLinear.test_auc <- sum(round(svmLinear.predict$b) == ground_truth[partition[[run_ind]]$test_index]) / length(svmLinear.predict$b)
    print(paste("SVM linear test acc:", svmLinear.test_auc))    
  }  
  
  svmLinear = list(    
    predict = svmLinear.predict$b
    , test_auc = svmLinear.test_auc
  )
  
  rm(list = ls(pattern="svmLinear."))
  rm(list = ls(pattern="temp."))
  
  # 6) SVM with radial kernel
  svmRadial.model <- input.models$svmRadial.model
  
  if (length(partition[[run_ind]]$test_index) > 1) {
    svmRadial.predict <- predict(svmRadial.model 
                                 , newdata = data[partition[[run_ind]]$test_index, selected_features]
                                 , type = "prob")
  } else {
    svmRadial.predict <- predict(svmRadial.model 
                                 , newdata = data[, selected_features]
                                 , type = "prob")
    svmRadial.predict <- svmRadial.predict[partition[[run_ind]]$test_index, ]
  }
  
  if (type_measure == "auc" && length(table(ground_truth[partition[[run_ind]]$test_index])) == 2) { 
    temp.predict <- prediction(svmRadial.predict$b, ground_truth[partition[[run_ind]]$test_index])  
    svmRadial.test_auc <- unlist(slot(performance(temp.predict, "auc"), "y.values"))  
    print(paste("SVM radial test AUC:", svmRadial.test_auc))     
  } else {
    svmRadial.test_auc <- sum(round(svmRadial.predict$b) == ground_truth[partition[[run_ind]]$test_index]) / length(svmRadial.predict$b)
    print(paste("SVM radial test acc:", svmRadial.test_auc))    
  }  
  
  svmRadial = list(    
    predict = svmRadial.predict$b
    , test_auc = svmRadial.test_auc
  )
  
  rm(list = ls(pattern="svmRadial."))
  rm(list = ls(pattern="temp."))
  
  return (list(
    elasticNet = elasticNet
    , lasso = lasso
    , ridge = ridge
    , rf = rf
    , svmLinear = svmLinear
    , svmRadial = svmRadial
  ))
}
