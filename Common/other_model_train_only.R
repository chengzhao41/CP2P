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
  
  # Setup ---------------------------------    
  # 1) Elastic Net Logistic Regression
  elastic_net.cv_error = vector()
  elastic_net.cv_model = list()
  elastic_net.ALPHA <- c(1:9) / 10
  
  temp.cv_error_matrix <- foreach (temp = 1:N_CV_REPEATS, .combine=rbind, .errorhandling="stop") %do% {      
    for (alpha_index in 1:length(elastic_net.ALPHA))
    {      
      elastic_net.cv_model[[alpha_index]] = cv.glmnet(data[partition[[run_ind]]$training_index, selected_features],
                                                      ground_truth[partition[[run_ind]]$training_index],
                                                      alpha = elastic_net.ALPHA[alpha_index]
                                                      , type.measure = type_measure.glmnet
                                                      , family = "binomial"
                                                      , standardize = FALSE 
                                                      , nfolds = NFOLDS
                                                      , nlambda = 100
                                                      , parallel = TRUE
      )
      elastic_net.cv_error[alpha_index] = min(elastic_net.cv_model[[alpha_index]]$cvm)
    }
    elastic_net.cv_error    
  }
  
  if (N_CV_REPEATS == 1) {
    temp.cv_error_mean = temp.cv_error_matrix
  } else {
    temp.cv_error_mean = apply(temp.cv_error_matrix, 2, mean)  
  }
  
  stopifnot(length(temp.cv_error_mean) == length(elastic_net.ALPHA))
  temp.best_alpha_index = which(min(temp.cv_error_mean) == temp.cv_error_mean)[length(which(min(temp.cv_error_mean) == temp.cv_error_mean))] 
  print(paste("Best ALPHA:", elastic_net.ALPHA[temp.best_alpha_index]))
  
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  while (temp.non_zero_coeff < 3) {
    elastic_net.cv_model = cv.glmnet(
      data[partition[[run_ind]]$training_index, selected_features]
      , ground_truth[partition[[run_ind]]$training_index]
      , alpha = elastic_net.ALPHA[temp.best_alpha_index]
      , type.measure = type_measure.glmnet
      , family = "binomial"
      , standardize=FALSE, 
      , nlambda = 100
      , nfolds = NFOLDS
      , parallel = TRUE
    )
    temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min)
    temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index]    
    temp.loop_count = temp.loop_count + 1
    as.numeric(Sys.time())-> t 
    set.seed((t - floor(t)) * 1e8 -> seed) 
    #print(paste0("seed: ", seed))
    if (temp.loop_count > 5) {
      print("diverged")
      temp.min_lambda_index = 50
      break
    }
  } 
  print(temp.non_zero_coeff)  
  
  elasticNet.model = glmnet(data[partition[[run_ind]]$training_index, selected_features], 
                            ground_truth[partition[[run_ind]]$training_index], 
                            alpha = elastic_net.ALPHA[temp.best_alpha_index],
                            standardize=FALSE,
                            nlambda = 100,
                            family = "binomial")

  # 2) Lasso Logistic Regression  
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  while (temp.non_zero_coeff < 3) {      
    temp.cv_model = cv.glmnet(data[partition[[run_ind]]$training_index, selected_features]
                              , ground_truth[partition[[run_ind]]$training_index]
                              , alpha = 1
                              , type.measure = type_measure.glmnet
                              , family = "binomial"
                              , standardize = FALSE
                              , nlambda = 100
                              , nfolds = NFOLDS
                              , parallel = TRUE
    )
    temp.min_lambda_index = which(temp.cv_model$lambda == temp.cv_model$lambda.min)
    temp.non_zero_coeff = temp.cv_model$nzero[temp.min_lambda_index]
    temp.loop_count = temp.loop_count + 1
    as.numeric(Sys.time())-> t 
    set.seed((t - floor(t)) * 1e8 -> seed) 
    #print(paste0("seed: ", seed))
    if (temp.loop_count > 5) {
      print("diverged")
      temp.min_lambda_index = 50
      break
    }    
  }  
  print(temp.non_zero_coeff)  
  
  lasso.model = glmnet(data[partition[[run_ind]]$training_index, selected_features], 
                       ground_truth[partition[[run_ind]]$training_index], 
                       alpha = 1,
                       standardize=FALSE,
                       nlambda = 100,
                       family = "binomial")
  
  
  # 3) Ridge Regularized Logistic Regression  
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  while (temp.non_zero_coeff < 3) {      
    temp.cv_model = cv.glmnet(data[partition[[run_ind]]$training_index, selected_features]
                              , ground_truth[partition[[run_ind]]$training_index]
                              , alpha = 0
                              , type.measure = type_measure.glmnet
                              , family = "binomial"
                              , standardize = FALSE
                              , nlambda = 100
                              , nfolds = NFOLDS
                              , parallel = FALSE
    )
    temp.min_lambda_index = which(temp.cv_model$lambda == temp.cv_model$lambda.min)
    temp.non_zero_coeff = temp.cv_model$nzero[temp.min_lambda_index]
    temp.loop_count = temp.loop_count + 1
    as.numeric(Sys.time())-> t 
    set.seed((t - floor(t)) * 1e8 -> seed) 
    print(paste0("seed: ", seed))
    print(temp.loop_count)
    if (temp.loop_count > 5) {
      print("diverged")
      temp.min_lambda_index = 50
      break
    }        
  } 
  print(temp.non_zero_coeff)  
  
  ridge.model = glmnet(data[partition[[run_ind]]$training_index, selected_features], 
                       ground_truth[partition[[run_ind]]$training_index], 
                       alpha = 0,
                       standardize=FALSE,
                       nlambda = 100,
                       family = "binomial")
  
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
