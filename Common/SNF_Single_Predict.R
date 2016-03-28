SNF_Single_Predict <- function(
  feature.sets, 
  parameters, 
  data, 
  partition, 
  ground_truth, 
  run_ind = 0, 
  NFOLDS = 10, 
  verbose = FALSE, 
  type_measure = "auc") {  
  
  require("foreach")
  require("doParallel")
  require("SNFtool")
  require("ROCR")
  
  if (is.null(feature.sets)) {
    feature.sets = 1:dim(data)[2]
  }
  
  ground_truth[which(ground_truth) == FALSE] = 0
  
  # Error handling
  stopifnot(typeof(feature.sets) == "integer")
  stopifnot(typeof(parameters) == "list")
  stopifnot(typeof(partition[[run_ind]]) == "list")  
  stopifnot(typeof(parameters) == "list")    
  stopifnot(range(ground_truth)[1] == 0)
  stopifnot(range(ground_truth)[2] == 1)
  stopifnot(length(parameters$K) > 1)
  stopifnot(run_ind > 0)
  
  f.ind <- feature.sets
  data.use = data[, f.ind]
  
  temp.repeat_fold <- TRUE
  temp.loop_count = 0
  while (temp.repeat_fold) {
    temp.folds = createFolds(partition[[run_ind]]$training_index, k = NFOLDS)
    temp.repeat_fold = FALSE
    for (i in 1:length(temp.folds)) {
      if (min(table(ground_truth[temp.folds[[i]]])) < 2) {
        temp.repeat_fold = TRUE
        break
      }
    }
    temp.loop_count = temp.loop_count + 1
    if (temp.loop_count > 100) {
      stop("create fold error in SNF predict")
    }
  }
  
  auc.training <- foreach (i = 1:length(parameters$K), .combine=rbind, .errorhandling="stop") %dopar% {    
    auc.training.one_split  <- foreach (cv_ind = 1:NFOLDS, .combine=rbind, .errorhandling="stop") %do% {
      
      temp.cv.index <- partition[[run_ind]]$training_index[temp.folds[[cv_ind]]]
      temp.training.index <- setdiff(partition[[run_ind]]$training_index, temp.cv.index)        
      
      temp.test_data = data.use[c(temp.cv.index, partition[[run_ind]]$test_index), ]
      temp.training_data = data.use[temp.training.index, ]
      #stopifnot(dim(temp.test_data)[1] + dim(temp.training_data)[1] == dim(data.use)[1])
      
      temp.ground_truth = ground_truth[temp.training.index]
      stopifnot(rownames(temp.training_data) == names(temp.ground_truth))      
      stopifnot(length(table(temp.ground_truth)) == 2)
      temp.predicted_labels <- SNF_LP_Cheng(train = temp.training_data, 
                                            test = temp.test_data, 
                                            groups = as.numeric(temp.ground_truth), 
                                            K = parameters$K[i], 
                                            alpha = 0.5, 
                                            t = 20, 
                                            method = 1)           
      
      temp.predicted_labels <- temp.predicted_labels[[2]]      
      stopifnot(length(temp.predicted_labels) == dim(temp.training_data)[1] + dim(temp.test_data)[1])      
      temp.predicted_labels <- temp.predicted_labels[-(1 : length(temp.training.index))]            
      
      temp.acc <- mean(round(temp.predicted_labels[(1 : length(temp.cv.index))]) == ground_truth[temp.cv.index])
      if (length(table(ground_truth[temp.cv.index])) == 2) {
        temp.prediction_results <- prediction(temp.predicted_labels[(1 : length(temp.cv.index))], ground_truth[temp.cv.index])
        temp.auc <- unlist(slot(performance(temp.prediction_results, "auc"), "y.values"))  
      } else {
        temp.auc <- temp.acc
      }
      c(temp.auc, temp.acc)
    }
    c(mean(auc.training.one_split[, 1]), mean(auc.training.one_split[, 2]))
  }
  
  print(type_measure)
  if (type_measure == "auc") {
    fs.training_auc <- max(auc.training[, 1])
    best_parameters.ind <- which(auc.training[, 1] == fs.training_auc, arr.ind=TRUE)
    best_parameters.ind <- best_parameters.ind[ceiling(length(best_parameters.ind)/2)]
    fs.training_acc <- auc.training[best_parameters.ind, 2]    
  } else {
    fs.training_acc <- max(auc.training[, 2])    
    best_parameters.ind <- which(auc.training[, 2] == fs.training_acc, arr.ind=TRUE)
    best_parameters.ind <- best_parameters.ind[ceiling(length(best_parameters.ind)/2)]
    fs.training_auc <- auc.training[best_parameters.ind, 1]
  }
  
  print(paste("Best K:", parameters$K[best_parameters.ind]))        
  if (verbose) {      
    print(paste("Training AUC: ", fs.training_auc))    
    print(paste("Training accuracy: ", fs.training_acc))        
  }
  
  temp.training_data = data.use[partition[[run_ind]]$training_index, ]
  temp.test_data = data.use[c(partition[[run_ind]]$ES_index, partition[[run_ind]]$test_index), ]
  
  #stopifnot(length(unique(c(rownames(temp.training_data), rownames(temp.test_data)))) == dim(data.use)[1])  
  stopifnot(rownames(temp.test_data) == names(c(partition[[run_ind]]$ES_index, partition[[run_ind]]$test_index)))
  
  temp.predicted_labels <- SNF_LP_Cheng(temp.training_data, temp.test_data, ground_truth[partition[[run_ind]]$training_index], 
                                        K = parameters$K[best_parameters.ind], alpha = 0.5, t = 20, method = 1)             
  temp.all_predict <- temp.predicted_labels[[2]][-(1:(length(partition[[run_ind]]$training_index)))] 
  fs.test_prediction <- tail(temp.all_predict, length(partition[[run_ind]]$test_index))
  fs.ES_prediction <- head(temp.all_predict, length(partition[[run_ind]]$ES_index))
  
  stopifnot(range(fs.test_prediction)[1] >= 0)
  stopifnot(range(fs.test_prediction)[2] <= 1)
  stopifnot(length(fs.test_prediction) == length(partition[[run_ind]]$test_index))  
  
  if (length(fs.test_prediction) > 1) {
    if (length(table(ground_truth[partition[[run_ind]]$test_index])) == 1) {
      fs.test_auc <- 0
      print(paste("Test AUC: ", fs.test_auc))    
    } else {    
    temp.test_predict <- prediction(fs.test_prediction, ground_truth[partition[[run_ind]]$test_index])
    fs.test_auc <- unlist(slot(performance(temp.test_predict, "auc"), "y.values"))
    print(paste("Test AUC: ", fs.test_auc))    
    }
  }
  fs.test_accuracy <- mean(round(fs.test_prediction) == ground_truth[partition[[run_ind]]$test_index])        
  
  if (verbose) {
    print(paste("Test accuracy: ", fs.test_accuracy))    
  }   
  
  if (exists("fs.test_auc")) {    
    return (
      list(
        training.auc = fs.training_auc
        , training.accuracy = fs.training_acc
        , ES.prediction = fs.ES_prediction
        , test.prediction = fs.test_prediction
        , test.auc = fs.test_auc
        , test.accuracy = fs.test_accuracy
        , K = parameters$K[best_parameters.ind]
      )
    )
  } else {
    return (
      list(
        training.auc = fs.training_auc
        , training.accuracy = fs.training_acc
        , ES.prediction = fs.ES_prediction
        , test.prediction = fs.test_prediction        
        , test.accuracy = fs.test_accuracy
        , K = parameters$K[best_parameters.ind]
      )
    )    
  }
}
