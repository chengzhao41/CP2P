# Load Libraries and Dependencies ---------------------------------
library("SNFtool")
library("igraph")
library("foreach")
library("sva")
library("ROCR")
library('doParallel')
library("caret")
library("glmnet")
library("randomForest")
library("kernlab")
library("pROC")
library(preprocessCore)
rm(list = ls(pattern="temp*"))

if (Sys.info()["nodename"] == "cheng-OptiPlex-990") {
  setwd("/home/cheng/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Script")  
  input_partition_BEGIN = 1
  input_partition_END = 1  
  registerDoParallel(8)
  load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Epirubicin/WS/epirubicin_data.RData")
} else {
  setwd("/home/zhaoche7/drp2/Epirubicin/WS")  
  load("epirubicin_data.RData")
  registerDoParallel(8)
  args <- commandArgs(trailingOnly = TRUE)
  input_partition_BEGIN = as.integer(args[1])
  input_partition_END = as.integer(args[2])
  print(paste("begin", input_partition_BEGIN, "end", input_partition_END))
  setwd("/home/zhaoche7/drp/Script")
}

if (args[3] == "slope") {
  input_data <- epirubicin$combined_slope.sva
  input_label <- epirubicin.labels$slope_combined
  input_partition <- partition$slope
  feature.l1000 <- feature.l1000$cp
} else if (args[3] == "auc") {
  input_data <- epirubicin$combined_AUC.sva
  input_label <- epirubicin.labels$AUC_combined
  input_partition <- partition$AUC
  feature.l1000 <- feature.l1000$cp
} else {
  stopifnot(FALSE)
}

if (length(args) == 4) {
  training_var_amount <- as.integer(args[4])
  print(paste0("epirubicin_cpp_", input_partition_BEGIN, "to", input_partition_END, "_", args[3], "_var", training_var_amount, ".RData"))  
  if (args[3] == "slope") {    
    input_partition$cpp <- partition_var$cpp[[training_var_amount]]$slope
  } else if (args[3] == "auc") {    
    input_partition$cpp <- partition_var$cpp[[training_var_amount]]$AUC
  } else {
    stop("arg 3 is wrong")
  }  
} else {
  stop("args != 4")
}

rm(epirubicin)

# Other models make predictions ---------------------------------
source('Other_Model_Predict.R')
source('SNF_Single_Predict.R')
source('SNF_LP_Cheng.R')
source('mRMR_getFeatures.R')

snf.parameter <- seq(from = 5, to = 30, by = 5)

# mRMR + SNF ---------------------------------

cpp.snf.single.mRMR1000 = list()
cpp.other_model.mRMR1000 = list()

for (temp.run_ind in input_partition_BEGIN:input_partition_END) {  
  gc()
  print(Sys.time())
  print(temp.run_ind)
  
  temp.mRMR_features <- mRMR_getFeatures(
    input_data[input_partition$cpp[[temp.run_ind]]$training_index.single, ], 
    as.ordered(input_label[input_partition$cpp[[temp.run_ind]]$training_index.single]), 
    feature_count = 1000, 
    solution_count = 1)
  print(paste("Number of mRMR features:", length(temp.mRMR_features)))
  temp.mRMR_features <- as.integer(temp.mRMR_features)
  
  res <- try(  
    cpp.snf.single.mRMR1000[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = temp.mRMR_features, 
                                                                  parameters = list(K = snf.parameter), 
                                                                  data = input_data, 
                                                                  partition = input_partition$cpp,
                                                                  ground_truth = input_label,
                                                                  run_ind = temp.run_ind,
                                                                  NFOLDS = 5,
                                                                  type_measure = "auc")
  )
  temp.error_count = 0  
  while (inherits(res, "try-error") && temp.error_count < 5) {
    temp.error_count = temp.error_count + 1    
    res <- try(  
      cpp.snf.single.mRMR1000[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = temp.mRMR_features, 
                                                                    parameters = list(K = snf.parameter), 
                                                                    data = input_data, 
                                                                    partition = input_partition$cpp,
                                                                    ground_truth = input_label,
                                                                    run_ind = temp.run_ind,
                                                                    NFOLDS = 5,
                                                                    type_measure = "auc")
    )
  }
  cpp.snf.single.mRMR1000[[temp.run_ind]]$feature.sets = temp.mRMR_features
  
  res <- try(    
    cpp.other_model.mRMR1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                                    ground_truth = input_label, 
                                                                    partition = input_partition$cpp,
                                                                    selected_features = temp.mRMR_features, 
                                                                    NFOLDS = 5, 
                                                                    N_CV_REPEATS = 1, 
                                                                    run_ind = temp.run_ind, 
                                                                    type_measure = "auc")
  )
  temp.error_count = 0  
  while (inherits(res, "try-error") && temp.error_count < 5) {
    res <- try(    
      cpp.other_model.mRMR1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                                      ground_truth = input_label, 
                                                                      partition = input_partition$cpp,
                                                                      selected_features = temp.mRMR_features, 
                                                                      NFOLDS = 5, 
                                                                      N_CV_REPEATS = 1, 
                                                                      run_ind = temp.run_ind, 
                                                                      type_measure = "auc")
    )   
  }
}

cpp.other_model.all = list()
for (temp.run_ind in input_partition_BEGIN:input_partition_END) {
  gc()
  print(Sys.time())
  print(temp.run_ind)
  
  res <- try( 
    cpp.other_model.all[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                               ground_truth = input_label, 
                                                               partition = input_partition$cpp,
                                                               selected_features = NULL, 
                                                               NFOLDS = 5, 
                                                               N_CV_REPEATS = 1, 
                                                               run_ind = temp.run_ind,
                                                               type_measure = "auc")
  )
  temp.error_count = 0
  while (inherits(res, "try-error") && temp.error_count < 5) {
    temp.error_count = temp.error_count + 1
    res <- try( 
      cpp.other_model.all[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                                 ground_truth = input_label, 
                                                                 partition = input_partition$cpp,
                                                                 selected_features = NULL, 
                                                                 NFOLDS = 5, 
                                                                 N_CV_REPEATS = 1, 
                                                                 run_ind = temp.run_ind,
                                                                 type_measure = "auc")
    )
  }
}

cpp.other_model.l1000 = list()
for (temp.run_ind in input_partition_BEGIN:input_partition_END) {
  gc()
  print(Sys.time())
  print(temp.run_ind)
  res <- try( 
    cpp.other_model.l1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                                 ground_truth = input_label, 
                                                                 partition = input_partition$cpp,
                                                                 selected_features = feature.l1000, 
                                                                 NFOLDS = 5, 
                                                                 N_CV_REPEATS = 1, 
                                                                 run_ind = temp.run_ind, 
                                                                 type_measure = "auc")
  )
  temp.error_count = 0
  while (inherits(res, "try-error") && temp.error_count < 5) {
    temp.error_count = temp.error_count + 1
    res <- try( 
      cpp.other_model.l1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                                   ground_truth = input_label, 
                                                                   partition = input_partition$cpp,
                                                                   selected_features = feature.l1000, 
                                                                   NFOLDS = 5, 
                                                                   N_CV_REPEATS = 1, 
                                                                   run_ind = temp.run_ind, 
                                                                   type_measure = "auc")
    )
  }
}

cpp.snf.single.all = list()

for (temp.run_ind in input_partition_BEGIN:input_partition_END) {  
  gc()
  print(Sys.time())
  print(temp.run_ind)
  res <- try( 
    cpp.snf.single.all[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = NULL, 
                                                             parameters = list(K = snf.parameter), 
                                                             data = input_data, 
                                                             partition = input_partition$cpp,
                                                             ground_truth = input_label,
                                                             run_ind = temp.run_ind,
                                                             NFOLDS = 5, 
                                                             type_measure = "auc")
  )
  temp.error_count = 0  
  while (inherits(res, "try-error") && temp.error_count < 5) {
    temp.error_count = temp.error_count + 1
    res <- try( 
      cpp.snf.single.all[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = NULL, 
                                                               parameters = list(K = snf.parameter), 
                                                               data = input_data, 
                                                               partition = input_partition$cpp,
                                                               ground_truth = input_label,
                                                               run_ind = temp.run_ind,
                                                               NFOLDS = 5, 
                                                               type_measure = "auc")      
    )
  }
}

cpp.snf.single.l1000 = list()
for (temp.run_ind in input_partition_BEGIN:input_partition_END) {  
  gc()
  print(Sys.time())
  print(temp.run_ind)
  
  res <- try(   
    cpp.snf.single.l1000[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = feature.l1000, 
                                                               parameters = list(K = snf.parameter), 
                                                               data = input_data, 
                                                               partition = input_partition$cpp,
                                                               ground_truth = input_label,
                                                               run_ind = temp.run_ind,
                                                               NFOLDS = 5,
                                                               type_measure = "auc")
  )
  temp.error_count = 0  
  while (inherits(res, "try-error") && temp.error_count < 5) {
    temp.error_count = temp.error_count + 1
    res <- try(   
      cpp.snf.single.l1000[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = feature.l1000, 
                                                                 parameters = list(K = snf.parameter), 
                                                                 data = input_data, 
                                                                 partition = input_partition$cpp,
                                                                 ground_truth = input_label,
                                                                 run_ind = temp.run_ind,
                                                                 NFOLDS = 5,
                                                                 type_measure = "auc")
    )
  }
}

print("completed!")
print(paste("BEGIN and END:", input_partition_BEGIN, input_partition_END))
setwd("/home/zhaoche7/")
rm(input_data)
if (exists("training_var_amount")) {
  save.image(paste0("epirubicin_cpp_", input_partition_BEGIN, "to", input_partition_END, "_", args[3], "_var", training_var_amount, ".RData"))
} else {  
  save.image(paste0("epirubicin_cpp_", input_partition_BEGIN, "to", input_partition_END, "_", args[3], ".RData"))
}

