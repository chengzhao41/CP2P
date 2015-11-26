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
  PARTITION_BEGIN = 1
  PARTITION_END = 1  
  registerDoParallel(7)
  load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Epirubicin/WS/epirubicin_data.RData")
} else {
  setwd("/home/zhaoche7/drp2/Epirubicin/WS")  
  load("epirubicin_data.RData")
  registerDoParallel(8)
  args <- commandArgs(trailingOnly = TRUE)
  PARTITION_BEGIN = as.integer(args[1])
  PARTITION_END = as.integer(args[2])
  print(paste("begin", PARTITION_BEGIN, "end", PARTITION_END))
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
} else if (args[3] == "ic50") {
  input_data <- epirubicin$combined_IC50.sva
  input_label <- epirubicin.labels$IC50_combined
  input_partition <- partition$IC50
  feature.l1000 <- feature.l1000$cp
} else {
  stop("invalid args[3]")
  input_data <- epirubicin$combined_AUC.sva
  input_label <- epirubicin.labels$AUC_combined
  input_partition <- partition$AUC
  feature.l1000 <- feature.l1000$cp
  input_partition <- list()
  input_partition$cp <- list(list(test_index = which(epirubicin.labels$slope_combined.source == "patient"), 
                                  training_index.single = which(epirubicin.labels$slope_combined.source != "patient")))
}

if (length(args) == 4) {
  training_var_amount <- as.integer(args[[4]])
  if (args[3] == "slope") {    
    input_partition$cp <- partition_var$cp[[training_var_amount]]$slope
  } else if (args[3] == "auc") {    
    input_partition$cp <- partition_var$cp[[training_var_amount]]$AUC
  } else if (args[3] == "ic50") {    
    input_partition$cp <- partition_var$cp[[training_var_amount]]$IC50
  } else {
    stop("arg 3 is wrong")
  }
  print(paste0("epirubicin_cp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))
} else {
  stop("args not 4")
}

rm(epirubicin)

# Other models make predictions ---------------------------------
source('Other_Model_Predict.R')
source('SNF_Single_Predict.R')
source('SNF_LP_Cheng.R')
source('mRMR_getFeatures.R')

cp.snf.single.mRMR1000 = list()
cp.other_model.mRMR1000 = list()
snf.parameter <- seq(from = 5, to = 30, by = 5)

# mRMR + SNF ---------------------------------

for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {  
  print(Sys.time())
  print(temp.run_ind)
  gc()
  temp.mRMR_features <- mRMR_getFeatures(
    input_data[input_partition$cp[[temp.run_ind]]$training_index.single, ], 
    as.ordered(input_label[input_partition$cp[[temp.run_ind]]$training_index.single]), 
    feature_count = 1000, 
    solution_count = 1)
  print(paste("Number of mRMR features:", length(temp.mRMR_features)))
  temp.mRMR_features <- as.integer(temp.mRMR_features)
  gc()
  res <- try(  
    cp.snf.single.mRMR1000[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = temp.mRMR_features, 
                                                                 parameters = list(K = snf.parameter), 
                                                                 data = input_data, 
                                                                 partition = input_partition$cp,
                                                                 ground_truth = input_label,
                                                                 run_ind = temp.run_ind,
                                                                 NFOLDS = 5,
                                                                 type_measure = "auc")
  )
  temp.error_count = 0  
  while (inherits(res, "try-error") && temp.error_count < 5) {
    temp.error_count = temp.error_count + 1    
    res <- try(  
      cp.snf.single.mRMR1000[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = temp.mRMR_features, 
                                                                   parameters = list(K = snf.parameter), 
                                                                   data = input_data, 
                                                                   partition = input_partition$cp,
                                                                   ground_truth = input_label,
                                                                   run_ind = temp.run_ind,
                                                                   NFOLDS = 5,
                                                                   type_measure = "auc")
    )
  }
  cp.snf.single.mRMR1000[[temp.run_ind]]$feature.sets = temp.mRMR_features
  
  res <- try(    
    cp.other_model.mRMR1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                                   ground_truth = input_label, 
                                                                   partition = input_partition$cp,
                                                                   selected_features = temp.mRMR_features, 
                                                                   NFOLDS = 5, 
                                                                   N_CV_REPEATS = 1, 
                                                                   run_ind = temp.run_ind, 
                                                                   type_measure = "auc")
  )
  temp.error_count = 0  
  while (inherits(res, "try-error") && temp.error_count < 5) {
    res <- try(    
      cp.other_model.mRMR1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                                     ground_truth = input_label, 
                                                                     partition = input_partition$cp,
                                                                     selected_features = temp.mRMR_features, 
                                                                     NFOLDS = 5, 
                                                                     N_CV_REPEATS = 1, 
                                                                     run_ind = temp.run_ind, 
                                                                     type_measure = "auc")
    )   
  }
}




cp.other_model.all = list()
for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {
  print(Sys.time())
  print(temp.run_ind)
  
  res <- try( 
  cp.other_model.all[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                            ground_truth = input_label, 
                                                            partition = input_partition$cp,
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
      cp.other_model.all[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                                ground_truth = input_label, 
                                                                partition = input_partition$cp,
                                                                selected_features = NULL, 
                                                                NFOLDS = 5, 
                                                                N_CV_REPEATS = 1, 
                                                                run_ind = temp.run_ind,
                                                                type_measure = "auc")
    )
  }
}

cp.other_model.l1000 = list()
for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {
  print(Sys.time())
  print(temp.run_ind)
  res <- try( 
  cp.other_model.l1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                              ground_truth = input_label, 
                                                              partition = input_partition$cp,
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
      cp.other_model.l1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                                  ground_truth = input_label, 
                                                                  partition = input_partition$cp,
                                                                  selected_features = feature.l1000, 
                                                                  NFOLDS = 5, 
                                                                  N_CV_REPEATS = 1, 
                                                                  run_ind = temp.run_ind, 
                                                                  type_measure = "auc")
    )
  }
}


cp.snf.single.all = list()

for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {  
  gc()
  print(Sys.time())
  print(temp.run_ind)
  res <- try( 
  cp.snf.single.all[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = NULL, 
                                                          parameters = list(K = snf.parameter), 
                                                          data = input_data, 
                                                          partition = input_partition$cp,
                                                          ground_truth = input_label,
                                                          run_ind = temp.run_ind,
                                                          NFOLDS = 5, 
                                                          type_measure = "auc")
  )
  temp.error_count = 0  
  while (inherits(res, "try-error") && temp.error_count < 5) {
    temp.error_count = temp.error_count + 1
    res <- try( 
      cp.snf.single.all[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = NULL, 
                                                              parameters = list(K = snf.parameter), 
                                                              data = input_data, 
                                                              partition = input_partition$cp,
                                                              ground_truth = input_label,
                                                              run_ind = temp.run_ind,
                                                              NFOLDS = 5, 
                                                              type_measure = "auc")      
    )
  }
}

cp.snf.single.l1000 = list()
for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {  
  gc()
  print(Sys.time())
  print(temp.run_ind)
  
  res <- try(   
  cp.snf.single.l1000[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = feature.l1000, 
                                                            parameters = list(K = snf.parameter), 
                                                            data = input_data, 
                                                            partition = input_partition$cp,
                                                            ground_truth = input_label,
                                                            run_ind = temp.run_ind,
                                                            NFOLDS = 5,
                                                            type_measure = "auc")
  )
  temp.error_count = 0  
  while (inherits(res, "try-error") && temp.error_count < 5) {
    temp.error_count = temp.error_count + 1
    res <- try(   
      cp.snf.single.l1000[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = feature.l1000, 
                                                                parameters = list(K = snf.parameter), 
                                                                data = input_data, 
                                                                partition = input_partition$cp,
                                                                ground_truth = input_label,
                                                                run_ind = temp.run_ind,
                                                                NFOLDS = 5,
                                                                type_measure = "auc")
    )
  }
}


print("completed!")
print(paste("BEGIN and END:", PARTITION_BEGIN, PARTITION_END))
setwd("/home/zhaoche7/")
rm(input_data)
if (exists("training_var_amount")) {
  save.image(paste0("epirubicin_cp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))  
} else {
  save.image(paste0("epirubicin_cp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], ".RData"))  
}
