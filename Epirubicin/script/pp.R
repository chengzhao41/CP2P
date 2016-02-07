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
  input_data <- scale(epirubicin$patient)
  input_label <- epirubicin.labels$patient 
  input_partition <- partition$slope$pp
}

if (length(args) == 3) {
  training_var_amount <- as.integer(args[[3]])  
  if (!is.na(training_var_amount)) {
    input_partition = partition_var$pp[[training_var_amount]]
  }  else {
    stop("arg 3 is wrong")
  }
  print(paste0("epirubicin_pp_", PARTITION_BEGIN, "to", PARTITION_END, "_var", training_var_amount, ".RData"))
} else {
  stop("args length != 3")
}

rm(epirubicin)
# mRMR + SNF ---------------------------------
source('mRMR_getFeatures.R')
source('Other_Model_Predict.R')
source('SNF_Single_Predict.R')
source('SNF_LP_Cheng.R')
snf.parameter <- seq(from = 5, to = 30, by = 5)

pp.snf.single.mRMR1000 = list()
pp.other_model.mRMR1000 = list()

for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {  
  gc()
  print(Sys.time())
  print(temp.run_ind)
  
  temp.mRMR_features <- mRMR_getFeatures(
    input_data[input_partition[[temp.run_ind]]$training_index.single, ], 
    as.ordered(input_label[input_partition[[temp.run_ind]]$training_index.single]), 
    feature_count = 1000, 
    solution_count = 1)
  print(paste("Number of mRMR features:", length(temp.mRMR_features)))
  temp.mRMR_features <- as.integer(temp.mRMR_features)
  gc()
  res <- try(  
    pp.snf.single.mRMR1000[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = temp.mRMR_features, 
                                                                 parameters = list(K = snf.parameter), 
                                                                 data = input_data, 
                                                                 partition = input_partition,
                                                                 ground_truth = input_label,
                                                                 run_ind = temp.run_ind,
                                                                 NFOLDS = 5,
                                                                 type_measure = "auc")
  )
  temp.error_count = 0  
  while (inherits(res, "try-error") && temp.error_count < 5) {
    temp.error_count = temp.error_count + 1    
    res <- try(  
      pp.snf.single.mRMR1000[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = temp.mRMR_features, 
                                                                   parameters = list(K = snf.parameter), 
                                                                   data = input_data, 
                                                                   partition = input_partition,
                                                                   ground_truth = input_label,
                                                                   run_ind = temp.run_ind,
                                                                   NFOLDS = 5,
                                                                   type_measure = "auc")
    )
  }
  pp.snf.single.mRMR1000[[temp.run_ind]]$feature.sets = temp.mRMR_features
  
  res <- try(    
    pp.other_model.mRMR1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                                   ground_truth = input_label, 
                                                                   partition = input_partition,
                                                                   selected_features = temp.mRMR_features, 
                                                                   NFOLDS = 5, 
                                                                   N_CV_REPEATS = 1, 
                                                                   run_ind = temp.run_ind, 
                                                                   type_measure = "auc")
  )
  temp.error_count = 0  
  while (inherits(res, "try-error") && temp.error_count < 5) {
    res <- try(    
      pp.other_model.mRMR1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                                     ground_truth = input_label, 
                                                                     partition = input_partition,
                                                                     selected_features = temp.mRMR_features, 
                                                                     NFOLDS = 5, 
                                                                     N_CV_REPEATS = 1, 
                                                                     run_ind = temp.run_ind, 
                                                                     type_measure = "auc")
    )   
  }
}

# Other models make predictions ---------------------------------
pp.other_model.all = list()
for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {
  gc()
  print(Sys.time())
  print(temp.run_ind)
  
  res <- try( 
  pp.other_model.all[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                            ground_truth = input_label, 
                                                            partition = input_partition,
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
      pp.other_model.all[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                                ground_truth = input_label, 
                                                                partition = input_partition,
                                                                selected_features = NULL, 
                                                                NFOLDS = 5, 
                                                                N_CV_REPEATS = 1, 
                                                                run_ind = temp.run_ind,
                                                                type_measure = "auc")
    )
  }
}

pp.other_model.l1000 = list()
for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {
  print(Sys.time())
  print(temp.run_ind)
  res <- try( 
  pp.other_model.l1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                              ground_truth = input_label, 
                                                              partition = input_partition,
                                                              selected_features = feature.l1000$pp, 
                                                              NFOLDS = 5, 
                                                              N_CV_REPEATS = 1, 
                                                              run_ind = temp.run_ind, 
                                                              type_measure = "auc")
  )
  temp.error_count = 0
  while (inherits(res, "try-error") && temp.error_count < 5) {
    temp.error_count = temp.error_count + 1
    res <- try( 
      pp.other_model.l1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                                  ground_truth = input_label, 
                                                                  partition = input_partition,
                                                                  selected_features = feature.l1000$pp, 
                                                                  NFOLDS = 5, 
                                                                  N_CV_REPEATS = 1, 
                                                                  run_ind = temp.run_ind, 
                                                                  type_measure = "auc")
    )
  }
}

pp.snf.single.all = list()

for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {  
  print(Sys.time())
  print(temp.run_ind)
  res <- try( 
  pp.snf.single.all[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = NULL, 
                                                          parameters = list(K = snf.parameter), 
                                                          data = input_data, 
                                                          partition = input_partition,
                                                          ground_truth = input_label,
                                                          run_ind = temp.run_ind,
                                                          NFOLDS = 5, 
                                                          type_measure = "auc")
  )
  temp.error_count = 0  
  while (inherits(res, "try-error") && temp.error_count < 5) {
    temp.error_count = temp.error_count + 1
    res <- try( 
      pp.snf.single.all[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = NULL, 
                                                              parameters = list(K = snf.parameter), 
                                                              data = input_data, 
                                                              partition = input_partition,
                                                              ground_truth = input_label,
                                                              run_ind = temp.run_ind,
                                                              NFOLDS = 5, 
                                                              type_measure = "auc")      
    )
  }
}

pp.snf.single.l1000 = list()
for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {  
  print(Sys.time())
  print(temp.run_ind)
  
  res <- try(   
  pp.snf.single.l1000[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = feature.l1000$pp, 
                                                            parameters = list(K = snf.parameter), 
                                                            data = input_data, 
                                                            partition = input_partition,
                                                            ground_truth = input_label,
                                                            run_ind = temp.run_ind,
                                                            NFOLDS = 5,
                                                            type_measure = "auc")
  )
  temp.error_count = 0  
  while (inherits(res, "try-error") && temp.error_count < 5) {
    temp.error_count = temp.error_count + 1
    res <- try(   
      pp.snf.single.l1000[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = feature.l1000$pp, 
                                                                parameters = list(K = snf.parameter), 
                                                                data = input_data, 
                                                                partition = input_partition,
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
  save.image(paste0("epirubicin_pp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))  
} else {
  save.image(paste0("epirubicin_pp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], ".RData"))  
}
