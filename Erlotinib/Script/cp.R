stop("Set working directory to current source file")
#setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Code/Erlotinib/Script")

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

### User supplies these values ###
args <- vector()

args[1] = 1
args[2] = 2
PARTITION_BEGIN = as.integer(args[1])
PARTITION_END = as.integer(args[2]) 
registerDoParallel(4)
load("../WS/erlotinib_data.RData")
args[3] = "slope"
args[4] = 3
### End ###

if (args[3] == "slope") {
  input_data <- erlotinib$all.ComBat.slope
  input_label <- erlotinib.labels$all.slope
  input_partition <- partition$all.slope
  feature.l1000 <- feature.l1000$cp
} else if ((args[3] == "ic50") || (args[3] == "auc")) {
  input_data <- erlotinib$all.ComBat.IC50
  input_label <- erlotinib.labels$all.IC50
  input_partition <- partition$all.IC50
  feature.l1000 <- feature.l1000$cp
} else if (args[3] == "lung_slope") {
  input_data <- erlotinib$lung_all.slope.ComBat
  input_label <- erlotinib.labels$lung_all.slope
  #input_partition <- partition$lung.slope
  input_partition <- list()
  input_partition$cp <- list(list(test_index = which(erlotinib.labels$lung_all.slope.source == "patient"), 
                                  training_index.single = which(erlotinib.labels$lung_all.slope.source != "patient")))
  feature.l1000 <- feature.l1000$cp
} else if (args[3] == "lung_ic50") {
  input_data <- erlotinib$lung_all.IC50.ComBat
  input_label <- erlotinib.labels$lung_all.IC50
  #input_partition <- partition$lung.IC50
  input_partition <- list()
  input_partition$cp <- list(list(test_index = which(erlotinib.labels$lung_all.IC50.source == "patient"), 
                                  training_index.single = which(erlotinib.labels$lung_all.IC50.source != "patient")))
  feature.l1000 <- feature.l1000$cp
} else if (args[3] == "GEO") {
  input_data <- erlotinib$combined_GEO
  input_label <- erlotinib.labels$GEO_IC50_combined
  input_partition <- partition$GEO
  feature.l1000 <- feature.l1000$geo_cp  
} else {
  stop("false args[3]")
  input_partition <- list()
  input_partition$cp <- list(list(test_index = which(erlotinib.labels$IC50_breast.source == "patient"), 
                                  training_index.single = which(erlotinib.labels$IC50_breast.source != "patient")))
}

if (length(args) == 4) {
  training_var_amount <- as.integer(args[[4]])
  if (args[3] == "slope") {    
    input_partition$cp <- partition_var$cp[[training_var_amount]]$slope
  } else if ((args[3] == "ic50") || (args[3] == "auc")) {    
    input_partition$cp <- partition_var$cp[[training_var_amount]]$IC50
  } else if (args[3] == "lung_slope") {    
    input_partition$cp <- partition_var$cp_lung_only[[training_var_amount]]$slope
  } else if ((args[3] == "lung_ic50") || (args[3] == "auc")) {    
    input_partition$cp <- partition_var$cp_lung_only[[training_var_amount]]$IC50
  } else {
    stop("arg 3 is wrong")
  }
  print(paste0("erlotinib_cp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))
} 

rm(erlotinib)

setwd("../../Common/")

# Other models make predictions ---------------------------------
source('Other_Model_Predict.R')
source('SNF_Single_Predict.R')
source('SNF_LP.R')

snf.parameter <- seq(from = 5, to = 30, by = 5)

# mRMR + SNF ---------------------------------
source('mRMR_getFeatures.R')
cp.snf.single.mRMR1000 = list()
cp.other_model.mRMR1000 = list()

for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {  
  print(Sys.time())
  print(temp.run_ind)
  
  temp.mRMR_features <- mRMR_getFeatures(
    input_data[input_partition$cp[[temp.run_ind]]$training_index.single, ], 
    as.ordered(input_label[input_partition$cp[[temp.run_ind]]$training_index.single]), 
    feature_count = 1000, 
    solution_count = 1)
  print(paste("Number of mRMR features:", length(temp.mRMR_features)))
  temp.mRMR_features <- as.integer(temp.mRMR_features)
  
  cp.snf.single.mRMR1000[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = temp.mRMR_features, 
                                                               parameters = list(K = snf.parameter), 
                                                               data = input_data, 
                                                               partition = input_partition$cp,
                                                               ground_truth = input_label,
                                                               run_ind = temp.run_ind,
                                                               NFOLDS = 5,
                                                               type_measure = "auc")
  
  cp.snf.single.mRMR1000[[temp.run_ind]]$feature.sets = temp.mRMR_features
  
  cp.other_model.mRMR1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                                 ground_truth = input_label, 
                                                                 partition = input_partition$cp,
                                                                 selected_features = temp.mRMR_features, 
                                                                 NFOLDS = 5, 
                                                                 N_CV_REPEATS = 1, 
                                                                 run_ind = temp.run_ind, 
                                                                 type_measure = "auc")
  
  
}

cp.other_model.all = list()
for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {
  print(Sys.time())
  print(temp.run_ind)
  
  
  cp.other_model.all[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                            ground_truth = input_label, 
                                                            partition = input_partition$cp,
                                                            selected_features = NULL, 
                                                            NFOLDS = 5, 
                                                            N_CV_REPEATS = 1, 
                                                            run_ind = temp.run_ind,
                                                            type_measure = "auc")
}

cp.other_model.l1000 = list()
for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {
  print(Sys.time())
  print(temp.run_ind)
  
  cp.other_model.l1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                              ground_truth = input_label, 
                                                              partition = input_partition$cp,
                                                              selected_features = feature.l1000, 
                                                              NFOLDS = 3, 
                                                              N_CV_REPEATS = 1, 
                                                              run_ind = temp.run_ind, 
                                                              type_measure = "auc")
}

cp.snf.single.all = list()

for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {  
  print(Sys.time())
  print(temp.run_ind)
  
  cp.snf.single.all[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = NULL, 
                                                          parameters = list(K = snf.parameter), 
                                                          data = input_data, 
                                                          partition = input_partition$cp,
                                                          ground_truth = input_label,
                                                          run_ind = temp.run_ind,
                                                          NFOLDS = 5, 
                                                          type_measure = "auc")
}

cp.snf.single.l1000 = list()
for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {  
  print(Sys.time())
  print(temp.run_ind)
  
  cp.snf.single.l1000[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = feature.l1000, 
                                                            parameters = list(K = snf.parameter), 
                                                            data = input_data, 
                                                            partition = input_partition$cp,
                                                            ground_truth = input_label,
                                                            run_ind = temp.run_ind,
                                                            NFOLDS = 5,
                                                            type_measure = "auc")
  
}

print("completed!")
print(paste("BEGIN and END:", PARTITION_BEGIN, PARTITION_END))
setwd("../Erlotinib/output_WS/")

rm(input_data)

if (exists("training_var_amount")) {
  save.image(paste0("erlotinib_cp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))  
} else {
  save.image(paste0("erlotinib_cp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], ".RData"))  
}

