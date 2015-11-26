stop("Set working directory to current source file")
#setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Code/Docetaxel/Script")

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
load("../WS/docetaxel_data.RData")
args[3] = "slope"
args[4] = 10 
### End ###

INPUT_NFOLDS = 5
if (args[3] == "slope") {
  input_data <- docetaxel$combined_slope.ComBat
  input_label <- docetaxel.labels$slope_combined
  input_partition <- partition$slope
} else if (args[3] == "auc") {
  input_data <- docetaxel$combined_AUC.ComBat
  input_label <- docetaxel.labels$AUC_combined
  input_partition <- partition$AUC
} else if (args[3] == "ic50") {
  input_data <- docetaxel$combined_IC50.ComBat
  input_label <- docetaxel.labels$IC50_combined
  input_partition <- partition$IC50
} else if (args[3] == "slope_breast") {
  input_data <- docetaxel$breast_slope.ComBat
  input_label <- docetaxel.labels$slope_breast
  input_partition <- partition$slope_breast
  INPUT_NFOLDS = 3
} else if (args[3] == "auc_breast") {
  input_data <- docetaxel$breast_AUC.ComBat
  input_label <- docetaxel.labels$AUC_breast
  input_partition <- partition$AUC_breast
  INPUT_NFOLDS = 3
} else if (args[3] == "ic50_breast") {
  input_data <- docetaxel$breast_IC50.ComBat
  input_label <- docetaxel.labels$IC50_breast
  input_partition <- partition$IC50_breast
  INPUT_NFOLDS = 3
} else {
  stopifnot(FALSE)
}

snf.parameter <- seq(from = 5, to = 30, by = 5)

if (length(args) == 4) {
  training_var_amount <- as.integer(args[[4]])
  print(paste0("docetaxel_cpp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))  
  if (args[3] == "slope") {    
    input_partition$cpp <- partition_var$cpp[[training_var_amount]]$slope
  } else if (args[3] == "auc") {    
    input_partition$cpp <- partition_var$cpp[[training_var_amount]]$AUC
  } else if (args[3] == "ic50") {    
    input_partition$cpp <- partition_var$cpp[[training_var_amount]]$IC50
  } else if (args[3] == "slope_breast") {    
    input_partition$cpp <- partition_var$cpp_breast_only[[training_var_amount]]$slope
    temp.max <- length(partition_var$cpp_breast_only[[1]]$IC50[[1]]$training_index.single)
    snf.parameter <- seq(from = 5, to = temp.max, by = 5)
  } else if (args[3] == "auc_breast") {    
    input_partition$cpp <- partition_var$cpp_breast_only[[training_var_amount]]$AUC
    temp.max <- length(partition_var$cpp_breast_only[[1]]$AUC[[1]]$training_index.single)
    snf.parameter <- seq(from = 5, to = temp.max, by = 5)
  } else if (args[3] == "ic50_breast") {    
    input_partition$cpp <- partition_var$cpp_breast_only[[training_var_amount]]$IC50
    temp.max <- length(partition_var$cpp_breast_only[[1]]$slope[[1]]$training_index.single)
    snf.parameter <- seq(from = 5, to = temp.max, by = 5)
  } else {
    stop("arg 3 is wrong")
  }  
}

feature.l1000 <- feature.l1000$cp

setwd("../../Common/")

# Other models make predictions ---------------------------------
source('Other_Model_Predict.R')

#1: original run
cpp.other_model.all = list()
for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {
  print(Sys.time())
  print(temp.run_ind)
  cpp.other_model.all[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                             ground_truth = input_label, 
                                                             partition = input_partition$cpp,
                                                             selected_features = NULL, 
                                                             NFOLDS = 5, 
                                                             N_CV_REPEATS = 1, 
                                                             run_ind = temp.run_ind,
                                                             type_measure = "acc")
}

cpp.other_model.l1000 = list()
for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {
  print(Sys.time())
  print(temp.run_ind)
  cpp.other_model.l1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                               ground_truth = input_label, 
                                                               partition = input_partition$cpp,
                                                               selected_features = feature.l1000, 
                                                               NFOLDS = 5, 
                                                               N_CV_REPEATS = 1, 
                                                               run_ind = temp.run_ind, 
                                                               type_measure = "acc"
  )
}

source('SNF_Single_Predict.R')
source('SNF_LP.R')
cpp.snf.single.all = list()

for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {  
  print(Sys.time())
  print(temp.run_ind)
  cpp.snf.single.all[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = NULL, 
                                                           parameters = list(K = snf.parameter), 
                                                           data = input_data, 
                                                           partition = input_partition$cpp,
                                                           ground_truth = input_label,
                                                           run_ind = temp.run_ind,
                                                           NFOLDS = 5, 
                                                           type_measure = "acc"                                                            
  )
}

cpp.snf.single.l1000 = list()
for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {  
  print(Sys.time())
  print(temp.run_ind)
  cpp.snf.single.l1000[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = feature.l1000, 
                                                             parameters = list(K = snf.parameter), 
                                                             data = input_data, 
                                                             partition = input_partition$cpp,
                                                             ground_truth = input_label,
                                                             run_ind = temp.run_ind,
                                                             NFOLDS = 5,
                                                             type_measure = "acc" 
  )
}

# mRMR + SNF ---------------------------------
source('mRMR_getFeatures.R')
cpp.snf.single.mRMR1000 = list()
cpp.other_model.mRMR1000 = list()

for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {  
  print(Sys.time())
  print(temp.run_ind)
  
  temp.mRMR_features <- mRMR_getFeatures(
    input_data[input_partition$cpp[[temp.run_ind]]$training_index.single, ], 
    as.ordered(input_label[input_partition$cpp[[temp.run_ind]]$training_index.single]), 
    feature_count = 1000, 
    solution_count = 1)
  print(paste("Number of mRMR features:", length(temp.mRMR_features)))
  temp.mRMR_features <- as.integer(temp.mRMR_features)
  
  cpp.snf.single.mRMR1000[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = temp.mRMR_features, 
                                                                parameters = list(K = snf.parameter), 
                                                                data = input_data, 
                                                                partition = input_partition$cpp,
                                                                ground_truth = input_label,
                                                                run_ind = temp.run_ind,
                                                                NFOLDS = 5,
                                                                type_measure = "acc")
  cpp.snf.single.mRMR1000[[temp.run_ind]]$feature.sets = temp.mRMR_features
  
  cpp.other_model.mRMR1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                                  ground_truth = input_label, 
                                                                  partition = input_partition$cpp,
                                                                  selected_features = temp.mRMR_features, 
                                                                  NFOLDS = 5, 
                                                                  N_CV_REPEATS = 1, 
                                                                  run_ind = temp.run_ind, 
                                                                  type_measure = "acc")
}

print("completed!")
print(paste("BEGIN and END:", PARTITION_BEGIN, PARTITION_END))
setwd("../Docetaxel/output_WS/")
remove(docetaxel)
remove(input_data)


if (exists("training_var_amount")) {
  save.image(paste0("docetaxel_cpp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))  
} else {
  save.image(paste0("docetaxel_cpp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], ".RData"))  
}
