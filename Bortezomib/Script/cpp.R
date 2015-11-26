stop("Set working directory to current source file")
#setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Code/Bortezomib/Script")

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
load("../WS/bortezomib_data.RData")
args[3] = "slope"
args[4] = 2 
### End ###

print(paste("begin", PARTITION_BEGIN, "end", PARTITION_END))
setwd("../../Common/")

if (args[3] == "slope") {
  input_data <- bortezomib$combined_slope.sva
  input_label <- bortezomib.labels$slope_combined
  input_partition <- partition$slope
} else if (args[3] == "auc") {
  input_data <- bortezomib$combined_AUC.sva
  input_label <- bortezomib.labels$AUC_combined
  input_partition <- partition$AUC
} else if (args[3] == "ic50") {
  input_data <- bortezomib$combined_IC50.sva
  input_label <- bortezomib.labels$IC50_combined
  input_partition <- partition$IC50
} else {
  stopifnot(FALSE)
}

if (length(args) == 4) {
  training_var_amount <- as.integer(args[[4]])
  print(paste0("bortezomib_cpp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))  
  if (args[3] == "slope") {    
    input_partition$cpp <- partition_var$cpp[[training_var_amount]]$slope
  } else if (args[3] == "auc") {    
    input_partition$cpp <- partition_var$cpp[[training_var_amount]]$AUC
  } else if (args[3] == "ic50") {    
    input_partition$cpp <- partition_var$cpp[[training_var_amount]]$IC50
  } else {
    stop("arg 3 is wrong")
  }  
}

feature.l1000 <- feature.l1000$cp

# Other models make predictions ---------------------------------
source('Other_Model_Predict.R')

snf.parameter <- seq(from = 5, to = 30, by = 5)

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
                                                             type_measure = "auc")
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
                                                               type_measure = "auc")
  
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
                                                           type_measure = "auc")
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
                                                             type_measure = "auc")
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
                                                                type_measure = "auc")
  cpp.snf.single.mRMR1000[[temp.run_ind]]$feature.sets = temp.mRMR_features
  
  cpp.other_model.mRMR1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                                  ground_truth = input_label, 
                                                                  partition = input_partition$cpp,
                                                                  selected_features = temp.mRMR_features, 
                                                                  NFOLDS = 5, 
                                                                  N_CV_REPEATS = 1, 
                                                                  run_ind = temp.run_ind, 
                                                                  type_measure = "auc")
  
}

print("completed!")
print(paste("BEGIN and END:", PARTITION_BEGIN, PARTITION_END))
setwd("../Bortezomib/output_WS/")
rm(bortezomib, input_data)

if (exists("training_var_amount")) {
  save.image(paste0("bortezomib_cpp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))
} else {  
  save.image(paste0("bortezomib_cpp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], ".RData"))
}

