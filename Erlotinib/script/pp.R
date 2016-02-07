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
args[2] = 1
PARTITION_BEGIN = as.integer(args[1])
PARTITION_END = as.integer(args[2]) 
registerDoParallel(4)
load("../WS/erlotinib_data.RData")
args[3] = 10 
### End ###

feature.l1000 <- feature.l1000$pp
INPUT_NFOLDS = 5
input_data <- scale(erlotinib$patient)
input_label <- erlotinib.labels$patient
input_partition <- list()

if (length(args) == 3) {
  training_var_amount <- as.integer(args[[3]])    
  if (!is.na(training_var_amount)) {
    input_partition$pp = partition_var$pp[[training_var_amount]]
  }  else {
    stop("arg 3 is wrong")
  }
  print(paste0("erlotinib_pp_", PARTITION_BEGIN, "to", PARTITION_END, "_var", training_var_amount, ".RData"))
} else {
  print("args length != 3")
}

rm(erlotinib)

print(paste("begin", PARTITION_BEGIN, "end", PARTITION_END))
setwd("../../Common/")
# Other models make predictions ---------------------------------
source('Other_Model_Predict.R')
source('SNF_Single_Predict.R')
source('SNF_LP.R')

snf.parameter <- seq(from = 2, to = 15, by = 3)

# mRMR + SNF ---------------------------------
source('mRMR_getFeatures.R')
pp.snf.single.mRMR1000 = list()
pp.other_model.mRMR1000 = list()

for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {  
  print(Sys.time())
  print(temp.run_ind)
  gc()
  temp.mRMR_features <- mRMR_getFeatures(
    input_data[input_partition$pp[[temp.run_ind]]$training_index.single, ], 
    as.ordered(input_label[input_partition$pp[[temp.run_ind]]$training_index.single]), 
    feature_count = 1000, 
    solution_count = 8)
  print(paste("Number of mRMR features:", length(temp.mRMR_features)))
  temp.mRMR_features <- as.integer(temp.mRMR_features)
  gc()
  
  pp.snf.single.mRMR1000[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = temp.mRMR_features, 
                                                               parameters = list(K = snf.parameter), 
                                                               data = input_data, 
                                                               partition = input_partition$pp,
                                                               ground_truth = input_label,
                                                               run_ind = temp.run_ind,
                                                               NFOLDS = INPUT_NFOLDS,
                                                               type_measure = "acc")
  pp.snf.single.mRMR1000[[temp.run_ind]]$feature.sets = temp.mRMR_features
  
  pp.other_model.mRMR1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                                 ground_truth = input_label, 
                                                                 partition = input_partition$pp,
                                                                 selected_features = temp.mRMR_features, 
                                                                 NFOLDS = INPUT_NFOLDS, 
                                                                 N_CV_REPEATS = 1, 
                                                                 run_ind = temp.run_ind, 
                                                                 type_measure = "acc")
}

pp.other_model.all = list()
for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {
  print(Sys.time())
  print(temp.run_ind)
  
  pp.other_model.all[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                            ground_truth = input_label, 
                                                            partition = input_partition$pp,
                                                            selected_features = NULL, 
                                                            NFOLDS = INPUT_NFOLDS, 
                                                            N_CV_REPEATS = 1, 
                                                            run_ind = temp.run_ind,
                                                            type_measure = "acc")
}

pp.other_model.l1000 = list()
for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {
  print(Sys.time())
  print(temp.run_ind)
  res <- try( 
    pp.other_model.l1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                                ground_truth = input_label, 
                                                                partition = input_partition$pp,
                                                                selected_features = feature.l1000, 
                                                                NFOLDS = INPUT_NFOLDS, 
                                                                N_CV_REPEATS = 1, 
                                                                run_ind = temp.run_ind, 
                                                                type_measure = "acc")
  )
  temp.error_count = 0
  while (inherits(res, "try-error") && temp.error_count < 5) {
    temp.error_count = temp.error_count + 1
    res <- try( 
      pp.other_model.l1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                                  ground_truth = input_label, 
                                                                  partition = input_partition$pp,
                                                                  selected_features = feature.l1000, 
                                                                  NFOLDS = INPUT_NFOLDS, 
                                                                  N_CV_REPEATS = 1, 
                                                                  run_ind = temp.run_ind, 
                                                                  type_measure = "acc")
    )
  }
}

pp.snf.single.all = list()

for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {  
  print(Sys.time())
  print(temp.run_ind)
  
  pp.snf.single.all[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = NULL, 
                                                          parameters = list(K = snf.parameter), 
                                                          data = input_data, 
                                                          partition = input_partition$pp,
                                                          ground_truth = input_label,
                                                          run_ind = temp.run_ind,
                                                          NFOLDS = INPUT_NFOLDS, 
                                                          type_measure = "acc")
  temp.error_count = 0  
}

pp.snf.single.l1000 = list()
for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {  
  print(Sys.time())
  print(temp.run_ind)
  
  pp.snf.single.l1000[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = feature.l1000, 
                                                            parameters = list(K = snf.parameter), 
                                                            data = input_data, 
                                                            partition = input_partition$pp,
                                                            ground_truth = input_label,
                                                            run_ind = temp.run_ind,
                                                            NFOLDS = INPUT_NFOLDS,
                                                            type_measure = "acc")
}

print("completed!")
print(paste("BEGIN and END:", PARTITION_BEGIN, PARTITION_END))
setwd("../Erlotinib/output_WS/")
rm(input_data)

if (exists("training_var_amount")) {
  save.image(paste0("erlotinib_pp_", PARTITION_BEGIN, "to", PARTITION_END, "_var", training_var_amount, ".RData"))  
} else {
  save.image(paste0("erlotinib_pp_", PARTITION_BEGIN, "to", PARTITION_END, ".RData"))  
}


