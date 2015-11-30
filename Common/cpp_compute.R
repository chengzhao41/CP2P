snf.parameter <- seq(from = 5, to = 30, by = 5)
# Other models make predictions ---------------------------------
source('Common/Other_Model_Predict.R')
source('Common/mRMR_getFeatures.R')
source('Common/SNF_Single_Predict.R')
source('Common/SNF_LP.R')

# mRMR + SNF ---------------------------------
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

print(paste("BEGIN and END:", PARTITION_BEGIN, PARTITION_END))
rm(input_data)
print("completed computation!")
