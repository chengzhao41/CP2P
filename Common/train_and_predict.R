stopifnot(!is.null(input_feature.l1000))
stopifnot(mean(input_data) < 0.1)
stopifnot(dim(input_data)[1] == length(input_label))
stopifnot(rownames(input_data) == names(input_label))

source('Common/Other_Model_Predict.R')
source('Common/mRMR_getFeatures.R')
source('Common/SNF_Single_Predict.R')
source('Common/SNF_LP.R')

print("Training and Predicting!")

# mRMR + SNF ---------------------------------
snf.single.mRMR1000 = list()
other_model.mRMR1000 = list()

for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {  
  print(Sys.time())
  print(temp.run_ind)
  
  temp.mRMR_features <- mRMR_getFeatures(
    input_data[input_partition[[temp.run_ind]]$training_index, ], 
    as.ordered(input_label[input_partition[[temp.run_ind]]$training_index]), 
    feature_count = 1000, 
    solution_count = 1)
  print(paste("Number of mRMR features:", length(temp.mRMR_features)))
  temp.mRMR_features <- as.integer(temp.mRMR_features)
  
  snf.single.mRMR1000[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = temp.mRMR_features, 
                                                               parameters = list(K = input_snf.parameter), 
                                                               data = input_data, 
                                                               partition = input_partition,
                                                               ground_truth = input_label,
                                                               run_ind = temp.run_ind,
                                                               NFOLDS = 5,
                                                               type_measure = input.type_measure)
  
  snf.single.mRMR1000[[temp.run_ind]]$feature.sets = temp.mRMR_features
  
  other_model.mRMR1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                                 ground_truth = input_label, 
                                                                 partition = input_partition,
                                                                 selected_features = temp.mRMR_features, 
                                                                 NFOLDS = 5, 
                                                                 N_CV_REPEATS = 1, 
                                                                 run_ind = temp.run_ind, 
                                                                 type_measure = input.type_measure)
}

# Other models make predictions ---------------------------------
other_model.all = list()
for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {
  print(Sys.time())
  print(temp.run_ind)
  
  other_model.all[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                            ground_truth = input_label, 
                                                            partition = input_partition,
                                                            selected_features = NULL, 
                                                            NFOLDS = 5, 
                                                            N_CV_REPEATS = 1, 
                                                            run_ind = temp.run_ind,
                                                            type_measure = input.type_measure)
}

other_model.l1000 = list()
for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {
  print(Sys.time())
  print(temp.run_ind)
  other_model.l1000[[temp.run_ind]] <- Other_Model_Predict(data = input_data, 
                                                              ground_truth = input_label, 
                                                              partition = input_partition,
                                                              selected_features = input_feature.l1000, 
                                                              NFOLDS = 5, 
                                                              N_CV_REPEATS = 1, 
                                                              run_ind = temp.run_ind, 
                                                              type_measure = input.type_measure)
}

snf.single.all = list()

for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {  
  print(Sys.time())
  print(temp.run_ind)
  
  snf.single.all[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = NULL, 
                                                          parameters = list(K = input_snf.parameter), 
                                                          data = input_data, 
                                                          partition = input_partition,
                                                          ground_truth = input_label,
                                                          run_ind = temp.run_ind,
                                                          NFOLDS = 5, 
                                                          type_measure = input.type_measure)
}

snf.single.l1000 = list()
for (temp.run_ind in PARTITION_BEGIN:PARTITION_END) {  
  print(Sys.time())
  print(temp.run_ind)
  
  snf.single.l1000[[temp.run_ind]] <- SNF_Single_Predict(feature.sets = input_feature.l1000, 
                                                            parameters = list(K = input_snf.parameter), 
                                                            data = input_data, 
                                                            partition = input_partition,
                                                            ground_truth = input_label,
                                                            run_ind = temp.run_ind,
                                                            NFOLDS = 5,
                                                            type_measure = input.type_measure)
}

print("completed!")
print(paste("BEGIN and END:", PARTITION_BEGIN, PARTITION_END))
rm(input_data)
