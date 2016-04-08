stopifnot(!is.null(input_feature.l1000))
stopifnot(mean(input_data) < 0.1)
stopifnot(dim(input_data)[1] == length(input_label))
stopifnot(rownames(input_data) == names(input_label))

source('Common/other_model_train_only.R')
source('Common/other_model_predict_only.R.R')
source('Common/mRMR_getFeatures.R')
source('Common/SNF_Single_Predict.R')
source('Common/SNF_LP.R')

print("Training and Predicting!")

# mRMR + SNF ---------------------------------
snf.single.mRMR1000 = list()
other_model.mRMR1000 = list()

already_run_mRMR = FALSE

for (loop_ind.parInd in PARTITION_BEGIN:PARTITION_END) {  
  print(Sys.time())
  print(loop_ind.parInd)
  gc()
  if (train_once == FALSE || already_run_mRMR == FALSE) {
    temp.mRMR_features <- mRMR_getFeatures(
      input_data[input_partition[[loop_ind.parInd]]$training_index, ], 
      as.ordered(input_label[input_partition[[loop_ind.parInd]]$training_index]), 
      feature_count = 1000, 
      solution_count = 1)
    print(paste("Number of mRMR features:", length(temp.mRMR_features)))
    temp.mRMR_features <- as.integer(temp.mRMR_features)
    already_run_mRMR = TRUE
    gc()
    
    input.models <- Other_Model_Train_Only(data = input_data, 
                                           ground_truth = input_label, 
                                           partition = input_partition,
                                           selected_features = temp.mRMR_features, 
                                           NFOLDS = INPUT.NFOLDS, 
                                           N_CV_REPEATS = 1, 
                                           run_ind = loop_ind.parInd, 
                                           type_measure = input.type_measure)
  }
  snf.single.mRMR1000[[loop_ind.parInd]] <- SNF_Single_Predict(feature.sets = temp.mRMR_features, 
                                                               parameters = list(K = input_snf.parameter), 
                                                               data = input_data, 
                                                               partition = input_partition,
                                                               ground_truth = input_label,
                                                               run_ind = loop_ind.parInd,
                                                               NFOLDS = INPUT.NFOLDS,
                                                               type_measure = input.type_measure)
  
  snf.single.mRMR1000[[loop_ind.parInd]]$feature.sets = temp.mRMR_features
  
  other_model.mRMR1000[[loop_ind.parInd]] <- Other_Model_Predict_Only(data = input_data, 
                                                                 ground_truth = input_label, 
                                                                 partition = input_partition,
                                                                 selected_features = temp.mRMR_features, 
                                                                 NFOLDS = INPUT.NFOLDS, 
                                                                 N_CV_REPEATS = 1, 
                                                                 run_ind = loop_ind.parInd, 
                                                                 type_measure = input.type_measure,
                                                                 input.models = input.models)
}

# Other models make predictions ---------------------------------
other_model.all = list()
already_trained = FALSE
for (loop_ind.parInd in PARTITION_BEGIN:PARTITION_END) {
  print(Sys.time())
  print(loop_ind.parInd)
  gc()
  
  if (train_once == FALSE || already_trained == FALSE) {
    input.models <- Other_Model_Train_Only(data = input_data, 
                                           ground_truth = input_label, 
                                           partition = input_partition,
                                           selected_features = NULL, 
                                           NFOLDS = INPUT.NFOLDS, 
                                           N_CV_REPEATS = 1, 
                                           run_ind = loop_ind.parInd,
                                           type_measure = input.type_measure)
    already_trained = TRUE
  }
  
  
  other_model.all[[loop_ind.parInd]] <- Other_Model_Predict_Only(data = input_data, 
                                                            ground_truth = input_label, 
                                                            partition = input_partition,
                                                            selected_features = NULL, 
                                                            NFOLDS = INPUT.NFOLDS, 
                                                            N_CV_REPEATS = 1, 
                                                            run_ind = loop_ind.parInd,
                                                            type_measure = input.type_measure,
                                                            input.models = input.models)
}

other_model.l1000 = list()
already_trained = FALSE
for (loop_ind.parInd in PARTITION_BEGIN:PARTITION_END) {
  print(Sys.time())
  print(loop_ind.parInd)
  
  if (train_once == FALSE || already_trained == FALSE) {
    input.models <- Other_Model_Train_Only(data = input_data, 
                                           ground_truth = input_label, 
                                           partition = input_partition,
                                           selected_features = input_feature.l1000, 
                                           NFOLDS = INPUT.NFOLDS, 
                                           N_CV_REPEATS = 1, 
                                           run_ind = loop_ind.parInd, 
                                           type_measure = input.type_measure)
    already_trained = TRUE
  }
  
  other_model.l1000[[loop_ind.parInd]] <- Other_Model_Predict_Only(data = input_data, 
                                                              ground_truth = input_label, 
                                                              partition = input_partition,
                                                              selected_features = input_feature.l1000, 
                                                              NFOLDS = INPUT.NFOLDS, 
                                                              N_CV_REPEATS = 1, 
                                                              run_ind = loop_ind.parInd, 
                                                              type_measure = input.type_measure,
                                                              input.models = input.models)
}

snf.single.all = list()

for (loop_ind.parInd in PARTITION_BEGIN:PARTITION_END) {  
  print(Sys.time())
  print(loop_ind.parInd)
  
  snf.single.all[[loop_ind.parInd]] <- SNF_Single_Predict(feature.sets = NULL, 
                                                          parameters = list(K = input_snf.parameter), 
                                                          data = input_data, 
                                                          partition = input_partition,
                                                          ground_truth = input_label,
                                                          run_ind = loop_ind.parInd,
                                                          NFOLDS = INPUT.NFOLDS, 
                                                          type_measure = input.type_measure)
}

snf.single.l1000 = list()
for (loop_ind.parInd in PARTITION_BEGIN:PARTITION_END) {  
  print(Sys.time())
  print(loop_ind.parInd)
  
  snf.single.l1000[[loop_ind.parInd]] <- SNF_Single_Predict(feature.sets = input_feature.l1000, 
                                                            parameters = list(K = input_snf.parameter), 
                                                            data = input_data, 
                                                            partition = input_partition,
                                                            ground_truth = input_label,
                                                            run_ind = loop_ind.parInd,
                                                            NFOLDS = INPUT.NFOLDS,
                                                            type_measure = input.type_measure)
}

print("completed!")
print(paste("BEGIN and END:", PARTITION_BEGIN, PARTITION_END))
rm(input_data)
