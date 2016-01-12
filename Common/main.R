rm(list = ls())
setwd("/home/zhaoche7/drp3")
source("Common/load_library.R")

### User supplies these values ###
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5)

PARTITION_BEGIN = as.integer(args[1])
PARTITION_END = as.integer(args[2])
registerDoParallel(8)

# args <- vector()
# registerDoParallel(4)
# args[1] = 1 # partition start
# args[2] = 1 # partition end
# args[3] = 10 # index for training sets
# args[4] = "bortezomib"
# args[5] = "c50p2p_slope"

### End ###

print(args)

training_var_amount <- as.integer(args[[3]])
  
if (args[4] == "bortezomib") {
  load("Bortezomib/WS/bortezomib_data.RData")
  
  if (args[5] == "p2p") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- bortezomib$patient.combat
    input_label <- bortezomib.labels$patient
    input_partition <- partition$cell_lines_all[[training_var_amount]]$p2p
    input_feature.l1000 <- feature.l1000$pp
  } else if (args[5] == "cp2p_slope") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- bortezomib$combined_slope.sva 
    input_label <- bortezomib.labels$slope_combined
    input_partition = partition$cell_lines_all[[training_var_amount]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- bortezomib$combined_AUC.sva
    input_label <- bortezomib.labels$AUC_combined
    input_partition = partition$cell_lines_all[[training_var_amount]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- bortezomib$combined_IC50.sva
    input_label <- bortezomib.labels$IC50_combined
    input_partition = partition$cell_lines_all[[training_var_amount]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50_p50") {
    stopifnot(training_var_amount <= length(partition$patient_50))
    
    input_data <- bortezomib$combined_IC50.sva
    input_label <- bortezomib.labels$IC50_combined
    input_partition = partition$patient_50[[training_var_amount]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_p50") {
    stopifnot(training_var_amount <= length(partition$patient_50))
    
    input_data <- bortezomib$combined_AUC.sva
    input_label <- bortezomib.labels$AUC_combined
    input_partition = partition$patient_50[[training_var_amount]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_slope_p50") {
    stopifnot(training_var_amount <= length(partition$patient_50))
    
    input_data <- bortezomib$combined_slope.sva
    input_label <- bortezomib.labels$slope_combined
    input_partition = partition$patient_50[[training_var_amount]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_ic50_p50") {
    stopifnot(training_var_amount <= length(partition$patient_50))
    
    input_data <- bortezomib$combined_IC50.sva
    input_label <- bortezomib.labels$IC50_combined
    input_partition = partition$patient_50[[training_var_amount]]$c2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_auc_p50") {
    stopifnot(training_var_amount <= length(partition$patient_50))
    
    input_data <- bortezomib$combined_AUC.sva
    input_label <- bortezomib.labels$AUC_combined
    input_partition = partition$patient_50[[training_var_amount]]$c2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_slope_p50") {
    stopifnot(training_var_amount <= length(partition$patient_50))
    
    input_data <- bortezomib$combined_slope.sva
    input_label <- bortezomib.labels$slope_combined
    input_partition = partition$patient_50[[training_var_amount]]$c2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_slope_p100") {
    stopifnot(training_var_amount <= length(partition$patient_100))
    
    input_data <- bortezomib$combined_slope.sva
    input_label <- bortezomib.labels$slope_combined
    input_partition = partition$patient_100[[training_var_amount]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_p100") {
    stopifnot(training_var_amount <= length(partition$patient_100))
    
    input_data <- bortezomib$combined_AUC.sva
    input_label <- bortezomib.labels$AUC_combined
    input_partition = partition$patient_100[[training_var_amount]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50_p100") {
    stopifnot(training_var_amount <= length(partition$patient_100))
    
    input_data <- bortezomib$combined_IC50.sva
    input_label <- bortezomib.labels$IC50_combined
    input_partition = partition$patient_100[[training_var_amount]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else {
    stop(paste("args[5]", args[5], "is invalid."))
  }
  
  input.type_measure = "auc"
  input_snf.parameter <- seq(from = 5, to = 30, by = 5)
  rm(bortezomib, bortezomib.labels)
} else if (args[4] == "docetaxel") {
  load("Docetaxel/WS/docetaxel_data.RData")
  
  if (args[5] == "p2p") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- docetaxel$patient.combat
    input_label <- docetaxel.labels$patient
    input_partition <- partition$cell_lines_all[[training_var_amount]]$p2p
    input_feature.l1000 <- feature.l1000$pp
  } else if (args[5] == "cp2p_slope") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- docetaxel$combined_slope.ComBat 
    input_label <- docetaxel.labels$slope_combined
    input_partition = partition$cell_lines_all[[training_var_amount]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- docetaxel$combined_AUC.ComBat
    input_label <- docetaxel.labels$AUC_combined
    input_partition = partition$cell_lines_all[[training_var_amount]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- docetaxel$combined_IC50.ComBat
    input_label <- docetaxel.labels$IC50_combined
    input_partition = partition$cell_lines_all[[training_var_amount]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_slope_breast") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- docetaxel$combined_slope.ComBat 
    input_label <- docetaxel.labels$slope_combined
    input_partition = partition$cell_lines_breast[[training_var_amount]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_breast") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- docetaxel$combined_AUC.ComBat
    input_label <- docetaxel.labels$AUC_combined
    input_partition = partition$cell_lines_breast[[training_var_amount]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50_breast") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- docetaxel$combined_IC50.ComBat
    input_label <- docetaxel.labels$IC50_combined
    input_partition = partition$cell_lines_breast[[training_var_amount]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50_p23") {
    stopifnot(training_var_amount <= length(partition$patient_23))
    
    input_data <- docetaxel$combined_IC50.ComBat
    input_label <- docetaxel.labels$IC50_combined
    input_partition = partition$patient_23[[training_var_amount]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_p23") {
    stopifnot(training_var_amount <= length(partition$patient_23))
    
    input_data <- docetaxel$combined_AUC.ComBat
    input_label <- docetaxel.labels$AUC_combined
    input_partition = partition$patient_23[[training_var_amount]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_slope_p23") {
    stopifnot(training_var_amount <= length(partition$patient_23))
    
    input_data <- docetaxel$combined_slope.ComBat
    input_label <- docetaxel.labels$slope_combined
    input_partition = partition$patient_23[[training_var_amount]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_ic50_p23") {
    stopifnot(training_var_amount <= length(partition$patient_23))
    
    input_data <- docetaxel$combined_IC50.ComBat
    input_label <- docetaxel.labels$IC50_combined
    input_partition = partition$patient_23[[training_var_amount]]$c2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_auc_p23") {
    stopifnot(training_var_amount <= length(partition$patient_23))
    
    input_data <- docetaxel$combined_AUC.ComBat
    input_label <- docetaxel.labels$AUC_combined
    input_partition = partition$patient_23[[training_var_amount]]$c2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_slope_p23") {
    stopifnot(training_var_amount <= length(partition$patient_23))
    
    input_data <- docetaxel$combined_slope.ComBat
    input_label <- docetaxel.labels$slope_combined
    input_partition = partition$patient_23[[training_var_amount]]$c2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else {
    stop(paste("args[5]", args[5], "is invalid."))
  }
  
  input.type_measure = "acc"
  input_snf.parameter <- seq(from = 5, to = 30, by = 5)
  rm(docetaxel, docetaxel.labels)
} else {
  stop(paste("args[4]", args[4], "is invalid."))
}

stopifnot(!is.null(input_data))
stopifnot(!is.null(input_label))
stopifnot(!is.null(input_partition))

rm(partition, feature.l1000)

# do the computation
print(paste0(args[4], "_", args[5], "_", PARTITION_BEGIN, "to", PARTITION_END, "_parInd", training_var_amount, ".RData"))
source("Common/train_and_predict.R")
##

save.image(paste0(args[4], "_", args[5], "_", PARTITION_BEGIN, "to", PARTITION_END, "_parInd", training_var_amount, ".RData"))
