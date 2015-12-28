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
# args[4] = "Bortezomib"
# args[5] = "P2P"

### End ###

print(args)

training_var_amount <- as.integer(args[[3]])
  
if (args[4] == "Bortezomib") {
  load("Bortezomib/WS/bortezomib_data.RData")
  
  if (args[5] == "p2p") {
    stopifnot(training_var_amount > length(partition$cell_lines_all))
    
    input_data <- bortezomib$patient.combat
    input_label <- bortezomib.labels$patient
    input_partition = partition$cell_lines_all[[training_var_amount]]$p2p
    feature.l1000 <- feature.l1000$pp
  } else if (args[5] == "cp2p_slope") {
    stopifnot(training_var_amount > length(partition$cell_lines_all))
    
    input_data <- bortezomib$combined_slope.sva 
    input_label <- bortezomib.labels$slope_combined
    input_partition = partition$cell_lines_all[[training_var_amount]]$cp2p.slope
    feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc") {
    stopifnot(training_var_amount > length(partition$cell_lines_all))
    
    input_data <- bortezomib$combined_AUC.sva
    input_label <- bortezomib.labels$AUC_combined
    input_partition = partition$cell_lines_all[[training_var_amount]]$cp2p.AUC
    feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50") {
    stopifnot(training_var_amount > length(partition$cell_lines_all))
    
    input_data <- bortezomib$combined_IC50.sva
    input_label <- bortezomib.labels$IC50_combined
    input_partition = partition$cell_lines_all[[training_var_amount]]$cp2p.IC50
    feature.l1000 <- feature.l1000$cp
  } else {
    stop(paste("args[5]", args[5], "is invalid."))
  }
  
  input.type_measure = "auc"
  snf.parameter <- seq(from = 5, to = 30, by = 5)
  rm(bortezomib)
} else {
  stop(paste("args[4]", args[4], "is invalid."))
}

# do the computation
print(paste0(args[4], "_", args[5], PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))
source("Common/train_and_predict.R")
##

save.image(paste0(args[4], "_", args[5], PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))
