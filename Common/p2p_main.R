rm(list = ls())
setwd("/home/zhaoche7/drp3")
source("Common/load_library.R")

### User supplies these values ###
args <- commandArgs(trailingOnly = TRUE)
print(length(args))
PARTITION_BEGIN = as.integer(args[1])
PARTITION_END = as.integer(args[2]) 
registerDoParallel(8)

#args <- vector()
#registerDoParallel(4)
#args[1] = 1 # partition start
#args[2] = 1 # partition end
#args[3] = 10 # index for training sets

### End ###
stopifnot(length(args) == 4)
print(args)

if (args[4] == "bortezomib") {
  load("Bortezomib/WS/bortezomib_data.RData")
  input_data <- bortezomib$patient.combat
  input_label <- bortezomib.labels$patient
  input.type_measure = "auc"
  snf.parameter <- seq(from = 5, to = 30, by = 5)
  rm(bortezomib)
} else {
  stop(paste("args[4]", args[4], "is invalid."))
}

training_var_amount <- as.integer(args[[3]])
input_partition = partition$cell_lines_all[[training_var_amount]]$p2p

# do the computation
print(paste0(args[4], "_p2p_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))
feature.l1000 <- feature.l1000$pp
source("common/train_and_predict.R")
##

save.image(paste0(args[4], "_p2p_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))
