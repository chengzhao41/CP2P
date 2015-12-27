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
stopifnot(length(args) == 3)
load("Bortezomib/WS/bortezomib_data.RData")

input_data <- bortezomib$patient.combat
input_label <- bortezomib.labels$patient  

training_var_amount <- as.integer(args[[3]])  
if (!is.na(training_var_amount)) {
  input_partition = partition$cell_lines_all[[training_var_amount]]$p2p
}  else {
  stop("arg 3 is wrong")
}
print(paste0("bortezomib_pp_", PARTITION_BEGIN, "to", PARTITION_END, "_var", training_var_amount, ".RData"))

rm(bortezomib)

print(paste("begin", PARTITION_BEGIN, "end", PARTITION_END))

feature.l1000 <- feature.l1000$pp
stopifnot(!is.null(feature.l1000))

# do the computation 
input.type_measure = "auc"
snf.parameter <- seq(from = 5, to = 30, by = 5)
source("common/train_and_predict.R")
##

if (exists("training_var_amount")) {
  save.image(paste0("bortezomib_p2p_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))  
} else {
  save.image(paste0("bortezomib_p2p_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], ".RData"))  
}
