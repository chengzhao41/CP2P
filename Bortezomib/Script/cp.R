source("Common/load_library.R")

### User supplies these values ###
args <- vector()

args[1] = 1 # PARTITION_BEGIN
args[2] = 1 # PARTITION_END
PARTITION_BEGIN = as.integer(args[1])
PARTITION_END = as.integer(args[2]) 
args[3] = "slope" # cell line response type
args[4] = 10 # partition_ind to use
### End ###

print(paste("begin", PARTITION_BEGIN, "end", PARTITION_END))
load("Bortezomib/WS/bortezomib_data.RData")

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
} else if (args[3] == "slope_once") {
  input_data <- bortezomib$combined_slope.sva
  input_label <- bortezomib.labels$slope_combined
  input_partition <- list()
  input_partition$cp <- list(list(test_index = which(bortezomib.labels$slope_combined.source == "patient"), 
                                  training_index.single = which(bortezomib.labels$slope_combined.source != "patient")))
} else if (args[3] == "auc_once") {
  input_data <- bortezomib$combined_AUC.sva
  input_label <- bortezomib.labels$AUC_combined
  input_partition <- list()
  input_partition$cp <- list(list(test_index = which(bortezomib.labels$AUC_combined.source == "patient"), 
                                  training_index.single = which(bortezomib.labels$AUC_combined.source != "patient")))
} else if (args[3] == "ic50_once") {
  input_data <- bortezomib$combined_IC50.sva
  input_label <- bortezomib.labels$IC50_combined
  input_partition <- list()
  input_partition$cp <- list(list(test_index = which(bortezomib.labels$IC50_combined.source == "patient"), 
                                  training_index.single = which(bortezomib.labels$IC50_combined.source != "patient")))
} else {
  stopifnot(FALSE)
  input_data <- bortezomib$combined_slope.sva
  input_label <- bortezomib.labels$slope_combined
  input_partition <- list(list(test_index = which(bortezomib.labels$slope_combined.source == "patient"), 
                               training_index.single = which(bortezomib.labels$slope_combined.source != "patient")))
}

if (length(args) == 4) {
  training_var_amount <- as.integer(args[[4]])
  if (args[3] == "slope") {    
    input_partition$cp <- partition_var$cp[[training_var_amount]]$slope
  } else if (args[3] == "auc") {    
    input_partition$cp <- partition_var$cp[[training_var_amount]]$AUC
  } else if (args[3] == "ic50") {    
    input_partition$cp <- partition_var$cp[[training_var_amount]]$IC50
  } else {
    stop("arg 3 is wrong")
  }
  print(paste0("bortezomib_cp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))
}

feature.l1000 <- feature.l1000$cp
remove(bortezomib)

# do the computations
snf.parameter <- seq(from = 5, to = 30, by = 5)
source("Common/cp_compute.R")

print("completed!")
print(paste("BEGIN and END:", PARTITION_BEGIN, PARTITION_END))

### output
rm(input_data)
OUTPUT_DIR = "Bortezomib/output_WS/"
if (exists("training_var_amount")) {
  save.image(file = paste0(OUTPUT_DIR, "bortezomib_cp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))  
} else {
  save.image(file = paste0(OUTPUT_DIR, "bortezomib_cp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], ".RData"))  
}
