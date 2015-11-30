source("Common/load_library.R")

### User supplies these values ###
args <- vector()

args[1] = 1
args[2] = 1
PARTITION_BEGIN = as.integer(args[1])
PARTITION_END = as.integer(args[2]) 
args[3] = "slope"
args[4] = 2 
### End ###

load("Bortezomib/WS/bortezomib_data.RData")

print(paste("begin", PARTITION_BEGIN, "end", PARTITION_END))

if (args[3] == "slope" || args[3] == "slope_p100" || args[3] == "slope_p50") {
  input_data <- bortezomib$combined_slope.sva
  input_label <- bortezomib.labels$slope_combined
  input_partition <- partition$slope
} else if (args[3] == "auc" || args[3] == "auc_p100" || args[3] == "auc_p50") {
  input_data <- bortezomib$combined_AUC.sva
  input_label <- bortezomib.labels$AUC_combined
  input_partition <- partition$AUC
} else if (args[3] == "ic50" || args[3] == "ic50_p100" || args[3] == "ic50_p50") {
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
  } else if (args[3] == "slope_p100") {    
    input_partition$cpp <- partition_var$cVar_p100_p[[training_var_amount]]$slope
  } else if (args[3] == "auc_p100") {    
    input_partition$cpp <- partition_var$cVar_p100_p[[training_var_amount]]$AUC
  } else if (args[3] == "ic50_p100") {    
    input_partition$cpp <- partition_var$cVar_p100_p[[training_var_amount]]$IC50
  } else if (args[3] == "slope_p50") {    
    input_partition$cpp <- partition_var$cVar_p50_p[[training_var_amount]]$slope
  } else if (args[3] == "auc_p50") {    
    input_partition$cpp <- partition_var$cVar_p50_p[[training_var_amount]]$AUC
  } else if (args[3] == "ic50_p50") {    
    input_partition$cpp <- partition_var$cVar_p50_p[[training_var_amount]]$IC50
  } else {
    stop("arg 3 is wrong")
  }  
}
remove (bortezomib)
feature.l1000 <- feature.l1000$cp

# do the computation 
source("common/cpp_compute.R")
##

OUTPUT_DIR = "Bortezomib/output_WS/"

if (exists("training_var_amount")) {
  save.image(paste0(OUTPUT_DIR, "bortezomib_cpp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))
} else {  
  save.image(paste0(OUTPUT_DIR, "bortezomib_cpp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], ".RData"))
}
