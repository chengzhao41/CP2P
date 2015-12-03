source("Common/load_library.R")

### User supplies these values ###
args <- vector()

args[1] = 1
args[2] = 1
PARTITION_BEGIN = as.integer(args[1])
PARTITION_END = as.integer(args[2]) 
args[3] = "slope"
args[4] = 3
### End ###

load("Epirubicin/WS/epirubicin_data.RData")

if (args[3] == "slope") {
  input_data <- epirubicin$combined_slope.sva
  input_label <- epirubicin.labels$slope_combined
  input_partition <- partition$slope
  feature.l1000 <- feature.l1000$cp
} else if (args[3] == "auc") {
  input_data <- epirubicin$combined_AUC.sva
  input_label <- epirubicin.labels$AUC_combined
  input_partition <- partition$AUC
  feature.l1000 <- feature.l1000$cp
} else {
  stopifnot(FALSE)
}

if (length(args) == 4) {
  training_var_amount <- as.integer(args[4])
  print(paste0("epirubicin_cpp_", input_partition_BEGIN, "to", input_partition_END, "_", args[3], "_var", training_var_amount, ".RData"))  
  if (args[3] == "slope") {    
    input_partition$cpp <- partition_var$cpp[[training_var_amount]]$slope
  } else if (args[3] == "auc") {    
    input_partition$cpp <- partition_var$cpp[[training_var_amount]]$AUC
  } else {
    stop("arg 3 is wrong")
  }  
} else {
  stop("args != 4")
}

rm(epirubicin)

# do the computation 
snf.parameter <- seq(from = 5, to = 30, by = 5)
source("common/cpp_compute.R")
##

OUTPUT_DIR = "Epirubicin/output_WS/"
if (exists("training_var_amount")) {
  save.image(paste0(OUTPUT_DIR, "epirubicin_cpp_", input_partition_BEGIN, "to", input_partition_END, "_", args[3], "_var", training_var_amount, ".RData"))
} else {  
  save.image(paste0(OUTPUT_DIR, "epirubicin_cpp_", input_partition_BEGIN, "to", input_partition_END, "_", args[3], ".RData"))
}
