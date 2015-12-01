source("Common/load_library.R")

### User supplies these values ###
args <- vector()

args[1] = 1
args[2] = 2
PARTITION_BEGIN = as.integer(args[1])
PARTITION_END = as.integer(args[2]) 
args[3] = "slope"
args[4] = 10 
### End ###

load("Docetaxel/WS/docetaxel_data.RData")

INPUT_NFOLDS = 5
if (args[3] == "slope" || args[3] == "slope_p20") {
  input_data <- docetaxel$combined_slope.ComBat
  input_label <- docetaxel.labels$slope_combined
  input_partition <- partition$slope
} else if (args[3] == "auc" || args[3] == "auc_p20") {
  input_data <- docetaxel$combined_AUC.ComBat
  input_label <- docetaxel.labels$AUC_combined
  input_partition <- partition$AUC
} else if (args[3] == "ic50" || args[3] == "ic50_p20") {
  input_data <- docetaxel$combined_IC50.ComBat
  input_label <- docetaxel.labels$IC50_combined
  input_partition <- partition$IC50
} else if (args[3] == "slope_breast") {
  input_data <- docetaxel$breast_slope.ComBat
  input_label <- docetaxel.labels$slope_breast
  input_partition <- partition$slope_breast
  INPUT_NFOLDS = 3
} else if (args[3] == "auc_breast") {
  input_data <- docetaxel$breast_AUC.ComBat
  input_label <- docetaxel.labels$AUC_breast
  input_partition <- partition$AUC_breast
  INPUT_NFOLDS = 3
} else if (args[3] == "ic50_breast") {
  input_data <- docetaxel$breast_IC50.ComBat
  input_label <- docetaxel.labels$IC50_breast
  input_partition <- partition$IC50_breast
  INPUT_NFOLDS = 3
} else {
  stopifnot(FALSE)
}

if (length(args) == 4) {
  training_var_amount <- as.integer(args[[4]])
  print(paste0("docetaxel_cpp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))  
  if (args[3] == "slope") {    
    input_partition$cpp <- partition_var$cpp[[training_var_amount]]$slope
  } else if (args[3] == "auc") {    
    input_partition$cpp <- partition_var$cpp[[training_var_amount]]$AUC
  } else if (args[3] == "ic50") {    
    input_partition$cpp <- partition_var$cpp[[training_var_amount]]$IC50
  } else if (args[3] == "slope_breast") {    
    input_partition$cpp <- partition_var$cpp_breast_only[[training_var_amount]]$slope
    temp.max <- length(partition_var$cpp_breast_only[[1]]$IC50[[1]]$training_index.single)
    snf.parameter <- seq(from = 5, to = temp.max, by = 5)
  } else if (args[3] == "auc_breast") {    
    input_partition$cpp <- partition_var$cpp_breast_only[[training_var_amount]]$AUC
    temp.max <- length(partition_var$cpp_breast_only[[1]]$AUC[[1]]$training_index.single)
    snf.parameter <- seq(from = 5, to = temp.max, by = 5)
  } else if (args[3] == "ic50_breast") {    
    input_partition$cpp <- partition_var$cpp_breast_only[[training_var_amount]]$IC50
    temp.max <- length(partition_var$cpp_breast_only[[1]]$slope[[1]]$training_index.single)
    snf.parameter <- seq(from = 5, to = temp.max, by = 5)
  } else if (args[3] == "slope_p20") {    
    input_partition$cpp <- partition_var$cVar_p20_p[[training_var_amount]]$slope
  } else if (args[3] == "auc_p20") {    
    input_partition$cpp <- partition_var$cVar_p20_p[[training_var_amount]]$AUC
  } else if (args[3] == "ic50_p20") {    
    input_partition$cpp <- partition_var$cVar_p20_p[[training_var_amount]]$IC50
  } else {
    stop("arg 3 is wrong")
  }  
}

feature.l1000 <- feature.l1000$cp
remove(docetaxel)

# do the computation 
snf.parameter <- seq(from = 5, to = 30, by = 5)
source("common/cpp_compute.R")
##

OUTPUT_DIR = "Docetaxel/output_WS/"
if (exists("training_var_amount")) {
  save.image(paste0(OUTPUT_DIR, "docetaxel_cpp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))  
} else {
  save.image(paste0(OUTPUT_DIR, "docetaxel_cpp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], ".RData"))  
}
