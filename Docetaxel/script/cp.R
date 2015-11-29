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
load("Docetaxel/WS/docetaxel_data.RData")

if (args[3] == "slope") {
  input_data <- docetaxel$combined_slope.ComBat
  input_label <- docetaxel.labels$slope_combined
  input_partition <- partition$slope
} else if (args[3] == "auc") {
  input_data <- docetaxel$combined_AUC.ComBat
  input_label <- docetaxel.labels$AUC_combined
  input_partition <- partition$AUC
} else if (args[3] == "ic50") {
  input_data <- docetaxel$combined_IC50.ComBat
  input_label <- docetaxel.labels$IC50_combined
  input_partition <- partition$IC50
} else if (args[3] == "slope_breast") {
  input_data <- docetaxel$breast_slope.ComBat
  input_label <- docetaxel.labels$slope_breast
  input_partition <- list()
  input_partition$cp <- list(list(test_index = which(docetaxel.labels$slope_breast.source == "patient"), 
                                  training_index.single = which(docetaxel.labels$slope_breast.source != "patient")))
} else if (args[3] == "auc_breast") {
  input_data <- docetaxel$breast_AUC.ComBat
  input_label <- docetaxel.labels$AUC_breast
  input_partition <- list()
  input_partition$cp <- list(list(test_index = which(docetaxel.labels$AUC_breast.source == "patient"), 
                                  training_index.single = which(docetaxel.labels$AUC_breast.source != "patient")))
} else if (args[3] == "ic50_breast") {
  input_data <- docetaxel$breast_IC50.ComBat
  input_label <- docetaxel.labels$IC50_breast
  input_partition <- list()
  input_partition$cp <- list(list(test_index = which(docetaxel.labels$IC50_breast.source == "patient"), 
                                  training_index.single = which(docetaxel.labels$IC50_breast.source != "patient")))
} else {
  stopifnot(FALSE)
  input_data <- docetaxel$breast_IC50.ComBat
  input_label <- docetaxel.labels$IC50_breast
  input_partition <- list()
  input_partition$cp <- list(list(test_index = which(docetaxel.labels$IC50_breast.source == "patient"), 
                                  training_index.single = which(docetaxel.labels$IC50_breast.source != "patient")))
  #AUC
  input_data <- docetaxel$breast_AUC.ComBat
  input_label <- docetaxel.labels$AUC_breast
  input_partition <- list()
  input_partition$cp <- list(list(test_index = which(docetaxel.labels$AUC_breast.source == "patient"), 
                                  training_index.single = which(docetaxel.labels$AUC_breast.source != "patient")))
}
remove(docetaxel)

if (length(args) == 4) {
  training_var_amount <- as.integer(args[[4]])
  if (args[3] == "slope") {    
    input_partition$cp <- partition_var$cp[[training_var_amount]]$slope
  } else if (args[3] == "auc") {    
    input_partition$cp <- partition_var$cp[[training_var_amount]]$AUC
  } else if (args[3] == "ic50") {    
    input_partition$cp <- partition_var$cp[[training_var_amount]]$IC50
  } else if (args[3] == "slope_breast") {    
    input_partition$cp <- partition_var$cp_breast_only[[training_var_amount]]$slope
    temp.max <- length(partition_var$cp_breast_only[[1]]$IC50[[1]]$training_index.single)
    snf.parameter <- seq(from = 5, to = temp.max, by = 5)
  } else if (args[3] == "auc_breast") {    
    input_partition$cp <- partition_var$cp_breast_only[[training_var_amount]]$AUC
    temp.max <- length(partition_var$cp_breast_only[[1]]$AUC[[1]]$training_index.single)
    snf.parameter <- seq(from = 5, to = temp.max, by = 5)
  } else if (args[3] == "ic50_breast") {    
    input_partition$cp <- partition_var$cp_breast_only[[training_var_amount]]$IC50
    temp.max <- length(partition_var$cp_breast_only[[1]]$slope[[1]]$training_index.single)
    snf.parameter <- seq(from = 5, to = temp.max, by = 5)
  }  else {
    stop("arg 3 is wrong")
  }
  print(paste0("docetaxel_cp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))
}

feature.l1000 <- feature.l1000$cp

# do the computations
source("Common/cp_compute.R")
#

OUTPUT_DIR = "Docetaxel/output_WS/"
if (exists("training_var_amount")) {
  save.image(paste0(OUTPUT_DIR, "docetaxel_cp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))  
} else {  
  save.image(paste0(OUTPUT_DIR, "docetaxel_cp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], ".RData"))  
}
