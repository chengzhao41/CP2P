source("Common/load_library.R")

### User supplies these values ###
args <- vector()

args[1] = 1 # PARTITION_BEGIN
args[2] = 2 # PARTITION_END
PARTITION_BEGIN = as.integer(args[1])
PARTITION_END = as.integer(args[2]) 
args[3] = "slope" # cell line response type
args[4] = 3 # partition_ind to use
### End ###

load("Erlotinib/WS/erlotinib_data.RData")

if (args[3] == "slope") {
  input_data <- erlotinib$all.ComBat.slope
  input_label <- erlotinib.labels$all.slope
  input_partition <- partition$all.slope
  feature.l1000 <- feature.l1000$cp
} else if ((args[3] == "ic50") || (args[3] == "auc")) {
  input_data <- erlotinib$all.ComBat.IC50
  input_label <- erlotinib.labels$all.IC50
  input_partition <- partition$all.IC50
  feature.l1000 <- feature.l1000$cp
} else if (args[3] == "lung_slope") {
  input_data <- erlotinib$lung_all.slope.ComBat
  input_label <- erlotinib.labels$lung_all.slope
  #input_partition <- partition$lung.slope
  input_partition <- list()
  input_partition$cp <- list(list(test_index = which(erlotinib.labels$lung_all.slope.source == "patient"), 
                                  training_index.single = which(erlotinib.labels$lung_all.slope.source != "patient")))
  feature.l1000 <- feature.l1000$cp
} else if (args[3] == "lung_ic50") {
  input_data <- erlotinib$lung_all.IC50.ComBat
  input_label <- erlotinib.labels$lung_all.IC50
  #input_partition <- partition$lung.IC50
  input_partition <- list()
  input_partition$cp <- list(list(test_index = which(erlotinib.labels$lung_all.IC50.source == "patient"), 
                                  training_index.single = which(erlotinib.labels$lung_all.IC50.source != "patient")))
  feature.l1000 <- feature.l1000$cp
} else if (args[3] == "GEO") {
  input_data <- erlotinib$combined_GEO
  input_label <- erlotinib.labels$GEO_IC50_combined
  input_partition <- partition$GEO
  feature.l1000 <- feature.l1000$geo_cp  
} else {
  stop("false args[3]")
  input_partition <- list()
  input_partition$cp <- list(list(test_index = which(erlotinib.labels$IC50_breast.source == "patient"), 
                                  training_index.single = which(erlotinib.labels$IC50_breast.source != "patient")))
}

if (length(args) == 4) {
  training_var_amount <- as.integer(args[[4]])
  if (args[3] == "slope") {    
    input_partition$cp <- partition_var$cp[[training_var_amount]]$slope
  } else if ((args[3] == "ic50") || (args[3] == "auc")) {    
    input_partition$cp <- partition_var$cp[[training_var_amount]]$IC50
  } else if (args[3] == "lung_slope") {    
    input_partition$cp <- partition_var$cp_lung_only[[training_var_amount]]$slope
  } else if ((args[3] == "lung_ic50") || (args[3] == "auc")) {    
    input_partition$cp <- partition_var$cp_lung_only[[training_var_amount]]$IC50
  } else {
    stop("arg 3 is wrong")
  }
  print(paste0("erlotinib_cp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))
} 

rm(erlotinib)

source("Common/cp_compute.R")

print("completed!")
print(paste("BEGIN and END:", PARTITION_BEGIN, PARTITION_END))

rm(input_data)
OUTPUT_DIR = "Erlotinib/output_WS/"
if (exists("training_var_amount")) {
  save.image(paste0(OUTPUT_DIR, "erlotinib_cp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))  
} else {
  save.image(paste0(OUTPUT_DIR, "erlotinib_cp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], ".RData"))  
}
