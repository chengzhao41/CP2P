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

print(paste("begin", PARTITION_BEGIN, "end", PARTITION_END))
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
} else if (args[3] == "ic50") {
  input_data <- epirubicin$combined_IC50.sva
  input_label <- epirubicin.labels$IC50_combined
  input_partition <- partition$IC50
  feature.l1000 <- feature.l1000$cp
} else {
  stop("invalid args[3]")
  input_data <- epirubicin$combined_AUC.sva
  input_label <- epirubicin.labels$AUC_combined
  input_partition <- partition$AUC
  feature.l1000 <- feature.l1000$cp
  input_partition <- list()
  input_partition$cp <- list(list(test_index = which(epirubicin.labels$slope_combined.source == "patient"), 
                                  training_index.single = which(epirubicin.labels$slope_combined.source != "patient")))
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
  print(paste0("epirubicin_cp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))
} else {
  stop("args not 4")
}
rm(epirubicin)

# do the computation
source("Common/cp_compute.R")
# end

OUTPUT_DIR = "Epirubicin/output_WS/"
if (exists("training_var_amount")) {
  save.image(paste0(OUTPUT_DIR, "epirubicin_cp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))  
} else {
  save.image(paste0(OUTPUT_DIR, "epirubicin_cp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], ".RData"))  
}
