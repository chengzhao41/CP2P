source("Common/load_library.R")

### User supplies these values ###
args <- vector()

args[1] = 1
args[2] = 2
PARTITION_BEGIN = as.integer(args[1])
PARTITION_END = as.integer(args[2]) 
args[3] = "slope"
args[4] = 3
### End ###

load("Erlotinib/WS/erlotinib_data.RData")

if (args[3] == "slope" || args[3] == "slope_p20" || args[3] == "slope_p20_lusc") {
  input_data <- erlotinib$all.ComBat.slope
  input_label <- erlotinib.labels$all.slope
  input_partition <- partition$all.slope
  feature.l1000 <- feature.l1000$cp
} else if (args[3] == "ic50" || args[3] == "auc" 
           || args[3] == "ic50_p20" || args[3] == "auc_p20"
           || args[3] == "ic50_p20_lusc" || args[3] == "auc_p20_lusc") {
  input_data <- erlotinib$all.ComBat.IC50
  input_label <- erlotinib.labels$all.IC50
  input_partition <- partition$all.IC50
  feature.l1000 <- feature.l1000$cp
} else if (args[3] == "lung_slope") {
  input_data <- erlotinib$lung_all.slope.ComBat
  input_label <- erlotinib.labels$lung_all.slope
  input_partition <- partition$lung.slope
  feature.l1000 <- feature.l1000$cp
} else if ((args[3] == "lung_ic50") || (args[3] == "lung_auc")) {
  input_data <- erlotinib$lung_all.IC50.ComBat
  input_label <- erlotinib.labels$lung_all.IC50
  input_partition <- partition$lung.IC50
  feature.l1000 <- feature.l1000$cp
} else if (args[3] == "GEO") {
  input_data <- erlotinib$combined_GEO
  input_label <- erlotinib.labels$GEO_IC50_combined
  input_partition <- partition$GEO
  feature.l1000 <- feature.l1000$geo_cp  
} else {
  stop("false args[3]")
}

if (length(args) == 4) {
  training_var_amount <- as.integer(args[[4]])
  print(paste0("erlotinib_cpp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))  
  if (args[3] == "slope") {    
    input_partition$cpp <- partition_var$cpp[[training_var_amount]]$slope
  } else if ((args[3] == "ic50") || (args[3] == "auc")) {    
    input_partition$cpp <- partition_var$cpp[[training_var_amount]]$IC50
  } else if (args[3] == "lung_slope") {    
    input_partition$cpp <- partition_var$cpp_lung_only[[training_var_amount]]$slope
  } else if ((args[3] == "lung_ic50") || (args[3] == "lung_auc")) {    
    input_partition$cpp <- partition_var$cpp_lung_only[[training_var_amount]]$IC50
  } else if (args[3] == "slope_p20") {    
    input_partition$cpp <- partition_var$cVar_p20_p[[training_var_amount]]$slope
  } else if (args[3] == "ic50_p20" || (args[3] == "auc_p20")) {    
    input_partition$cpp <- partition_var$cVar_p20_p[[training_var_amount]]$IC50
  } else if (args[3] == "slope_p20_lusc") {    
    input_partition$cpp <- partition_var$cVar_p20_p_lusc[[training_var_amount]]$slope
  } else if (args[3] == "ic50_p20_lusc" || (args[3] == "auc_p20_lusc")) {    
    input_partition$cpp <- partition_var$cVar_p20_p_lusc[[training_var_amount]]$IC50
  } else {
    stop("arg 3 is wrong")
  }
  stopifnot(!is.null(input_partition$cpp))
} else {
  stop("not 4 args")
}
stopifnot(exists("training_var_amount"))
rm(erlotinib)

# do the computation 
snf.parameter <- seq(from = 5, to = 30, by = 5)
source("common/cpp_compute.R")
##

OUTPUT_DIR = "Erlotinib/output_WS/"
if (exists("training_var_amount")) {
  save.image(paste0(OUTPUT_DIR, "erlotinib_cpp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], "_var", training_var_amount, ".RData"))  
} else {
  save.image(paste0(OUTPUT_DIR, "erlotinib_cpp_", PARTITION_BEGIN, "to", PARTITION_END, "_", args[3], ".RData"))  
}
