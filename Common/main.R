rm(list = ls())
setwd("/home/zhaoche7/drp3")
source("Common/load_library.R")
registerDoParallel(8)

### User supplies these values ###
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5)

# args <- vector()
# registerDoParallel(4)
# args[1] = "1" # partition start
# args[2] = "2" # partition end
# args[3] = "4" # index for training sets
# args[4] = "docetaxel"
# args[5] = "cp2p_slope_breast"

PARTITION_BEGIN = as.integer(args[1])
PARTITION_END = as.integer(args[2])

### End ###

print(args)

training_var_amount <- as.integer(args[[3]])
INPUT.NFOLDS = 5  
if (args[4] == "bortezomib") {
  load("Bortezomib/WS/bortezomib_data.RData")
  
  input.type_measure = "auc"
  input_snf.parameter <- seq(from = 5, to = 30, by = 5)
  
  if (args[5] == "p2p") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- bortezomib$patient.combat
    input_label <- bortezomib.labels$patient
    input_partition <- partition$cell_lines_all[[training_var_amount]]$p2p
    input_feature.l1000 <- feature.l1000$pp
  } else if (args[5] == "cp2p_slope") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- bortezomib$slope_combined.sva 
    input_label <- bortezomib.labels$slope_combined
    input_partition = partition$cell_lines_all[[training_var_amount]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- bortezomib$AUC_combined.sva
    input_label <- bortezomib.labels$AUC_combined
    input_partition = partition$cell_lines_all[[training_var_amount]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- bortezomib$IC50_combined.sva
    input_label <- bortezomib.labels$IC50_combined
    input_partition = partition$cell_lines_all[[training_var_amount]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50_p50") {
    stopifnot(training_var_amount <= length(partition$patient_50))
    
    input_data <- bortezomib$IC50_combined.sva
    input_label <- bortezomib.labels$IC50_combined
    input_partition = partition$patient_50[[training_var_amount]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_p50") {
    stopifnot(training_var_amount <= length(partition$patient_50))
    
    input_data <- bortezomib$AUC_combined.sva
    input_label <- bortezomib.labels$AUC_combined
    input_partition = partition$patient_50[[training_var_amount]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_slope_p50") {
    stopifnot(training_var_amount <= length(partition$patient_50))
    
    input_data <- bortezomib$slope_combined.sva
    input_label <- bortezomib.labels$slope_combined
    input_partition = partition$patient_50[[training_var_amount]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_ic50_p50") {
    stopifnot(training_var_amount <= length(partition$patient_50))
    
    input_data <- bortezomib$IC50_combined.sva
    input_label <- bortezomib.labels$IC50_combined
    input_partition = partition$patient_50[[training_var_amount]]$c2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_auc_p50") {
    stopifnot(training_var_amount <= length(partition$patient_50))
    
    input_data <- bortezomib$AUC_combined.sva
    input_label <- bortezomib.labels$AUC_combined
    input_partition = partition$patient_50[[training_var_amount]]$c2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_slope_p50") {
    stopifnot(training_var_amount <= length(partition$patient_50))
    
    input_data <- bortezomib$slope_combined.sva
    input_label <- bortezomib.labels$slope_combined
    input_partition = partition$patient_50[[training_var_amount]]$c2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_slope_p100") {
    stopifnot(training_var_amount <= length(partition$patient_100))
    
    input_data <- bortezomib$slope_combined.sva
    input_label <- bortezomib.labels$slope_combined
    input_partition = partition$patient_100[[training_var_amount]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_p100") {
    stopifnot(training_var_amount <= length(partition$patient_100))
    
    input_data <- bortezomib$AUC_combined.sva
    input_label <- bortezomib.labels$AUC_combined
    input_partition = partition$patient_100[[training_var_amount]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50_p100") {
    stopifnot(training_var_amount <= length(partition$patient_100))
    
    input_data <- bortezomib$IC50_combined.sva
    input_label <- bortezomib.labels$IC50_combined
    input_partition = partition$patient_100[[training_var_amount]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else {
    stop(paste("args[5]", args[5], "is invalid."))
  }
  rm(bortezomib, bortezomib.labels)
} else if (args[4] == "docetaxel") {
  load("Docetaxel/WS/docetaxel_data.RData")
  
  input.type_measure = "acc"
  input_snf.parameter <- seq(from = 5, to = 30, by = 5)
  
  if (args[5] == "p2p") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- docetaxel$patient
    input_label <- docetaxel.labels$patient
    input_partition <- partition$cell_lines_all[[training_var_amount]]$p2p
    input_feature.l1000 <- feature.l1000$pp
    input_snf.parameter <- seq(from = 5, to = length(input_partition[[1]]$training_index), by = 3)
  } else if (args[5] == "cp2p_slope") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- docetaxel$slope_combined.ComBat 
    input_label <- docetaxel.labels$slope_combined
    input_partition = partition$cell_lines_all[[training_var_amount]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- docetaxel$AUC_combined.ComBat
    input_label <- docetaxel.labels$AUC_combined
    input_partition = partition$cell_lines_all[[training_var_amount]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- docetaxel$IC50_combined.ComBat
    input_label <- docetaxel.labels$IC50_combined
    input_partition = partition$cell_lines_all[[training_var_amount]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_slope_breast") {
    stopifnot(training_var_amount <= length(partition$cell_lines_breast))
    
    input_data <- docetaxel$slope_breast.ComBat 
    input_label <- docetaxel.labels$slope_breast.combined
    input_partition = partition$cell_lines_breast[[training_var_amount]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_breast") {
    stopifnot(training_var_amount <= length(partition$cell_lines_breast))
    
    input_data <- docetaxel$AUC_breast.ComBat
    input_label <- docetaxel.labels$AUC_breast.combined
    input_partition = partition$cell_lines_breast[[training_var_amount]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50_breast") {
    stopifnot(training_var_amount <= length(partition$cell_lines_breast))
    
    input_data <- docetaxel$IC50_breast.ComBat
    input_label <- docetaxel.labels$IC50_breast.combined
    input_partition = partition$cell_lines_breast[[training_var_amount]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50_p23") {
    stopifnot(training_var_amount <= length(partition$patient_23))
    
    input_data <- docetaxel$IC50_combined.ComBat
    input_label <- docetaxel.labels$IC50_combined
    input_partition = partition$patient_23[[training_var_amount]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_p23") {
    stopifnot(training_var_amount <= length(partition$patient_23))
    
    input_data <- docetaxel$AUC_combined.ComBat
    input_label <- docetaxel.labels$AUC_combined
    input_partition = partition$patient_23[[training_var_amount]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_slope_p23") {
    stopifnot(training_var_amount <= length(partition$patient_23))
    
    input_data <- docetaxel$slope_combined.ComBat
    input_label <- docetaxel.labels$slope_combined
    input_partition = partition$patient_23[[training_var_amount]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_ic50_p23") {
    stopifnot(training_var_amount <= length(partition$patient_23))
    
    input_data <- docetaxel$IC50_combined.ComBat
    input_label <- docetaxel.labels$IC50_combined
    input_partition = partition$patient_23[[training_var_amount]]$c2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_auc_p23") {
    stopifnot(training_var_amount <= length(partition$patient_23))
    
    input_data <- docetaxel$AUC_combined.ComBat
    input_label <- docetaxel.labels$AUC_combined
    input_partition = partition$patient_23[[training_var_amount]]$c2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_slope_p23") {
    stopifnot(training_var_amount <= length(partition$patient_23))
    
    input_data <- docetaxel$slope_combined.ComBat
    input_label <- docetaxel.labels$slope_combined
    input_partition = partition$patient_23[[training_var_amount]]$c2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else {
    stop(paste("args[5]", args[5], "is invalid."))
  }
  rm(docetaxel, docetaxel.labels)
} else if (args[4] == "erlotinib_gdsc" || args[4] == "erlotinib_ccle") {
  
  if (args[4] == "erlotinib_gdsc") {
    load("Erlotinib/WS/erlotinib_homogenized_data_gdsc.RData")
  } else {
    load("Erlotinib/WS/erlotinib_homogenized_data_ccle.RData")
  }
  
  input.type_measure = "acc"
  input_snf.parameter <- seq(from = 5, to = 30, by = 5)
  
  if (args[5] == "p2p") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- erlotinib$patient
    input_label <- erlotinib.labels$patient
    input_partition <- partition$cell_lines_all[[training_var_amount]]$p2p
    input_feature.l1000 <- feature.l1000$pp
    input_snf.parameter <- seq(from = 5, to = length(input_partition[[1]]$training_index), by = 3)
  } else if (args[5] == "cp2p_slope") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- erlotinib$slope_combined.ComBat 
    input_label <- erlotinib.labels$slope_combined
    input_partition = partition$cell_lines_all[[training_var_amount]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- erlotinib$AUC_combined.ComBat
    input_label <- erlotinib.labels$AUC_combined
    input_partition = partition$cell_lines_all[[training_var_amount]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- erlotinib$IC50_combined.ComBat
    input_label <- erlotinib.labels$IC50_combined
    input_partition = partition$cell_lines_all[[training_var_amount]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_slope_lung") {
    stopifnot(training_var_amount <= length(partition$cell_lines_lung))
    
    input_data <- erlotinib$slope_lung.ComBat 
    input_label <- erlotinib.labels$slope_lung
    input_partition = partition$cell_lines_lung[[training_var_amount]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_lung") {
    stopifnot(training_var_amount <= length(partition$cell_lines_lung))
    
    input_data <- erlotinib$AUC_lung.ComBat
    input_label <- erlotinib.labels$AUC_lung
    input_partition = partition$cell_lines_lung[[training_var_amount]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50_lung") {
    stopifnot(training_var_amount <= length(partition$cell_lines_lung))
    
    input_data <- erlotinib$IC50_lung.ComBat
    input_label <- erlotinib.labels$IC50_lung
    input_partition = partition$cell_lines_lung[[training_var_amount]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50_p24") {
    stopifnot(training_var_amount <= length(partition$patient_24))
    
    input_data <- erlotinib$IC50_combined.ComBat
    input_label <- erlotinib.labels$IC50_combined
    input_partition = partition$patient_24[[training_var_amount]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_p24") {
    stopifnot(training_var_amount <= length(partition$patient_24))
    
    input_data <- erlotinib$AUC_combined.ComBat
    input_label <- erlotinib.labels$AUC_combined
    input_partition = partition$patient_24[[training_var_amount]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_slope_p24") {
    stopifnot(training_var_amount <= length(partition$patient_24))
    
    input_data <- erlotinib$slope_combined.ComBat
    input_label <- erlotinib.labels$slope_combined
    input_partition = partition$patient_24[[training_var_amount]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_ic50_p24") {
    stopifnot(training_var_amount <= length(partition$patient_24))
    
    input_data <- erlotinib$IC50_combined.ComBat
    input_label <- erlotinib.labels$IC50_combined
    input_partition = partition$patient_24[[training_var_amount]]$c2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_auc_p24") {
    stopifnot(training_var_amount <= length(partition$patient_24))
    
    input_data <- erlotinib$AUC_combined.ComBat
    input_label <- erlotinib.labels$AUC_combined
    input_partition = partition$patient_24[[training_var_amount]]$c2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_slope_p24") {
    stopifnot(training_var_amount <= length(partition$patient_24))
    
    input_data <- erlotinib$slope_combined.ComBat
    input_label <- erlotinib.labels$slope_combined
    input_partition = partition$patient_24[[training_var_amount]]$c2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else {
    stop(paste("args[5]", args[5], "is invalid."))
  }
  rm(erlotinib, erlotinib.labels)
} else if (args[4] == "epirubicin") {
  load("Epirubicin/WS/epirubicin_data.RData")
  
  input.type_measure = "auc"
  input_snf.parameter <- seq(from = 5, to = 30, by = 5)
  
  if (args[5] == "p2p") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- epirubicin$patient
    input_label <- epirubicin.labels$patient
    input_partition <- partition$cell_lines_all[[training_var_amount]]$p2p
    input_feature.l1000 <- feature.l1000$pp
    INPUT.NFOLDS = 3 # too imbalanced to do 5
  } else if (args[5] == "cp2p_slope") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- epirubicin$slope_combined.sva 
    input_label <- epirubicin.labels$slope_combined
    input_partition = partition$cell_lines_all[[training_var_amount]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- epirubicin$AUC_combined.sva
    input_label <- epirubicin.labels$AUC_combined
    input_partition = partition$cell_lines_all[[training_var_amount]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
    INPUT.NFOLDS = 3 # too imbalanced to do 5
  } else if (args[5] == "cp2p_ic50") {
    stopifnot(training_var_amount <= length(partition$cell_lines_all))
    
    input_data <- epirubicin$IC50_combined.sva
    input_label <- epirubicin.labels$IC50_combined
    input_partition = partition$cell_lines_all[[training_var_amount]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50_p40") {
    stopifnot(training_var_amount <= length(partition$patient_40))
    
    input_data <- epirubicin$IC50_combined.sva
    input_label <- epirubicin.labels$IC50_combined
    input_partition = partition$patient_40[[training_var_amount]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_p40") {
    stopifnot(training_var_amount <= length(partition$patient_40))
    
    input_data <- epirubicin$AUC_combined.sva
    input_label <- epirubicin.labels$AUC_combined
    input_partition = partition$patient_40[[training_var_amount]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_slope_p40") {
    stopifnot(training_var_amount <= length(partition$patient_40))
    
    input_data <- epirubicin$slope_combined.sva
    input_label <- epirubicin.labels$slope_combined
    input_partition = partition$patient_40[[training_var_amount]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_ic50_p40") {
    stopifnot(training_var_amount <= length(partition$patient_40))
    
    input_data <- epirubicin$IC50_combined.sva
    input_label <- epirubicin.labels$IC50_combined
    input_partition = partition$patient_40[[training_var_amount]]$c2p.IC50
    input_feature.l1000 <- feature.l1000$cp
    INPUT.NFOLDS = 3 # too imbalanced to do 5
    stop("c2p_ic50_p40 labels are too imbalanced!")
  } else if (args[5] == "c2p_auc_p40") {
    stopifnot(training_var_amount <= length(partition$patient_40))
    
    input_data <- epirubicin$AUC_combined.sva
    input_label <- epirubicin.labels$AUC_combined
    input_partition = partition$patient_40[[training_var_amount]]$c2p.AUC
    input_feature.l1000 <- feature.l1000$cp
    INPUT.NFOLDS = 3 # too imbalanced to do 5
  } else if (args[5] == "c2p_slope_p40") {
    stopifnot(training_var_amount <= length(partition$patient_40))
    
    input_data <- epirubicin$slope_combined.sva
    input_label <- epirubicin.labels$slope_combined
    input_partition = partition$patient_40[[training_var_amount]]$c2p.slope
    input_feature.l1000 <- feature.l1000$cp
    INPUT.NFOLDS = 3 # too imbalanced to do 5
  } else if (args[5] == "cp2p_slope_p80") {
    stopifnot(training_var_amount <= length(partition$patient_80))
    
    input_data <- epirubicin$slope_combined.sva
    input_label <- epirubicin.labels$slope_combined
    input_partition = partition$patient_80[[training_var_amount]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_p80") {
    stopifnot(training_var_amount <= length(partition$patient_80))
    
    input_data <- epirubicin$AUC_combined.sva
    input_label <- epirubicin.labels$AUC_combined
    input_partition = partition$patient_80[[training_var_amount]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50_p80") {
    stopifnot(training_var_amount <= length(partition$patient_80))
    
    input_data <- epirubicin$IC50_combined.sva
    input_label <- epirubicin.labels$IC50_combined
    input_partition = partition$patient_80[[training_var_amount]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else {
    stop(paste("args[5]", args[5], "is invalid."))
  }
  rm(epirubicin, epirubicin.labels)
} else {
  stop(paste("args[4]", args[4], "is invalid."))
}

stopifnot(!is.null(input_data))
stopifnot(!is.null(input_label))
stopifnot(!is.null(input_partition))

rm(partition, feature.l1000)

# do the computation
print(paste0(args[4], "_", args[5], "_", PARTITION_BEGIN, "to", PARTITION_END, "_parInd", training_var_amount, ".RData"))
source("Common/train_and_predict.R")
##

save.image(paste0(args[4], "_", args[5], "_", PARTITION_BEGIN, "to", PARTITION_END, "_parInd", training_var_amount, ".RData"))
