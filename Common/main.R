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
# args[2] = "1" # partition end
# args[3] = "1" # index for training sets
# args[4] = "epirubicin_ComBat"
# args[5] = "cp2p_slope_once"

PARTITION_BEGIN = as.integer(args[1])
PARTITION_END = as.integer(args[2])

### End ###

print(args)

parInd <- as.integer(args[[3]])
INPUT.NFOLDS = 5  
train_once = FALSE
input.type_measure.test <- NULL
if (args[4] == "bortezomib") {
  load("Bortezomib/WS/bortezomib_data.RData")
  
  input.type_measure = "auc"
  input_snf.parameter <- seq(from = 5, to = 30, by = 5)
  
  if (args[5] == "p2p") {
    stopifnot(parInd <= length(partition$cell_lines_all))
    
    input_data <- bortezomib$patient.combat
    input_label <- bortezomib.labels$patient
    input_partition <- partition$cell_lines_all[[parInd]]$p2p
    input_feature.l1000 <- feature.l1000$pp
  } else if (args[5] == "cp2p_slope") {
    stopifnot(parInd <= length(partition$cell_lines_all))
    
    input_data <- bortezomib$slope_combined.sva 
    input_label <- bortezomib.labels$slope_combined
    input_partition = partition$cell_lines_all[[parInd]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc") {
    stopifnot(parInd <= length(partition$cell_lines_all))
    
    input_data <- bortezomib$AUC_combined.sva
    input_label <- bortezomib.labels$AUC_combined
    input_partition = partition$cell_lines_all[[parInd]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50") {
    stopifnot(parInd <= length(partition$cell_lines_all))
    
    input_data <- bortezomib$IC50_combined.sva
    input_label <- bortezomib.labels$IC50_combined
    input_partition = partition$cell_lines_all[[parInd]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50_p20") {
    stopifnot(parInd <= length(partition$patient_20))
    
    input_data <- bortezomib$IC50_combined.sva
    input_label <- bortezomib.labels$IC50_combined
    input_partition = partition$patient_20[[parInd]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_p20") {
    stopifnot(parInd <= length(partition$patient_20))
    
    input_data <- bortezomib$AUC_combined.sva
    input_label <- bortezomib.labels$AUC_combined
    input_partition = partition$patient_20[[parInd]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_slope_p20") {
    stopifnot(parInd <= length(partition$patient_20))
    
    input_data <- bortezomib$slope_combined.sva
    input_label <- bortezomib.labels$slope_combined
    input_partition = partition$patient_20[[parInd]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_ic50_p20") {
    stopifnot(parInd <= length(partition$patient_20))
    
    train_once = TRUE
    input_data <- bortezomib$IC50_combined.sva
    input_label <- bortezomib.labels$IC50_combined
    input_partition = partition$patient_20[[parInd]]$c2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_auc_p20") {
    stopifnot(parInd <= length(partition$patient_20))
    
    train_once = TRUE
    input_data <- bortezomib$AUC_combined.sva
    input_label <- bortezomib.labels$AUC_combined
    input_partition = partition$patient_20[[parInd]]$c2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_slope_p20") {
    stopifnot(parInd <= length(partition$patient_20))
    
    train_once = TRUE
    input_data <- bortezomib$slope_combined.sva
    input_label <- bortezomib.labels$slope_combined
    input_partition = partition$patient_20[[parInd]]$c2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_slope_p40") {
    stopifnot(parInd <= length(partition$patient_40))
    
    input_data <- bortezomib$slope_combined.sva
    input_label <- bortezomib.labels$slope_combined
    input_partition = partition$patient_40[[parInd]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_p40") {
    stopifnot(parInd <= length(partition$patient_40))
    
    input_data <- bortezomib$AUC_combined.sva
    input_label <- bortezomib.labels$AUC_combined
    input_partition = partition$patient_40[[parInd]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50_p40") {
    stopifnot(parInd <= length(partition$patient_40))
    
    input_data <- bortezomib$IC50_combined.sva
    input_label <- bortezomib.labels$IC50_combined
    input_partition = partition$patient_40[[parInd]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_ic50_p40") {
    stopifnot(parInd <= length(partition$patient_40))
    
    train_once = TRUE
    input_data <- bortezomib$IC50_combined.sva
    input_label <- bortezomib.labels$IC50_combined
    input_partition = partition$patient_40[[parInd]]$c2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_auc_p40") {
    stopifnot(parInd <= length(partition$patient_40))
    
    train_once = TRUE
    input_data <- bortezomib$AUC_combined.sva
    input_label <- bortezomib.labels$AUC_combined
    input_partition = partition$patient_40[[parInd]]$c2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_slope_p40") {
    stopifnot(parInd <= length(partition$patient_40))
    
    train_once = TRUE
    input_data <- bortezomib$slope_combined.sva
    input_label <- bortezomib.labels$slope_combined
    input_partition = partition$patient_40[[parInd]]$c2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else {
    stop(paste("args[5]", args[5], "is invalid."))
  }
  rm(bortezomib, bortezomib.labels)
} else if (args[4] == "bortezomib_no_homo") {
  load("Bortezomib/WS/bortezomib_data.RData")
  
  input.type_measure = "auc"
  input_snf.parameter <- seq(from = 5, to = 30, by = 5)
  
  if (args[5] == "p2p") {
    stopifnot(parInd <= length(partition$cell_lines_all))
    
    input_data <- bortezomib$patient.combat
    input_label <- bortezomib.labels$patient
    input_partition <- partition$cell_lines_all[[parInd]]$p2p
    input_feature.l1000 <- feature.l1000$pp
  } else if (args[5] == "cp2p_slope") {
    stopifnot(parInd <= length(partition$cell_lines_all))
    
    input_data <- bortezomib$slope_combined
    input_label <- bortezomib.labels$slope_combined
    input_partition = partition$cell_lines_all[[parInd]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc") {
    stopifnot(parInd <= length(partition$cell_lines_all))
    
    input_data <- bortezomib$AUC_combined
    input_label <- bortezomib.labels$AUC_combined
    input_partition = partition$cell_lines_all[[parInd]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50") {
    stopifnot(parInd <= length(partition$cell_lines_all))
    
    input_data <- bortezomib$IC50_combined
    input_label <- bortezomib.labels$IC50_combined
    input_partition = partition$cell_lines_all[[parInd]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50_p20") {
    stopifnot(parInd <= length(partition$patient_20_no_homogenization))
    
    input_data <- bortezomib$IC50_combined
    input_label <- bortezomib.labels$IC50_combined
    input_partition = partition$patient_20_no_homogenization[[parInd]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_p20") {
    stopifnot(parInd <= length(partition$patient_20_no_homogenization))
    
    input_data <- bortezomib$AUC_combined
    input_label <- bortezomib.labels$AUC_combined
    input_partition = partition$patient_20_no_homogenization[[parInd]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_slope_p20") {
    stopifnot(parInd <= length(partition$patient_20_no_homogenization))
    
    input_data <- bortezomib$slope_combined
    input_label <- bortezomib.labels$slope_combined
    input_partition = partition$patient_20_no_homogenization[[parInd]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_ic50_p20") {
    stopifnot(parInd <= length(partition$patient_20))
    
    train_once = TRUE
    input_data <- bortezomib$IC50_combined
    input_label <- bortezomib.labels$IC50_combined
    input_partition = partition$patient_20_no_homogenization[[parInd]]$c2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_auc_p20") {
    stopifnot(parInd <= length(partition$patient_20))
    
    train_once = TRUE
    input_data <- bortezomib$AUC_combined
    input_label <- bortezomib.labels$AUC_combined
    input_partition = partition$patient_20_no_homogenization[[parInd]]$c2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_slope_p20") {
    stopifnot(parInd <= length(partition$patient_20_no_homogenization))
    
    train_once = TRUE
    input_data <- bortezomib$slope_combined
    input_label <- bortezomib.labels$slope_combined
    input_partition = partition$patient_20_no_homogenization[[parInd]]$c2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_slope_p40") {
    stopifnot(parInd <= length(partition$patient_40_no_homogenization))
    
    input_data <- bortezomib$slope_combined
    input_label <- bortezomib.labels$slope_combined
    input_partition = partition$patient_40_no_homogenization[[parInd]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_p40") {
    stopifnot(parInd <= length(partition$patient_40_no_homogenization))
    
    input_data <- bortezomib$AUC_combined
    input_label <- bortezomib.labels$AUC_combined
    input_partition = partition$patient_40_no_homogenization[[parInd]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50_p40") {
    stopifnot(parInd <= length(partition$patient_40_no_homogenization))
    
    input_data <- bortezomib$IC50_combined
    input_label <- bortezomib.labels$IC50_combined
    input_partition = partition$patient_40_no_homogenization[[parInd]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_ic50_p40") {
    stopifnot(parInd <= length(partition$patient_40_no_homogenization))
    
    train_once = TRUE
    input_data <- bortezomib$IC50_combined
    input_label <- bortezomib.labels$IC50_combined
    input_partition = partition$patient_40_no_homogenization[[parInd]]$c2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_auc_p40") {
    stopifnot(parInd <= length(partition$patient_40_no_homogenization))
    
    train_once = TRUE
    input_data <- bortezomib$AUC_combined
    input_label <- bortezomib.labels$AUC_combined
    input_partition = partition$patient_40_no_homogenization[[parInd]]$c2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_slope_p40") {
    stopifnot(parInd <= length(partition$patient_40_no_homogenization))
    
    train_once = TRUE
    input_data <- bortezomib$slope_combined
    input_label <- bortezomib.labels$slope_combined
    input_partition = partition$patient_40_no_homogenization[[parInd]]$c2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_slope_p80") {
    stopifnot(parInd <= length(partition$patient_80_no_homogenization))
    
    input_data <- bortezomib$slope_combined
    input_label <- bortezomib.labels$slope_combined
    input_partition = partition$patient_80_no_homogenization[[parInd]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_p80") {
    stopifnot(parInd <= length(partition$patient_80_no_homogenization))
    
    input_data <- bortezomib$AUC_combined
    input_label <- bortezomib.labels$AUC_combined
    input_partition = partition$patient_80_no_homogenization[[parInd]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50_p80") {
    stopifnot(parInd <= length(partition$patient_80_no_homogenization))
    
    input_data <- bortezomib$IC50_combined
    input_label <- bortezomib.labels$IC50_combined
    input_partition = partition$patient_80_no_homogenization[[parInd]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_ic50_p80") {
    stopifnot(parInd <= length(partition$patient_80_no_homogenization))
    
    train_once = TRUE
    input_data <- bortezomib$IC50_combined
    input_label <- bortezomib.labels$IC50_combined
    input_partition = partition$patient_80_no_homogenization[[parInd]]$c2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_auc_p80") {
    stopifnot(parInd <= length(partition$patient_80_no_homogenization))
    
    train_once = TRUE
    input_data <- bortezomib$AUC_combined
    input_label <- bortezomib.labels$AUC_combined
    input_partition = partition$patient_80_no_homogenization[[parInd]]$c2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_slope_p80") {
    stopifnot(parInd <= length(partition$patient_80_no_homogenization))
    
    train_once = TRUE
    input_data <- bortezomib$slope_combined
    input_label <- bortezomib.labels$slope_combined
    input_partition = partition$patient_80_no_homogenization[[parInd]]$c2p.slope
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
    stopifnot(parInd <= length(partition$cell_lines_all))
    
    input_data <- docetaxel$patient
    input_label <- docetaxel.labels$patient
    input_partition <- partition$cell_lines_all[[parInd]]$p2p
    input_feature.l1000 <- feature.l1000$pp
    input_snf.parameter <- seq(from = 5, to = length(input_partition[[1]]$training_index), by = 3)
  } else if (args[5] == "cp2p_slope") {
    stopifnot(parInd <= length(partition$cell_lines_all))
    
    input_data <- docetaxel$slope_combined 
    input_label <- docetaxel.labels$slope_combined
    input_partition = partition$cell_lines_all[[parInd]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc") {
    stopifnot(parInd <= length(partition$cell_lines_all))
    
    input_data <- docetaxel$AUC_combined
    input_label <- docetaxel.labels$AUC_combined
    input_partition = partition$cell_lines_all[[parInd]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50") {
    stopifnot(parInd <= length(partition$cell_lines_all))
    
    input_data <- docetaxel$IC50_combined
    input_label <- docetaxel.labels$IC50_combined
    input_partition = partition$cell_lines_all[[parInd]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_slope_breast") {
    stopifnot(parInd <= length(partition$cell_lines_breast))
    
    input_data <- docetaxel$slope_breast 
    input_label <- docetaxel.labels$slope_breast.combined
    input_partition = partition$cell_lines_breast[[parInd]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_breast") {
    stopifnot(parInd <= length(partition$cell_lines_breast))
    
    input_data <- docetaxel$AUC_breast
    input_label <- docetaxel.labels$AUC_breast.combined
    input_partition = partition$cell_lines_breast[[parInd]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50_breast") {
    stopifnot(parInd <= length(partition$cell_lines_breast))
    
    input_data <- docetaxel$IC50_breast
    input_label <- docetaxel.labels$IC50_breast.combined
    input_partition = partition$cell_lines_breast[[parInd]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50_p23") {
    stopifnot(parInd <= length(partition$patient_23))
    
    input_data <- docetaxel$IC50_combined
    input_label <- docetaxel.labels$IC50_combined
    input_partition = partition$patient_23[[parInd]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_p23") {
    stopifnot(parInd <= length(partition$patient_23))
    
    input_data <- docetaxel$AUC_combined
    input_label <- docetaxel.labels$AUC_combined
    input_partition = partition$patient_23[[parInd]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_slope_p23") {
    stopifnot(parInd <= length(partition$patient_23))
    
    input_data <- docetaxel$slope_combined
    input_label <- docetaxel.labels$slope_combined
    input_partition = partition$patient_23[[parInd]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_ic50_p23") {
    stopifnot(parInd <= length(partition$patient_23))
    
    train_once = TRUE
    input_data <- docetaxel$IC50_combined
    input_label <- docetaxel.labels$IC50_combined
    input_partition = partition$patient_23[[parInd]]$c2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_auc_p23") {
    stopifnot(parInd <= length(partition$patient_23))
    
    train_once = TRUE
    input_data <- docetaxel$AUC_combined
    input_label <- docetaxel.labels$AUC_combined
    input_partition = partition$patient_23[[parInd]]$c2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_slope_p23") {
    stopifnot(parInd <= length(partition$patient_23))
    
    train_once = TRUE
    input_data <- docetaxel$slope_combined
    input_label <- docetaxel.labels$slope_combined
    input_partition = partition$patient_23[[parInd]]$c2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else {
    stop(paste("args[5]", args[5], "is invalid."))
  }
  rm(docetaxel, docetaxel.labels)
} else if (args[4] == "erlotinib_gdsc" || args[4] == "erlotinib_ccle" || args[4] == "erlotinib_all") {
  
  if (args[4] == "erlotinib_gdsc") {
    load("Erlotinib/WS/erlotinib_homogenized_data_gdsc.RData")
  } else if (args[4] == "erlotinib_ccle") {
    load("Erlotinib/WS/erlotinib_homogenized_data_ccle.RData")
  } else if (args[4] == "erlotinib_all") {
    load("Erlotinib/WS/erlotinib_homogenized_data_all.RData")
  }
  
  input.type_measure = "acc"
  input_snf.parameter <- seq(from = 5, to = 30, by = 5)
  
  if (args[5] == "p2p") {
    stopifnot(parInd <= length(partition$cell_lines_all))
    
    input_data <- erlotinib$patient
    input_label <- erlotinib.labels$patient
    input_partition <- partition$cell_lines_all[[parInd]]$p2p
    input_feature.l1000 <- feature.l1000$pp
    input_snf.parameter <- seq(from = 5, to = length(input_partition[[1]]$training_index), by = 3)
  } else if (args[5] == "cp2p_slope") {
    stopifnot(parInd <= length(partition$cell_lines_all))
    
    input_data <- erlotinib$slope_combined.ComBat
    input_label <- erlotinib.labels$slope_combined
    input_partition = partition$cell_lines_all[[parInd]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc") {
    stopifnot(parInd <= length(partition$cell_lines_all))
    
    input_data <- erlotinib$AUC_combined.ComBat
    input_label <- erlotinib.labels$AUC_combined
    input_partition = partition$cell_lines_all[[parInd]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50") {
    stopifnot(parInd <= length(partition$cell_lines_all))
    
    input_data <- erlotinib$IC50_combined.ComBat
    input_label <- erlotinib.labels$IC50_combined
    input_partition = partition$cell_lines_all[[parInd]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_slope_lung") {
    stopifnot(parInd <= length(partition$cell_lines_lung))
    
    input_data <- erlotinib$slope_lung.ComBat 
    input_label <- erlotinib.labels$slope_lung
    input_partition = partition$cell_lines_lung[[parInd]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_lung") {
    stopifnot(parInd <= length(partition$cell_lines_lung))
    
    input_data <- erlotinib$AUC_lung.ComBat
    input_label <- erlotinib.labels$AUC_lung
    input_partition = partition$cell_lines_lung[[parInd]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50_lung") {
    stopifnot(parInd <= length(partition$cell_lines_lung))
    
    input_data <- erlotinib$IC50_lung.ComBat
    input_label <- erlotinib.labels$IC50_lung
    input_partition = partition$cell_lines_lung[[parInd]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_ic50_p24") {
    stopifnot(parInd <= length(partition$patient_24))
    
    input_data <- erlotinib$IC50_combined.ComBat
    input_label <- erlotinib.labels$IC50_combined
    input_partition = partition$patient_24[[parInd]]$cp2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_p24") {
    stopifnot(parInd <= length(partition$patient_24))
    
    input_data <- erlotinib$AUC_combined.ComBat
    input_label <- erlotinib.labels$AUC_combined
    input_partition = partition$patient_24[[parInd]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_slope_p24") {
    stopifnot(parInd <= length(partition$patient_24))
    
    input_data <- erlotinib$slope_combined.ComBat
    input_label <- erlotinib.labels$slope_combined
    input_partition = partition$patient_24[[parInd]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_ic50_p24") {
    stopifnot(parInd <= length(partition$patient_24))
    
    train_once = TRUE
    input_data <- erlotinib$IC50_combined.ComBat
    input_label <- erlotinib.labels$IC50_combined
    input_partition = partition$patient_24[[parInd]]$c2p.IC50
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_auc_p24") {
    stopifnot(parInd <= length(partition$patient_24))
    
    train_once = TRUE
    input_data <- erlotinib$AUC_combined.ComBat
    input_label <- erlotinib.labels$AUC_combined
    input_partition = partition$patient_24[[parInd]]$c2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_slope_p24") {
    stopifnot(parInd <= length(partition$patient_24))
    
    train_once = TRUE
    input_data <- erlotinib$slope_combined.ComBat
    input_label <- erlotinib.labels$slope_combined
    input_partition = partition$patient_24[[parInd]]$c2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else {
    stop(paste("args[5]", args[5], "is invalid."))
  }
  rm(erlotinib, erlotinib.labels)
} else if (args[4] == "epirubicin" 
           || args[4] == "epirubicin_ComBat_no_label"
           || args[4] == "epirubicin_ComBat") {
  
  if (args[4] == "epirubicin") {
    load("Epirubicin/WS/epirubicin_data_ComBat.RData")
    epirubicin$slope_combined.ComBat <- epirubicin$slope_combined
    epirubicin$AUC_combined.ComBat <- epirubicin$AUC_combined
    
  } else if (args[4] == "epirubicin_ComBat_no_label") {
    load("Epirubicin/WS/epirubicin_data_ComBat_no_labels.RData")
  } else if (args[4] == "epirubicin_ComBat") {
    load("Epirubicin/WS/epirubicin_data_ComBat.RData")
  } else {
    stop("args[4] is wrong!")
  }
  
  input.type_measure = "auc"
  input_snf.parameter <- seq(from = 5, to = 30, by = 5)
  
  if (args[5] == "p2p") {
    stopifnot(parInd <= length(partition$cell_lines_all))
    
    input_data <- epirubicin$patient
    input_label <- epirubicin.labels$patient
    input_partition <- partition$cell_lines_all[[parInd]]$p2p
    input_feature.l1000 <- feature.l1000$pp
  } else if (args[5] == "cp2p_slope") {
    stopifnot(parInd <= length(partition$cell_lines_all))
    
    input_data <- epirubicin$slope_combined 
    input_label <- epirubicin.labels$slope_combined
    input_partition = partition$cell_lines_all[[parInd]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc") {
    stopifnot(parInd <= length(partition$cell_lines_all))
    
    input_data <- epirubicin$AUC_combined
    input_label <- epirubicin.labels$AUC_combined
    input_partition = partition$cell_lines_all[[parInd]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_auc_p20") {
    stopifnot(parInd <= length(partition$patient_20))
    
    input_data <- epirubicin$AUC_combined
    input_label <- epirubicin.labels$AUC_combined
    input_partition = partition$patient_20[[parInd]]$cp2p.AUC
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "cp2p_slope_p20") {
    stopifnot(parInd <= length(partition$patient_20))
    
    input_data <- epirubicin$slope_combined
    input_label <- epirubicin.labels$slope_combined
    input_partition = partition$patient_20[[parInd]]$cp2p.slope
    input_feature.l1000 <- feature.l1000$cp
  } else if (args[5] == "c2p_auc_p20") {
    stopifnot(parInd <= length(partition$patient_20))
    
    train_once = TRUE
    input_data <- epirubicin$AUC_combined
    input_label <- epirubicin.labels$AUC_combined
    input_partition = partition$patient_20[[parInd]]$c2p.AUC
    input_feature.l1000 <- feature.l1000$cp
    input.type_measure = "acc"
    input.type_measure.test = "auc"
    INPUT.NFOLDS = 3 # needed or else the training fails
  } else if (args[5] == "c2p_slope_p20") {
    stopifnot(parInd <= length(partition$patient_20))
    
    train_once = TRUE
    input_data <- epirubicin$slope_combined
    input_label <- epirubicin.labels$slope_combined
    input_partition = partition$patient_20[[parInd]]$c2p.slope
    input_feature.l1000 <- feature.l1000$cp
    input.type_measure = "acc"
    input.type_measure.test = "auc"
    INPUT.NFOLDS = 3 # needed or else the training fails
  } else if (args[5] == "c2p_auc_p40") {
    stopifnot(parInd <= length(partition$patient_40))
    
    train_once = TRUE
    input_data <- epirubicin$AUC_combined
    input_label <- epirubicin.labels$AUC_combined
    input_partition = partition$patient_40[[parInd]]$c2p.AUC
    input_feature.l1000 <- feature.l1000$cp
    input.type_measure = "acc"
    input.type_measure.test = "auc"
    INPUT.NFOLDS = 3 # needed or else the training fails
  } else if (args[5] == "c2p_slope_p40") {
    stopifnot(parInd <= length(partition$patient_40))
    
    train_once = TRUE
    input_data <- epirubicin$slope_combined
    input_label <- epirubicin.labels$slope_combined
    input_partition = partition$patient_40[[parInd]]$c2p.slope
    input_feature.l1000 <- feature.l1000$cp
    input.type_measure = "acc"
    input.type_measure.test = "auc"
    INPUT.NFOLDS = 3 # needed or else the training fails
  } else if (args[5] == "c2p_slope_once") {
    stopifnot(parInd <= length(partition$patient_40))
    
    train_once = TRUE
    input_data <- epirubicin$slope_combined
    input_label <- epirubicin.labels$slope_combined
    input_partition = partition$patient_40[[parInd]]$c2p.slope
    input_feature.l1000 <- feature.l1000$cp
    
    input_partition[[1]]$training_index = 
      which(epirubicin.labels$slope_combined.source != "patient")
    input_partition[[1]]$test_index = 
      which(epirubicin.labels$slope_combined.source == "patient")
    
    INPUT.NFOLDS = 3 # needed or else the training fails
  } else if (args[5] == "c2p_auc_once") {
    stopifnot(parInd <= length(partition$patient_40))
    
    train_once = TRUE
    input_data <- epirubicin$auc_combined
    input_label <- epirubicin.labels$auc_combined
    input_partition = partition$patient_40[[parInd]]$c2p.auc
    input_feature.l1000 <- feature.l1000$cp
    
    input_partition[[1]]$training_index = 
      which(epirubicin.labels$auc_combined.source != "patient")
    input_partition[[1]]$test_index = 
      which(epirubicin.labels$auc_combined.source == "patient")
    
    INPUT.NFOLDS = 3 # needed or else the training fails
  } else {
    stop(paste("args[5]", args[5], "is invalid."))
  }
  rm(epirubicin, epirubicin.labels)
} else {
  stop(paste("args[4]", args[4], "is invalid."))
}

if (is.null(input.type_measure.test)) {
  input.type_measure.test = input.type_measure
}

stopifnot(!is.null(input_data))
stopifnot(!is.null(input_label))
stopifnot(!is.null(input_partition))

rm(partition, feature.l1000)

# do the computation
if (train_once == TRUE) {
  print("train_once = TRUE!!")
}
print(paste0(args[4], "_", args[5], "_", PARTITION_BEGIN, "to", PARTITION_END, "_parInd", parInd, ".RData"))
source("Common/train_and_predict.R")
##

save.image(paste0(args[4], "_", args[5], "_", PARTITION_BEGIN, "to", PARTITION_END, "_parInd", parInd, ".RData"))
