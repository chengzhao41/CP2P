# Load Libraries ----------------------------------------------------------
library("doParallel")
library("sva")

source('Common/preparing_data_helper.R')
source('Common/drug_cut/callingWaterfall.R')
source('Common/drug_cut/distancePointLine.R')
source('Common/drug_cut/distancePointSegment.R')
source('Common/comGENE.R')
source("Common/generate_random_partition.R")

# Load CGP Cell line data --------------------------------------------
docetaxel.labels <- list()
load("CGP/cdrug2_cgp_ccle_all.RData")

# Get Slope Summary Statistics --------------------------------------------
cgp_sensitivity_1 <- read.csv("CGP/cgp_sensitivity_1.csv")

temp.docetaxel_ind <- which(cgp_sensitivity_1$drug.name == "Docetaxel")
cgp_sensitivity_1 <- cgp_sensitivity_1[temp.docetaxel_ind, ]
dim(cgp_sensitivity_1)

temp.na_ind <- which(!is.na(cgp_sensitivity_1$slope0.sensitivity.call))
length(temp.na_ind)
temp.response_slope <- cgp_sensitivity_1$slope0.sensitivity.call[temp.na_ind]
length(temp.response_slope)
stopifnot(sum(is.na(temp.response_slope)) == 0)
temp.labels_slope <- temp.response_slope == 1

names(temp.labels_slope) <- cgp_sensitivity_1$cellline[temp.na_ind]
temp.ind <- rownames(data.ge.cgp) %in% names(temp.labels_slope)
temp.slope_ind <- match(rownames(data.ge.cgp)[temp.ind], names(temp.labels_slope))
stopifnot(names(temp.labels_slope)[temp.slope_ind] == rownames(data.ge.cgp)[temp.ind])
stopifnot(sum(is.na(temp.labels_slope)) == 0)

docetaxel.labels$slope <- temp.labels_slope[temp.slope_ind]
docetaxel.labels$slope_ind <- which(rownames(data.ge.cgp) %in% names(temp.labels_slope))
length(docetaxel.labels$slope_ind)

stopifnot(names(docetaxel.labels$slope) == rownames(data.ge.cgp)[docetaxel.labels$slope_ind])

# Get IC50 and AUC Summary Statistics -------------------------------------
temp.drug_ind <- which(druginfo.cgp$drug.name == "DOCETAXEL")

temp.ic50 <- drugpheno.cgp$IC50[, temp.drug_ind]
temp.ic50_ind <- which(!is.na(temp.ic50))
temp.ic50 <- temp.ic50[temp.ic50_ind]
length(temp.ic50)

temp.auc <- drugpheno.cgp$AUC[, temp.drug_ind]
temp.auc_ind <- which(!is.na(temp.auc))
temp.auc <- temp.auc[temp.auc_ind]
length(temp.auc)

temp.response <- callingWaterfall(temp.ic50, type="IC50")
temp.response_auc <- callingWaterfall(temp.auc, type="AUC")
table(temp.response)

temp.label_IC50 <- temp.response != "resistant"
names(temp.label_IC50) <- names(temp.response)

temp.label_auc <- temp.response_auc != "resistant"
names(temp.label_auc) <- names(temp.response_auc)

# check agreement between AUC and IC50
stopifnot(names(temp.label_auc) == names(temp.label_IC50))
print(paste(sum(temp.label_auc == temp.label_IC50), "/", length(temp.label_IC50)))
## 548 / 663 for AUC and IC50

# Getting the docetaxel labels from IC50 ----------------------------------
docetaxel.labels$IC50_ind <- which(rownames(data.ge.cgp) %in% names(temp.label_IC50))
docetaxel.labels$AUC_ind <- which(rownames(data.ge.cgp) %in% names(temp.label_auc))

stopifnot(names(temp.label_IC50) == rownames(data.ge.cgp)[docetaxel.labels$IC50_ind])
docetaxel.labels$IC50 <- temp.label_IC50
stopifnot(names(temp.label_auc) == rownames(data.ge.cgp)[docetaxel.labels$AUC_ind])
docetaxel.labels$AUC <- temp.label_auc

## check how much of the labels agree with each other for IC50 and slope
temp.ind <- names(docetaxel.labels$IC50) %in% names(docetaxel.labels$slope)
stopifnot(sum(temp.ind == FALSE) == length(docetaxel.labels$IC50) - length(docetaxel.labels$slope))
length(temp.ind)
length(docetaxel.labels$IC50)
length(docetaxel.labels$IC50[temp.ind])
length(docetaxel.labels$slope)
temp.slope_include <- which(!is.na(docetaxel.labels$slope))
stopifnot(names(docetaxel.labels$slope[temp.slope_include]) == names(docetaxel.labels$IC50[temp.ind]))

print(paste(sum(docetaxel.labels$slope[temp.slope_include] == docetaxel.labels$IC50[temp.ind]), 
            "/", length(docetaxel.labels$IC50[temp.ind])))
# 458 / 650 in with IC50 and slope

## check how much of the labels agree with each other for IC50 and slope
temp.ind <- names(docetaxel.labels$AUC) %in% names(docetaxel.labels$slope)
stopifnot(sum(temp.ind == FALSE) == length(docetaxel.labels$AUC) - length(docetaxel.labels$slope))
length(temp.ind)
length(docetaxel.labels$AUC)
length(docetaxel.labels$AUC[temp.ind])
length(docetaxel.labels$slope)
temp.slope_include <- which(!is.na(docetaxel.labels$slope))

stopifnot(names(docetaxel.labels$slope[temp.slope_include]) == names(docetaxel.labels$AUC[temp.ind]))

print(paste(sum(docetaxel.labels$slope[temp.slope_include] == docetaxel.labels$AUC[temp.ind]),
      "/", length(docetaxel.labels$AUC[temp.ind])))
# 567 / 650 with AUC and slope

# get the data for slopes -------------------------------------------------
docetaxel <- list()
docetaxel$cgp_slope <- data.ge.cgp[docetaxel.labels$slope_ind,  ]
docetaxel$cgp_IC50 <- data.ge.cgp[docetaxel.labels$IC50_ind,  ]
docetaxel$cgp_AUC <- data.ge.cgp[docetaxel.labels$AUC_ind,  ]

temp <- remove_features(input_data = docetaxel$cgp_slope, input_ind = docetaxel.labels$slope_ind,
                        input_label = docetaxel.labels$slope)
docetaxel$cgp_slope <- temp$data
docetaxel.labels$slope_ind <- temp$ind
docetaxel.labels$slope <- temp$label

rm(temp)
temp <- remove_features(input_data = docetaxel$cgp_IC50, input_ind = docetaxel.labels$IC50_ind,
                        input_label = docetaxel.labels$IC50)
docetaxel$cgp_IC50 <- temp$data
docetaxel.labels$IC50_ind <- temp$ind
docetaxel.labels$IC50 <- temp$label

rm(temp)
temp <- remove_features(input_data = docetaxel$cgp_AUC, input_ind = docetaxel.labels$AUC_ind,
                        input_label = docetaxel.labels$AUC)
docetaxel$cgp_AUC <- temp$data
docetaxel.labels$AUC_ind <- temp$ind
docetaxel.labels$AUC <- temp$label

## check agreement
temp.ind <- match(names(docetaxel.labels$slope), names(docetaxel.labels$AUC))
stopifnot(names(docetaxel.labels$AUC)[temp.ind] == names(docetaxel.labels$slope))

print(paste(sum(docetaxel.labels$AUC[temp.ind] == docetaxel.labels$slope), "/", sum(!is.na(temp.ind))))
# 528 / 605 for slope vs. AUC

temp.ind <- match(names(docetaxel.labels$slope), names(docetaxel.labels$IC50))
print(paste(sum(docetaxel.labels$IC50[temp.ind] == docetaxel.labels$slope), "/", sum(!is.na(temp.ind))))
# 427 / 605 for slope vs. IC50

temp.ind <- match(names(docetaxel.labels$AUC), names(docetaxel.labels$IC50))
print(paste(sum(docetaxel.labels$IC50[temp.ind] == docetaxel.labels$AUC), "/", sum(!is.na(temp.ind))))
# 511 / 618 for AUC vs. IC50 

# get patient data --------------------------------------------------------
load("Docetaxel/WS/pp.RData")
docetaxel$patient <- docetaxel.patient
docetaxel.labels$patient <- pp.ground_truth == 1
names(docetaxel.labels$patient) <- rownames(docetaxel.patient)

show_pca(input_data = docetaxel$cgp_slope, label = docetaxel.labels$slope)
show_pca(input_data = docetaxel$cgp_AUC, label = docetaxel.labels$AUC)
show_pca(input_data = docetaxel$cgp_IC50, label = docetaxel.labels$IC50)

# Using sva to harmonize patients and cell lines - Slopes -----------------
docetaxel.labels$slope_combined <- c(docetaxel.labels$patient, docetaxel.labels$slope)
stopifnot(substring(colnames(docetaxel$cgp_slope), 8) == annot.ge.cgp$EntrezGene.ID)
colnames(docetaxel$cgp_slope) <- annot.ge.cgp$symbol

temp.data <- comGENE(docetaxel$patient, scale(docetaxel$cgp_slope))
mean(temp.data[[1]]) #-3.25391e-18
mean(temp.data[[2]]) #1.134142e-18

docetaxel.labels$slope_combined.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("cgp", dim(temp.data[[2]])[1]))
docetaxel$combined_slope.ComBat <- ComBat_combine(batch = docetaxel.labels$slope_combined.source,
                                                  label = docetaxel.labels$slope_combined,
                                                  input_data = rbind(temp.data[[1]], temp.data[[2]]))
mean(docetaxel$combined_slope.ComBat) #-1.32049e-07
show_pca(input_data = docetaxel$combined_slope.ComBat, label = docetaxel.labels$slope_combined)
show_pca(input_data = docetaxel$combined_slope.ComBat, 
         label = c(rep("patient", dim(temp.data[[1]])[1]), 
                   sampleinfo.cgp$tissue.type[docetaxel.labels$slope_ind]))
rm(temp.data)

# Using sva to harmonize patients and cell lines - IC50 -----------------
docetaxel.labels$IC50_combined <- c(docetaxel.labels$patient, docetaxel.labels$IC50)
stopifnot(substring(colnames(docetaxel.labels$IC50), 8) == annot.ge.cgp$EntrezGene.ID)
colnames(docetaxel$cgp_IC50) <- annot.ge.cgp$symbol

temp.data <- comGENE(docetaxel$patient, scale(docetaxel$cgp_IC50))
mean(temp.data[[1]]) #-3.25391e-18
mean(temp.data[[2]]) #-3.049527e-18

docetaxel.labels$IC50_combined.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("cgp", dim(temp.data[[2]])[1]))
docetaxel$combined_IC50 <- rbind(temp.data[[1]], temp.data[[2]])
show_pca(input_data = docetaxel$combined_IC50, label = c(rep("patient", dim(temp.data[[1]])[1]), rep("cgp", dim(temp.data[[2]])[1])))

docetaxel$combined_IC50.ComBat <- ComBat_combine(batch = docetaxel.labels$IC50_combined.source,
                                                 label = docetaxel.labels$IC50_combined,
                                                 input_data = docetaxel$combined_IC50)
mean(docetaxel$combined_IC50.ComBat) #-9.550229e-07
show_pca(input_data = docetaxel$combined_IC50.ComBat, label = docetaxel.labels$IC50_combined)
show_pca(input_data = docetaxel$combined_IC50.ComBat, label = docetaxel.labels$IC50_combined.source)

rm(temp.data)

# Using sva to harmonize patients and cell lines - AUC -----------------
docetaxel.labels$AUC_combined <- c(docetaxel.labels$patient, docetaxel.labels$AUC)
stopifnot(substring(colnames(docetaxel.labels$AUC), 8) == annot.ge.cgp$EntrezGene.ID)
colnames(docetaxel$cgp_AUC) <- annot.ge.cgp$symbol

temp.data <- comGENE(docetaxel$patient, scale(docetaxel$cgp_AUC))
mean(temp.data[[1]]) #-3.049506e-18
mean(temp.data[[2]]) #-3.253883e-18

docetaxel.labels$AUC_combined.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("cgp", dim(temp.data[[2]])[1]))
docetaxel$combined_AUC <- rbind(temp.data[[1]], temp.data[[2]])

docetaxel$combined_AUC.ComBat <- ComBat_combine(batch = docetaxel.labels$AUC_combined.source,
                                                label = docetaxel.labels$AUC_combined, 
                                                input_data = docetaxel$combined_AUC)
mean(docetaxel$combined_AUC.ComBat) #-5.052369e-07
show_pca(input_data = docetaxel$combined_AUC.ComBat, label = docetaxel.labels$AUC_combined)
show_pca(input_data = docetaxel$combined_AUC.ComBat, label = docetaxel.labels$AUC_combined.source)

rm(temp.data)

# Breast Cell lines only --------------------------------------------------
## IC50
temp.ind <- which(sampleinfo.cgp$tissue.type[docetaxel.labels$IC50_ind] == "breast")
temp.data <- comGENE(docetaxel$patient, scale(docetaxel$cgp_IC50[temp.ind, ]))
mean(temp.data[[1]]) #-3.25391e-18
mean(temp.data[[2]]) #-6.222348e-18

docetaxel.labels$IC50_breast.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("cgp", dim(temp.data[[2]])[1]))
docetaxel.labels$IC50_breast_only <- docetaxel.labels$IC50[temp.ind]
docetaxel.labels$IC50_breast <- c(docetaxel.labels$patient, docetaxel.labels$IC50[temp.ind])

docetaxel$breast_IC50.ComBat <- ComBat_combine(batch = docetaxel.labels$IC50_breast.source,
                                               label = docetaxel.labels$IC50_breast,
                                               input_data = rbind(temp.data[[1]], temp.data[[2]]))
mean(docetaxel$breast_IC50.ComBat) #-1.961629e-05
show_pca(input_data = docetaxel$breast_IC50.ComBat, label = docetaxel.labels$IC50_breast)
show_pca(input_data = docetaxel$breast_IC50.ComBat, label = docetaxel.labels$IC50_breast.source)
rm(temp.data)
getwd()

# AUC
temp.ind <- which(sampleinfo.cgp$tissue.type[docetaxel.labels$AUC_ind] == "breast")
temp.data <- comGENE(docetaxel$patient, scale(docetaxel$cgp_AUC[temp.ind, ]))
mean(temp.data[[1]]) #-3.253883e-18
mean(temp.data[[2]]) #-6.222309e-18

docetaxel.labels$AUC_breast.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("cgp", dim(temp.data[[2]])[1]))
docetaxel.labels$AUC_breast_only <- docetaxel.labels$AUC[temp.ind]
docetaxel.labels$AUC_breast <- c(docetaxel.labels$patient, docetaxel.labels$AUC[temp.ind])

docetaxel$breast_AUC.ComBat <- ComBat_combine(batch = docetaxel.labels$AUC_breast.source,
                                              label = docetaxel.labels$AUC_breast,
                                              input_data = rbind(temp.data[[1]], temp.data[[2]]))
mean(docetaxel$breast_AUC.ComBat) # -1.174055e-05
show_pca(input_data = docetaxel$breast_AUC.ComBat, label = docetaxel.labels$AUC_breast)
show_pca(input_data = docetaxel$breast_AUC.ComBat, label = docetaxel.labels$AUC_breast.source)
rm(temp.data)
getwd()

# Slope
temp.ind <- which(sampleinfo.cgp$tissue.type[docetaxel.labels$slope_ind] == "breast")
temp.data <- comGENE(docetaxel$patient, scale(docetaxel$cgp_slope[temp.ind, ]))
mean(temp.data[[1]]) #-3.253883e-18
mean(temp.data[[2]]) #-1.43896e-17

docetaxel.labels$slope_breast.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("cgp", dim(temp.data[[2]])[1]))
docetaxel.labels$slope_breast_only <- docetaxel.labels$slope[temp.ind]
docetaxel.labels$slope_breast <- c(docetaxel.labels$patient, docetaxel.labels$slope[temp.ind])

docetaxel$breast_slope.ComBat <- ComBat_combine(batch = docetaxel.labels$slope_breast.source,
                                                label = docetaxel.labels$slope_breast,
                                                input_data = rbind(temp.data[[1]], temp.data[[2]]))
mean(docetaxel$breast_slope.ComBat) #-7.523195e-06
show_pca(input_data = docetaxel$breast_slope.ComBat, label = docetaxel.labels$slope_breast)
show_pca(input_data = docetaxel$breast_slope.ComBat, label = docetaxel.labels$slope_breast.source)
rm(temp.data)

# get l1000 features ------------------------------------------------------
Landmark_Genes_n978 <- read.csv("Common/Landmark_Genes_n978.csv")
stopifnot(colnames(docetaxel$cgp_IC50) == colnames(docetaxel$cgp_AUC))
stopifnot(colnames(docetaxel$cgp_slope) == colnames(docetaxel$cgp_AUC))
feature.l1000 <- list()
feature.l1000$cp <- which(colnames(docetaxel$combined_slope.ComBat) %in% Landmark_Genes_n978$Gene.Symbol)
feature.l1000$pp <- which(colnames(docetaxel$patient) %in% Landmark_Genes_n978$Gene.Symbol)
length(feature.l1000$cp)
length(feature.l1000$pp)
getwd()

# find the common cell lines that exists for all 3 response labels --------
docetaxel.labels$IC50_ind_common = intersect(
  which(names(docetaxel.labels$IC50) %in% names(docetaxel.labels$AUC)),
  which(names(docetaxel.labels$IC50) %in% names(docetaxel.labels$slope)))
docetaxel.labels$AUC_ind_common = intersect(
  which(names(docetaxel.labels$AUC) %in% names(docetaxel.labels$IC50)),
  which(names(docetaxel.labels$AUC) %in% names(docetaxel.labels$slope)))
docetaxel.labels$slope_ind_common = intersect(
  which(names(docetaxel.labels$slope) %in% names(docetaxel.labels$IC50)),
  which(names(docetaxel.labels$slope) %in% names(docetaxel.labels$AUC)))

stopifnot(names(docetaxel.labels$IC50[docetaxel.labels$IC50_ind_common]) 
          == names(docetaxel.labels$AUC[docetaxel.labels$AUC_ind_common]))
stopifnot(names(docetaxel.labels$IC50[docetaxel.labels$IC50_ind_common]) 
          == names(docetaxel.labels$slope[docetaxel.labels$slope_ind_common]))

# for breast only
docetaxel.labels$IC50_breast_ind_common = intersect(
  which(names(docetaxel.labels$IC50_breast) %in% names(docetaxel.labels$AUC_breast)),
  which(names(docetaxel.labels$IC50_breast) %in% names(docetaxel.labels$slope_breast)))
docetaxel.labels$AUC_breast_ind_common = intersect(
  which(names(docetaxel.labels$AUC_breast) %in% names(docetaxel.labels$IC50_breast)),
  which(names(docetaxel.labels$AUC_breast) %in% names(docetaxel.labels$slope_breast)))
docetaxel.labels$slope_breast_ind_common = intersect(
  which(names(docetaxel.labels$slope_breast) %in% names(docetaxel.labels$IC50_breast)),
  which(names(docetaxel.labels$slope_breast) %in% names(docetaxel.labels$AUC_breast)))

stopifnot(names(docetaxel.labels$IC50_breast[docetaxel.labels$IC50_breast_ind_common]) 
          == names(docetaxel.labels$AUC_breast[docetaxel.labels$AUC_breast_ind_common]))
stopifnot(names(docetaxel.labels$IC50_breast[docetaxel.labels$IC50_breast_ind_common]) 
          == names(docetaxel.labels$slope_breast[docetaxel.labels$slope_breast_ind_common]))

# create partitions -------------------------------------------------------
partition <- list()

input.labels_cell_lines = list()
input.labels_cell_lines$slope = docetaxel.labels$slope[docetaxel.labels$slope_ind_common]
input.labels_cell_lines$IC50 = docetaxel.labels$IC50[docetaxel.labels$IC50_ind_common]
input.labels_cell_lines$AUC = docetaxel.labels$AUC[docetaxel.labels$AUC_ind_common]

input.cell_line_order = list()
input.cell_line_order$slope = 1:length(input.labels_cell_lines$slope)
input.cell_line_order$IC50 = 1:length(input.labels_cell_lines$IC50)
input.cell_line_order$AUC = 1:length(input.labels_cell_lines$AUC)

stopifnot(input.cell_line_order$slope == input.cell_line_order$IC50)  
stopifnot(input.cell_line_order$AUC == input.cell_line_order$IC50)

cell_lines_all <- foreach(input.training_amount.p = seq(from = 14, to = 23, by = 1)
                          , .errorhandling = "stop") %dopar% {
                generate_random_partition(labels_cell_lines = input.labels_cell_lines, 
                                          labels_patient = docetaxel.labels$patient,
                                          training_amount.p = input.training_amount.p,
                                          cell_line_order = input.cell_line_order,
                                          acc_training = TRUE,
                                          training_amount.c = length(input.cell_line_order$AUC))
                          }

partition$cell_lines_all <- cell_lines_all

# Breast Only -------------------------------------------------------------
input.labels_cell_lines = list()
input.labels_cell_lines$slope = docetaxel.labels$slope_breast[docetaxel.labels$slope_breast_ind_common]
input.labels_cell_lines$IC50 = docetaxel.labels$IC50_breast[docetaxel.labels$IC50_breast_ind_common]
input.labels_cell_lines$AUC = docetaxel.labels$AUC_breast[docetaxel.labels$AUC_breast_ind_common]

input.cell_line_order = list()
input.cell_line_order$slope = 1:length(input.labels_cell_lines$slope)
input.cell_line_order$IC50 = 1:length(input.labels_cell_lines$IC50)
input.cell_line_order$AUC = 1:length(input.labels_cell_lines$AUC)

stopifnot(input.cell_line_order$slope == input.cell_line_order$IC50)  
stopifnot(input.cell_line_order$AUC == input.cell_line_order$IC50)

cell_lines_breast <- foreach(input.training_amount.p = seq(from = 14, to = 23, by = 1)
                          , .errorhandling = "stop") %dopar% {
    generate_random_partition(labels_cell_lines = input.labels_cell_lines, 
                              labels_patient = docetaxel.labels$patient,
                              training_amount.p = input.training_amount.p,
                              cell_line_order = input.cell_line_order,
                              acc_training = TRUE,
                              training_amount.c = length(input.cell_line_order$AUC),
                              input_partition = cell_lines_all[[input.training_amount.p - 13]])
                          }

partition$cell_lines_breast <- cell_lines_breast

# using cell lines based on brca similarities -----------------------------
load("CGP/cosmic.tcga.RData")

input.labels_cell_lines = list()
input.labels_cell_lines$slope = docetaxel.labels$slope[docetaxel.labels$slope_ind_common]
input.labels_cell_lines$IC50 = docetaxel.labels$IC50[docetaxel.labels$IC50_ind_common]
input.labels_cell_lines$AUC = docetaxel.labels$AUC[docetaxel.labels$AUC_ind_common]

input.cell_line_order = list()
input.cell_line_order$slope <- order(match(names(input.labels_cell_lines$slope), brca_ordered$cell.line))
input.cell_line_order$IC50 <- order(match(names(input.labels_cell_lines$IC50), brca_ordered$cell.line))
input.cell_line_order$AUC <- order(match(names(input.labels_cell_lines$AUC), brca_ordered$cell.line))

stopifnot(!is.unsorted(match(brca_ordered$cell.line, 
                             names(input.labels_cell_lines$slope)[input.cell_line_order$slope] ), na.rm = TRUE))
stopifnot( length(table(input.labels_cell_lines$slope[input.cell_line_order$slope][1:30])) == 2 )
stopifnot( length(table(input.labels_cell_lines$AUC[input.cell_line_order$AUC][1:30])) == 2 )
stopifnot( length(table(input.labels_cell_lines$IC50[input.cell_line_order$IC50][1:30])) == 2 )

patient_23 <- foreach(input.training_amount.c = seq(from = 30, to = 605, by = 30), .errorhandling = "stop") %dopar% {
  generate_random_partition(labels_cell_lines = input.labels_cell_lines, 
                            labels_patient = docetaxel.labels$patient,
                            training_amount.p = 23,
                            acc_training = TRUE,
                            cell_line_order = input.cell_line_order,
                            training_amount.c = input.training_amount.c)
}
partition$patient_23 <- patient_23

# save worksapce ----------------------------------------------------------
save(docetaxel, docetaxel.labels, sampleinfo.cgp, partition, feature.l1000, partition, 
     file = "Docetaxel/WS/docetaxel_data.RData")
