library("doParallel")
library("sva")

source('Docetaxel/Script/generate_random_partition.R')
source('Common/preparing_data_helper.R')
source('Common/drug_cut/callingWaterfall.R')
source('Common/drug_cut/distancePointLine.R')
source('Common/drug_cut/distancePointSegment.R')
source('Common/comGENE.R')

###

erlotinib.labels <- list()
load("CGP/cdrug2_cgp_ccle_all.RData")
#### Get the new Sensitivity Data
cgp_sensitivity_1 <- read.csv("CGP/cgp_sensitivity_1.csv")

temp.erlotinib_ind <- which(cgp_sensitivity_1$drug.name == "Erlotinib")
cgp_sensitivity_1 <- cgp_sensitivity_1[temp.erlotinib_ind, ]
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

erlotinib.labels$slope <- temp.labels_slope[temp.slope_ind]
erlotinib.labels$slope_ind <- which(rownames(data.ge.cgp) %in% names(temp.labels_slope))
length(erlotinib.labels$slope_ind)

stopifnot(names(erlotinib.labels$slope) == rownames(data.ge.cgp)[erlotinib.labels$slope_ind])

### END
###
temp.drug_ind <- which(druginfo.cgp$drug.name == "ERLOTINIB")

temp.ic50 <- drugpheno.cgp$IC50[, temp.drug_ind]
temp.ic50_ind <- which(!is.na(temp.ic50))
temp.ic50 <- temp.ic50[temp.ic50_ind]
length(temp.ic50)

temp.auc <- drugpheno.cgp$AUC[, temp.drug_ind]
temp.auc_ind <- which(!is.na(temp.auc))
temp.auc <- temp.auc[temp.auc_ind]
length(temp.auc)

source('Common/drug_cut/callingWaterfall.R')
source('Common/drug_cut/distancePointLine.R')
source('Common/drug_cut/distancePointSegment.R')

temp.response <- callingWaterfall(temp.ic50, type="IC50")
temp.response_auc <- callingWaterfall(temp.auc, type="AUC")

temp.label_IC50 <- temp.response != "resistant"
names(temp.label_IC50) <- names(temp.response)

temp.label_auc <- temp.response != "resistant"
names(temp.label_auc) <- names(temp.response_auc)

# check agreement between AUC and IC50
stopifnot(names(temp.label_auc) == names(temp.label_IC50))
sum(temp.label_auc == temp.label_IC50)
## 314 / 323 for AUC and IC50

#### Getting the erlotinib labels from IC50
erlotinib.labels$IC50_ind <- which(rownames(data.ge.cgp) %in% names(temp.label_IC50))
erlotinib.labels$AUC_ind <- which(rownames(data.ge.cgp) %in% names(temp.label_auc))

stopifnot(names(temp.label_IC50) == rownames(data.ge.cgp)[erlotinib.labels$IC50_ind])
erlotinib.labels$IC50 <- temp.label_IC50
stopifnot(names(temp.label_auc) == rownames(data.ge.cgp)[erlotinib.labels$AUC_ind])
erlotinib.labels$AUC <- temp.label_auc

## check how much of the labels agree with each other for IC50 and slope
temp.ind <- names(erlotinib.labels$IC50) %in% names(erlotinib.labels$slope)
stopifnot(sum(temp.ind == FALSE) == length(erlotinib.labels$IC50) - length(erlotinib.labels$slope))
length(temp.ind)
length(erlotinib.labels$IC50)
length(erlotinib.labels$IC50[temp.ind])
length(erlotinib.labels$slope)
temp.slope_include <- which(!is.na(erlotinib.labels$slope))
stopifnot(names(erlotinib.labels$slope[temp.slope_include]) == names(erlotinib.labels$IC50[temp.ind]))
sum(erlotinib.labels$slope[temp.slope_include] == erlotinib.labels$IC50[temp.ind])
# 297 / 309 in with IC50 and slope

## check how much of the labels agree with each other for IC50 and slope
temp.ind <- names(erlotinib.labels$AUC) %in% names(erlotinib.labels$slope)
stopifnot(sum(temp.ind == FALSE) == length(erlotinib.labels$AUC) - length(erlotinib.labels$slope))
length(temp.ind)
length(erlotinib.labels$AUC)
length(erlotinib.labels$AUC[temp.ind])
length(erlotinib.labels$slope)
temp.slope_include <- which(!is.na(erlotinib.labels$slope))

stopifnot(names(erlotinib.labels$slope[temp.slope_include]) == names(erlotinib.labels$AUC[temp.ind]))
sum(erlotinib.labels$slope[temp.slope_include] == erlotinib.labels$AUC[temp.ind])
# 304 / 348 with AUC and slope

###################### 
source('Common/preparing_data_helper.R')
## get the data for slopes
erlotinib <- list()
erlotinib$cgp_slope <- data.ge.cgp[erlotinib.labels$slope_ind,  ]
erlotinib$cgp_IC50 <- data.ge.cgp[erlotinib.labels$IC50_ind,  ]
erlotinib$cgp_AUC <- data.ge.cgp[erlotinib.labels$AUC_ind,  ]

temp <- remove_features(input_data = erlotinib$cgp_slope, input_ind = erlotinib.labels$slope_ind,
                        input_label = erlotinib.labels$slope)
erlotinib$cgp_slope <- temp$data
erlotinib.labels$slope_ind <- temp$ind
erlotinib.labels$slope <- temp$label

rm(temp)
temp <- remove_features(input_data = erlotinib$cgp_IC50, input_ind = erlotinib.labels$IC50_ind,
                        input_label = erlotinib.labels$IC50)
erlotinib$cgp_IC50 <- temp$data
erlotinib.labels$IC50_ind <- temp$ind
erlotinib.labels$IC50 <- temp$label

rm(temp)
temp <- remove_features(input_data = erlotinib$cgp_AUC, input_ind = erlotinib.labels$AUC_ind,
                        input_label = erlotinib.labels$AUC)
erlotinib$cgp_AUC <- temp$data
erlotinib.labels$AUC_ind <- temp$ind
erlotinib.labels$AUC <- temp$label

## check agreement
temp.ind <- match(names(erlotinib.labels$slope), names(erlotinib.labels$AUC))
stopifnot(names(erlotinib.labels$AUC)[temp.ind] == names(erlotinib.labels$slope))
sum(erlotinib.labels$AUC[temp.ind] == erlotinib.labels$slope)
sum(!is.na(temp.ind))
# 218 / 273 for slope vs. AUC
temp.ind <- match(names(erlotinib.labels$slope), names(erlotinib.labels$IC50))
sum(erlotinib.labels$IC50[temp.ind] == erlotinib.labels$slope)
sum(!is.na(temp.ind))
# 261  /  273 for slope vs. IC50
temp.ind <- match(names(erlotinib.labels$AUC), names(erlotinib.labels$IC50))
sum(erlotinib.labels$IC50[temp.ind] == erlotinib.labels$AUC)
sum(!is.na(temp.ind))
# 278 /  287 for AUC vs. IC50 


# Using sva to harmonize across different tissue types
source('Common/comGENE.R')
#load("erlotinib_data.RData")

show_pca(input_data = erlotinib$cgp_slope, label = erlotinib.labels$slope)
show_pca(input_data = erlotinib$cgp_AUC, label = erlotinib.labels$AUC)
show_pca(input_data = erlotinib$cgp_IC50, label = erlotinib.labels$IC50)

##
# get patient data 
load("Erlotinib/WS/clinical_data.RData")
temp.patient <- erlotinib$patient
temp.label <- erlotinib$patient_label
names(temp.label) <- rownames(erlotinib$patient)

## get the GEO data as well
temp.geo_cell_line <- erlotinib$cell_line.geo
temp.geo_cell_line.label <- erlotinib$cell_line.geo.label == 1
names(temp.geo_cell_line.label) <- rownames(erlotinib$cell_line.geo)
#load("Erlotinib/WS/erlotinib_data.RData")
erlotinib$GEO_IC50 <- temp.geo_cell_line
erlotinib.labels$GEO_IC50 <- temp.geo_cell_line.label
erlotinib.labels$GEO_IC50_combined <- c(erlotinib.labels$GEO_IC50, erlotinib.labels$patient)
erlotinib.labels$GEO_IC50_combined.source <- c(rep("GEO", length(erlotinib.labels$GEO_IC50)), rep("patient", length(erlotinib.labels$patient)))
## combine with SVA
## doing the GEO
stopifnot(colnames(erlotinib$GEO_IC50) == colnames(erlotinib$patient))  
mean(erlotinib$GEO_IC50)
mean(erlotinib$patient)
dim(erlotinib$GEO_IC50)
dim(erlotinib$patient)

erlotinib$combined_GEO <- rbind(erlotinib$GEO_IC50, erlotinib$patient)
show_pca(input_data = erlotinib$combined_GEO, label = erlotinib.labels$GEO_IC50_combined)
show_pca(input_data = rbind(erlotinib$GEO_IC50, erlotinib$patient), label = erlotinib.labels$GEO_IC50_combined.source)

## GEO and CGP lung only 
temp.cgp_lung_only.IC50 <- which(sampleinfo.cgp$tissue.type[erlotinib.labels$IC50_ind] == "lung")
temp.data <- comGENE(erlotinib$GEO_IC50, scale(erlotinib$cgp_IC50[temp.cgp_lung_only.IC50, ]))
mean(temp.data[[1]])
mean(temp.data[[2]])

erlotinib.labels$lung_all.IC50 <- c(erlotinib.labels$GEO_IC50, erlotinib.labels$IC50[temp.cgp_lung_only.IC50], erlotinib.labels$patient)
erlotinib.labels$lung_all.IC50.source <- c(rep("GEO", dim(temp.data[[1]])[1]), rep("CGP", dim(temp.data[[2]])[1]), rep("patient", length(erlotinib.labels$patient)))

stopifnot(names(erlotinib.labels$lung_all.IC50) == rownames(rbind(temp.data[[1]], temp.data[[2]], erlotinib$patient)))

erlotinib$lung_all.IC50.ComBat <- ComBat_combine(input_data = rbind(temp.data[[1]], temp.data[[2]], erlotinib$patient), 
                                                 label = erlotinib.labels$lung_all.IC50, 
                                                 batch = erlotinib.labels$lung_all.IC50.source)

dim(erlotinib$lung_all.IC50)
length(erlotinib.labels$lung_all.IC50.source)
length(erlotinib.labels$lung_all.IC50)

show_pca(input_data = erlotinib$lung_all.IC50.ComBat, label = erlotinib.labels$lung_all.IC50)
show_pca(input_data = erlotinib$lung_all.IC50.ComBat, label = erlotinib.labels$lung_all.IC50.source)

# slope
temp.cgp_lung_only.slope <- which(sampleinfo.cgp$tissue.type[erlotinib.labels$slope_ind] == "lung")
temp.data <- comGENE(erlotinib$GEO_IC50, scale(erlotinib$cgp_slope[temp.cgp_lung_only.slope, ]))
mean(temp.data[[1]])
mean(temp.data[[2]])

erlotinib.labels$lung_all.slope <- c(erlotinib.labels$GEO_IC50, erlotinib.labels$slope[temp.cgp_lung_only.slope], erlotinib.labels$patient)
erlotinib.labels$lung_all.slope.source <- c(rep("GEO", dim(temp.data[[1]])[1]), rep("CGP", dim(temp.data[[2]])[1]), rep("patient", length(erlotinib.labels$patient)))

stopifnot(names(erlotinib.labels$lung_all.slope) == rownames(rbind(temp.data[[1]], temp.data[[2]], erlotinib$patient)))


erlotinib$lung_all.slope.ComBat <- ComBat_combine(input_data = rbind(temp.data[[1]], temp.data[[2]], erlotinib$patient), 
                                                  label = erlotinib.labels$lung_all.slope, 
                                                  batch = erlotinib.labels$lung_all.slope.source)

dim(erlotinib$lung_all.slope.ComBat)
length(erlotinib.labels$lung_all.slope.source)
length(erlotinib.labels$lung_all.slope)

show_pca(input_data = erlotinib$lung_all.slope.ComBat, label = erlotinib.labels$lung_all.slope)
show_pca(input_data = erlotinib$lung_all.slope.ComBat, label = erlotinib.labels$lung_all.slope.source)

## all of CGP + GEO
# IC50
temp.data <- comGENE(erlotinib$GEO_IC50, scale(erlotinib$cgp_IC50))
mean(temp.data[[1]])
mean(temp.data[[2]])
erlotinib.labels$all.source.IC50 <- c(rep("GEO", dim(temp.data[[1]])[1]), rep("CGP", dim(temp.data[[2]])[1]), rep("patient", length(erlotinib.labels$patient)))
erlotinib.labels$all.IC50 <- c(erlotinib.labels$GEO_IC50, erlotinib.labels$IC50, erlotinib.labels$patient)

erlotinib$all.ComBat.IC50 <- ComBat_combine(batch = erlotinib.labels$all.source.IC50,
                                            label = erlotinib.labels$all.IC50, 
                                            input_data = rbind(temp.data[[1]], temp.data[[2]], erlotinib$patient), par.prior = FALSE)
show_pca(input_data = erlotinib$all.ComBat.IC50, label = erlotinib.labels$all.source.IC50)

# slope
temp.data <- comGENE(erlotinib$GEO_IC50, scale(erlotinib$cgp_slope))
mean(temp.data[[1]])
mean(temp.data[[2]])
erlotinib.labels$all.source.slope <- c(rep("GEO", dim(temp.data[[1]])[1]), rep("CGP", dim(temp.data[[2]])[1]), rep("patient", length(erlotinib.labels$patient)))
erlotinib.labels$all.slope <- c(erlotinib.labels$GEO_IC50, erlotinib.labels$slope, erlotinib.labels$patient)

erlotinib$all.ComBat.slope <- ComBat_combine(batch = erlotinib.labels$all.source.slope,
                                             label = erlotinib.labels$all.slope, 
                                             input_data = rbind(temp.data[[1]], temp.data[[2]], erlotinib$patient))
show_pca(input_data = erlotinib$all.ComBat.slope, label = erlotinib.labels$all.source.slope)


### get the partitioning
temp.test_amount <- list(cc = round(0.2 * length(erlotinib.labels$slope)))

partition <- list()

# GEO
temp.test_amount <- list(cc = 11)
partition$GEO <- generate_random_partition(input_labels_cell_lines = erlotinib.labels$GEO_IC50, 
                                           input_labels_patient = erlotinib.labels$patient, 
                                           test_amount = temp.test_amount)

# lung only IC50
temp.test_amount <- list(cc = round(0.2 * length(which(erlotinib.labels$lung_all.IC50.source != "patient"))))
partition$lung.IC50 <- generate_random_partition(input_labels_cell_lines = erlotinib.labels$lung_all.IC50[which(erlotinib.labels$lung_all.IC50.source != "patient")], 
                                                 input_labels_patient = erlotinib.labels$patient, 
                                                 test_amount = temp.test_amount)

# lung only slope
temp.test_amount <- list(cc = round(0.2 * length(which(erlotinib.labels$lung_all.slope.source != "patient"))))
partition$lung.slope <- generate_random_partition(input_labels_cell_lines = erlotinib.labels$lung_all.slope[which(erlotinib.labels$lung_all.slope.source != "patient")], 
                                                  input_labels_patient = erlotinib.labels$patient, 
                                                  test_amount = temp.test_amount)

# all IC50
temp.test_amount <- list(cc = round(0.2 * length(which(erlotinib.labels$all.source.IC50 != "patient"))))
partition$all.IC50 <- generate_random_partition(input_labels_cell_lines = erlotinib.labels$all.IC50[which(erlotinib.labels$all.source.IC50 != "patient")], 
                                                input_labels_patient = erlotinib.labels$patient, 
                                                test_amount = temp.test_amount)

# all slope
temp.test_amount <- list(cc = round(0.2 * length(which(erlotinib.labels$all.source.slope != "patient"))))
partition$all.slope <- generate_random_partition(input_labels_cell_lines = erlotinib.labels$all.slope[which(erlotinib.labels$all.source.slope != "patient")], 
                                                 input_labels_patient = erlotinib.labels$patient, 
                                                 test_amount = temp.test_amount)

### get l1000 features
stopifnot(colnames(erlotinib$cgp_IC50) == colnames(erlotinib$cgp_slope))

Landmark_Genes_n978 <- read.csv("Common/Landmark_Genes_n978.csv")

feature.l1000 <- list()
feature.l1000$cc <- which(colnames(erlotinib$cgp_slope) %in% Landmark_Genes_n978$Gene.Symbol)
feature.l1000$cp <- which(colnames(erlotinib$all.ComBat.slope) %in% Landmark_Genes_n978$Gene.Symbol)
length(feature.l1000$cc)
length(feature.l1000$cp)

feature.l1000$geo_cc <- which(colnames(erlotinib$GEO_IC50) %in% Landmark_Genes_n978$Gene.Symbol)
#feature.l1000$geo_cp <- which(colnames(erlotinib$) %in% Landmark_Genes_n978$Gene.Symbol)
length(feature.l1000$geo_cc)
length(feature.l1000$geo_cp)

feature.l1000$pp <- which(colnames(erlotinib$patient) %in% Landmark_Genes_n978$Gene.Symbol)
length(feature.l1000$pp)

getwd()
setwd("Erlotinib/WS/")
#save(erlotinib, erlotinib.labels, sampleinfo.cgp, partition, feature.l1000, file = "erlotinib_data.RData")


### create partitions for varying number of patients
partition_var <- list()

temp.cp <- foreach (training_amount.cp = seq(from = 30, to = 300, by = 30)) %do% {    
  
  temp.slope <- generate_random_partition.cp_var(input_labels_cell_lines = c(erlotinib.labels$GEO_IC50, erlotinib.labels$slope), 
                                                 input_labels_patient = erlotinib.labels$patient, 
                                                 training_amount = training_amount.cp)
  
  temp.IC50 <- generate_random_partition.cp_var(input_labels_cell_lines = c(erlotinib.labels$GEO_IC50, erlotinib.labels$IC50), 
                                                input_labels_patient = erlotinib.labels$patient, 
                                                training_amount = training_amount.cp)  
  
  list(slope = temp.slope, IC50 = temp.IC50)
}

# double check
for (temp.ind in 1:length(seq(from = 30, to = 300, by = 30))) {
  stopifnot(length(temp.cp[[temp.ind]]$slope[[100]]$test_index) == length(erlotinib.labels$patient))
  stopifnot(length(temp.cp[[temp.ind]]$slope[[100]]$training_index.single) == temp.ind * 30)
  
  temp.maxIC50 = 0
  temp.maxSlope = 0
  for (temp.ind2 in 1:100) {
    stopifnot(min(table(erlotinib.labels$all.slope[temp.cp[[temp.ind]]$slope[[temp.ind2]]$training_index.single])) >= 5)
    stopifnot(min(table(erlotinib.labels$all.IC50[temp.cp[[temp.ind]]$IC50[[temp.ind2]]$training_index.single])) >= 5)
    temp.maxIC50 = max(temp.max = 0, temp.cp[[temp.ind]]$IC50[[temp.ind2]]$training_index.single)
    temp.maxSlope = max(temp.max = 0, temp.cp[[temp.ind]]$slope[[temp.ind2]]$training_index.single)
  }
  print(temp.maxIC50)
  print(temp.maxSlope)  
}

partition_var$cp <- temp.cp

##################### 

temp.cp <- foreach (training_amount.cp = seq(from = 10, to = 80, by = 10)) %do% {    
  
  temp.slope <- generate_random_partition.cp_var(input_labels_cell_lines = erlotinib.labels$lung_all.slope[erlotinib.labels$lung_all.slope.source != "patient"], 
                                                 input_labels_patient = erlotinib.labels$patient, 
                                                 training_amount = training_amount.cp)
  
  temp.IC50 <- generate_random_partition.cp_var(input_labels_cell_lines = erlotinib.labels$lung_all.IC50[erlotinib.labels$lung_all.IC50.source != "patient"], 
                                                input_labels_patient = erlotinib.labels$patient, 
                                                training_amount = training_amount.cp)  
  
  list(slope = temp.slope, IC50 = temp.IC50)
}

# double check
for (temp.ind in 1:length(seq(from = 10, to = 80, by = 10))) {
  stopifnot(length(temp.cp[[temp.ind]]$slope[[100]]$test_index) == length(erlotinib.labels$patient))
  stopifnot(length(temp.cp[[temp.ind]]$slope[[100]]$training_index.single) == temp.ind * 10)
}

partition_var$cp_lung_only <- temp.cp

#####################

temp.cpp <- foreach (training_amount.cpp = seq(from = 1, to = 23, by = 1)) %do% {    
  
  temp.slope <- generate_random_partition.cpp_var(input_labels_cell_lines = c(erlotinib.labels$GEO_IC50, erlotinib.labels$slope), 
                                                  input_labels_patient = erlotinib.labels$patient, 
                                                  training_amount = training_amount.cpp)
  
  temp.IC50 <- generate_random_partition.cpp_var(input_labels_cell_lines = c(erlotinib.labels$GEO_IC50, erlotinib.labels$IC50), 
                                                 input_labels_patient = erlotinib.labels$patient, 
                                                 training_amount = training_amount.cpp)  
  
  list(slope = temp.slope, IC50 = temp.IC50)
}

# double check
for (temp.ind in 1:length(seq(from = 1, to = 23, by = 1))) {
  stopifnot(length(temp.cpp[[temp.ind]]$slope[[100]]$test_index) == length(erlotinib.labels$patient) - seq(from = 1, to = 23, by = 1)[temp.ind])
  stopifnot(length(temp.cpp[[temp.ind]]$slope[[100]]$training_index.single) == -length(temp.cpp[[temp.ind]]$slope[[100]]$test_index) + length(erlotinib.labels$all.slope))
  
  temp.maxIC50 = 0
  temp.maxSlope = 0
  for (temp.ind2 in 1:100) {
    stopifnot(min(table(erlotinib.labels$all.slope[temp.cpp[[temp.ind]]$slope[[temp.ind2]]$training_index.single])) >= 5)
    stopifnot(min(table(erlotinib.labels$all.slope[temp.cpp[[temp.ind]]$slope[[temp.ind2]]$training_index.single])) >= 5)
    temp.maxIC50 = max(temp.max = 0, temp.cpp[[temp.ind]]$IC50[[temp.ind2]]$training_index.single)
    temp.maxSlope = max(temp.max = 0, temp.cpp[[temp.ind]]$slope[[temp.ind2]]$training_index.single)
  }
  print(temp.maxIC50)
  print(temp.maxSlope) 
}

partition_var$cpp <- temp.cpp

########################

temp.cpp <- foreach (training_amount.cpp = seq(from = 1, to = 23, by = 1)) %do% {    
  
  temp.slope <- generate_random_partition.cpp_var(input_labels_cell_lines = erlotinib.labels$lung_all.slope[erlotinib.labels$lung_all.slope.source != "patient"], 
                                                  input_labels_patient = erlotinib.labels$patient, 
                                                  training_amount = training_amount.cpp)
  
  temp.IC50 <- generate_random_partition.cpp_var(input_labels_cell_lines = erlotinib.labels$lung_all.IC50[erlotinib.labels$lung_all.IC50.source != "patient"], 
                                                 input_labels_patient = erlotinib.labels$patient, 
                                                 training_amount = training_amount.cpp)  
  
  list(slope = temp.slope, IC50 = temp.IC50)
}


# double check
for (temp.ind in 1:length(seq(from = 1, to = 23, by = 1))) {
  stopifnot(length(temp.cpp[[temp.ind]]$slope[[100]]$test_index) == length(erlotinib.labels$patient) - seq(from = 1, to = 23, by = 1)[temp.ind])
  #stopifnot(length(temp.cpp[[temp.ind]]$slope[[100]]$training_index.single) == -length(temp.cpp[[temp.ind]]$slope[[100]]$test_index) + length(erlotinib.labels$lung_all.slope))
}

partition_var$cpp_lung_only <- temp.cpp


########################

## patients
temp.pp <- foreach (training_amount.cpp = seq(from = 14, to = 23, by = 1)) %do% {      
  generate_random_partition.pp_var(input_labels_patient = erlotinib.labels$patient, 
                                   training_amount = training_amount.cpp)
}

# double check
for (temp.ind in 14:23) {
  stopifnot(length(temp.pp[[temp.ind - 13]][[100]]$test_index) == length(erlotinib.labels$patient) - temp.ind)
  stopifnot(length(temp.pp[[temp.ind - 13]][[100]]$training_index.single) == temp.ind)
}

partition_var$pp <- temp.pp

#######################

temp.cpp <- foreach (cell_line_training_amount = seq(from = 20, to = 310, by = 20)) %do% {  
  
  temp.slope <- generate_random_partition.cpp_var2(input_labels_cell_lines = c(erlotinib.labels$GEO_IC50, erlotinib.labels$slope), 
                                                   input_labels_patient = erlotinib.labels$patient, 
                                                   cell_line_training_amount = cell_line_training_amount,
                                                   patient_training = 20)
  
  temp.IC50 <- generate_random_partition.cpp_var2(input_labels_cell_lines = c(erlotinib.labels$GEO_IC50, erlotinib.labels$IC50), 
                                                  input_labels_patient = erlotinib.labels$patient, 
                                                  cell_line_training_amount = cell_line_training_amount,
                                                  patient_training = 20)
  
  
  list(slope = temp.slope, IC50 = temp.IC50)
}

partition_var$cVar_p20_p <- temp.cpp



setwd("Erlotinib/WS")
#save(erlotinib, erlotinib.labels, sampleinfo.cgp, partition, feature.l1000, partition_var, file = "erlotinib_data.RData")
