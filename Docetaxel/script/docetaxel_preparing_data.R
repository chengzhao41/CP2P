stop("Set working directory to current source file")
#setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Code/Docetaxel/Script")

docetaxel.labels <- list()
load("../../CGP/cdrug2_cgp_ccle_all.RData")
#### Get the new Sensitivity Data
cgp_sensitivity_1 <- read.csv("../../CGP/cgp_sensitivity_1.csv")

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

### END
###
temp.drug_ind <- which(druginfo.cgp$drug.name == "DOCETAXEL")

temp.ic50 <- drugpheno.cgp$IC50[, temp.drug_ind]
temp.ic50_ind <- which(!is.na(temp.ic50))
temp.ic50 <- temp.ic50[temp.ic50_ind]
length(temp.ic50)

temp.auc <- drugpheno.cgp$AUC[, temp.drug_ind]
temp.auc_ind <- which(!is.na(temp.auc))
temp.auc <- temp.auc[temp.auc_ind]
length(temp.auc)

source('../../Common/drug_cut/callingWaterfall.R')
source('../../Common/drug_cut/distancePointLine.R')
source('../../Common/drug_cut/distancePointSegment.R')

temp.response <- callingWaterfall(temp.ic50, type="IC50")
temp.response_auc <- callingWaterfall(temp.auc, type="AUC")
table(temp.response)

temp.label_IC50 <- temp.response != "resistant"
names(temp.label_IC50) <- names(temp.response)

temp.label_auc <- temp.response_auc != "resistant"
names(temp.label_auc) <- names(temp.response_auc)

# check agreement between AUC and IC50
stopifnot(names(temp.label_auc) == names(temp.label_IC50))
sum(temp.label_auc == temp.label_IC50)
## 548 / 663 for AUC and IC50

#### Getting the docetaxel labels from IC50
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
sum(docetaxel.labels$slope[temp.slope_include] == docetaxel.labels$IC50[temp.ind])
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
sum(docetaxel.labels$slope[temp.slope_include] == docetaxel.labels$AUC[temp.ind])
# 567 / 650 with AUC and slope

###################### 
source('../../Common/preparing_data_helper.R')

## get the data for slopes
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
sum(docetaxel.labels$AUC[temp.ind] == docetaxel.labels$slope)
sum(!is.na(temp.ind))
# 528 / 605 for slope vs. AUC
temp.ind <- match(names(docetaxel.labels$slope), names(docetaxel.labels$IC50))
sum(docetaxel.labels$IC50[temp.ind] == docetaxel.labels$slope)
sum(!is.na(temp.ind))
# 427 / 605 for slope vs. IC50
temp.ind <- match(names(docetaxel.labels$AUC), names(docetaxel.labels$IC50))
sum(docetaxel.labels$IC50[temp.ind] == docetaxel.labels$AUC)
sum(!is.na(temp.ind))
# 511 / 618 for AUC vs. IC50 

#getwd()
#save(docetaxel, docetaxel.labels, sampleinfo.cgp, file = "../WS/docetaxel_data.RData")

# Using sva to harmonize across different tissue types
library("sva")
source('../../Common/comGENE.R')

# get patient data 
load("../WS/pp.RData")
docetaxel$patient <- docetaxel.patient
docetaxel.labels$patient <- pp.ground_truth == 1

show_pca(input_data = docetaxel$cgp_slope, label = docetaxel.labels$slope)
show_pca(input_data = docetaxel$cgp_AUC, label = docetaxel.labels$AUC)
show_pca(input_data = docetaxel$cgp_IC50, label = docetaxel.labels$IC50)

## Slopes
# Using sva to harmonize patients and cell lines
docetaxel.labels$slope_combined <- c(docetaxel.labels$slope, docetaxel.labels$patient)
stopifnot(substring(colnames(docetaxel$cgp_slope), 8) == annot.ge.cgp$EntrezGene.ID)
colnames(docetaxel$cgp_slope) <- annot.ge.cgp$symbol

temp.data <- comGENE(scale(docetaxel$cgp_slope), docetaxel$patient)
mean(temp.data[[1]])
mean(temp.data[[2]])

docetaxel.labels$slope_combined.source <- c(rep("cgp", dim(temp.data[[1]])[1]), rep("patient", dim(temp.data[[2]])[1]))

docetaxel$combined_slope.ComBat <- ComBat_combine(batch = docetaxel.labels$slope_combined.source,
                                                  label = docetaxel.labels$slope_combined,
                                                  input_data = rbind(temp.data[[1]], temp.data[[2]]))
show_pca(input_data = docetaxel$combined_slope.ComBat, label = docetaxel.labels$slope_combined)
show_pca(input_data = docetaxel$combined_slope.ComBat, label = c(sampleinfo.cgp$tissue.type[docetaxel.labels$slope_ind], rep("patient", dim(temp.data[[2]])[1])) == "patient")         
rm(temp.data)

## IC50
# Using sva to harmonize patients and cell lines
docetaxel.labels$IC50_combined <- c(docetaxel.labels$IC50, docetaxel.labels$patient)
stopifnot(substring(colnames(docetaxel.labels$IC50), 8) == annot.ge.cgp$EntrezGene.ID)
colnames(docetaxel$cgp_IC50) <- annot.ge.cgp$symbol

source('~/Dropbox/SNF_DRUG_PROJECT/Script/comGENE.R')
temp.data <- comGENE(scale(docetaxel$cgp_IC50), docetaxel$patient)
mean(temp.data[[1]])
mean(temp.data[[2]])

docetaxel.labels$IC50_combined.source <- c(rep("cgp", dim(temp.data[[1]])[1]), rep("patient", dim(temp.data[[2]])[1]))
docetaxel$combined_IC50 <- rbind(temp.data[[1]], temp.data[[2]])
show_pca(input_data = docetaxel$combined_IC50, label = c(rep("cgp", dim(temp.data[[1]])[1]), rep("patient", dim(temp.data[[2]])[1])))

docetaxel$combined_IC50.ComBat <- ComBat_combine(batch = docetaxel.labels$IC50_combined.source,
                                                 label = docetaxel.labels$IC50_combined,
                                                 input_data = docetaxel$combined_IC50)
show_pca(input_data = docetaxel$combined_IC50.ComBat, label = docetaxel.labels$IC50_combined)
show_pca(input_data = docetaxel$combined_IC50.ComBat, label = docetaxel.labels$IC50_combined.source)

rm(temp.data)
getwd()

## AUC
# Using sva to harmonize patients and cell lines
docetaxel.labels$AUC_combined <- c(docetaxel.labels$AUC, docetaxel.labels$patient)
stopifnot(substring(colnames(docetaxel.labels$AUC), 8) == annot.ge.cgp$EntrezGene.ID)
colnames(docetaxel$cgp_AUC) <- annot.ge.cgp$symbol

source('~/Dropbox/SNF_DRUG_PROJECT/Script/comGENE.R')
temp.data <- comGENE(docetaxel$cgp_AUC, docetaxel$patient)
mean(temp.data[[1]])
mean(temp.data[[2]])

docetaxel$combined_AUC.ComBat <- ComBat_combine(batch = docetaxel.labels$AUC_combined.source,
                                                label = docetaxel.labels$AUC_combined, 
                                                input_data = rbind(temp.data[[1]], temp.data[[2]]))
show_pca(input_data = docetaxel$combined_AUC.ComBat, label = docetaxel.labels$AUC_combined)
show_pca(input_data = docetaxel$combined_AUC.ComBat, label = docetaxel.labels$AUC_combined.source == "patient")

rm(temp.data)

### Breast Cell lines only
## IC50
temp.ind <- which(sampleinfo.cgp$tissue.type[docetaxel.labels$IC50_ind] == "breast")
temp.data <- comGENE(scale(docetaxel$cgp_IC50[temp.ind, ]), docetaxel$patient)
mean(temp.data[[1]])
mean(temp.data[[2]])

docetaxel.labels$IC50_breast.source <- c(rep("cgp", dim(temp.data[[1]])[1]), rep("patient", dim(temp.data[[2]])[1]))
docetaxel.labels$IC50_breast_only <- docetaxel.labels$IC50[temp.ind]
docetaxel.labels$IC50_breast <- c(docetaxel.labels$IC50[temp.ind], docetaxel.labels$patient)

docetaxel$breast_IC50.ComBat <- ComBat_combine(batch = docetaxel.labels$IC50_breast.source,
                                               label = docetaxel.labels$IC50_breast,
                                               input_data = rbind(temp.data[[1]], temp.data[[2]]))
show_pca(input_data = docetaxel$breast_IC50.ComBat, label = docetaxel.labels$IC50_breast)
show_pca(input_data = docetaxel$breast_IC50.ComBat, label = docetaxel.labels$IC50_breast.source)
rm(temp.data)
getwd()

# AUC
temp.ind <- which(sampleinfo.cgp$tissue.type[docetaxel.labels$AUC_ind] == "breast")
temp.data <- comGENE(scale(docetaxel$cgp_AUC[temp.ind, ]), docetaxel$patient)
mean(temp.data[[1]])
mean(temp.data[[2]])

docetaxel.labels$AUC_breast.source <- c(rep("cgp", dim(temp.data[[1]])[1]), rep("patient", dim(temp.data[[2]])[1]))
docetaxel.labels$AUC_breast_only <- docetaxel.labels$AUC[temp.ind]
docetaxel.labels$AUC_breast <- c(docetaxel.labels$AUC[temp.ind], docetaxel.labels$patient)

docetaxel$breast_AUC.ComBat <- ComBat_combine(batch = docetaxel.labels$AUC_breast.source,
                                              label = docetaxel.labels$AUC_breast,
                                              input_data = rbind(temp.data[[1]], temp.data[[2]]))
show_pca(input_data = docetaxel$breast_AUC.ComBat, label = docetaxel.labels$AUC_breast)
show_pca(input_data = docetaxel$breast_AUC.ComBat, label = docetaxel.labels$AUC_breast.source)
rm(temp.data)
getwd()

# Slope
temp.ind <- which(sampleinfo.cgp$tissue.type[docetaxel.labels$slope_ind] == "breast")
temp.data <- comGENE(scale(docetaxel$cgp_slope[temp.ind, ]), docetaxel$patient)
mean(temp.data[[1]])
mean(temp.data[[2]])

docetaxel.labels$slope_breast.source <- c(rep("cgp", dim(temp.data[[1]])[1]), rep("patient", dim(temp.data[[2]])[1]))
docetaxel.labels$slope_breast_only <- docetaxel.labels$slope[temp.ind]
docetaxel.labels$slope_breast <- c(docetaxel.labels$slope[temp.ind], docetaxel.labels$patient)

docetaxel$breast_slope.ComBat <- ComBat_combine(batch = docetaxel.labels$slope_breast.source,
                                                label = docetaxel.labels$slope_breast,
                                                input_data = rbind(temp.data[[1]], temp.data[[2]]))
show_pca(input_data = docetaxel$breast_slope.ComBat, label = docetaxel.labels$slope_breast)
show_pca(input_data = docetaxel$breast_slope.ComBat, label = docetaxel.labels$slope_breast.source)
rm(temp.data)


### get the partitioning
source('../Script/generate_random_partition.R')

partition <- list()

partition$slope <- generate_random_partition(input_labels_cell_lines = docetaxel.labels$slope, 
                                             input_labels_patient = docetaxel.labels$patient, 
                                             test_amount = temp.test_amount)

temp.test_amount <- list(cc = round(0.2 * length(docetaxel.labels$IC50)))
partition$IC50 <- generate_random_partition(input_labels_cell_lines = docetaxel.labels$IC50, 
                                            input_labels_patient = docetaxel.labels$patient, 
                                            test_amount = temp.test_amount)

temp.test_amount <- list(cc = round(0.2 * length(docetaxel.labels$AUC)))
partition$AUC <- generate_random_partition(input_labels_cell_lines = docetaxel.labels$AUC, 
                                           input_labels_patient = docetaxel.labels$patient, 
                                           test_amount = temp.test_amount)

temp.test_amount <- list(cc = 10)
partition$slope_breast <- generate_random_partition(input_labels_cell_lines = docetaxel.labels$slope_breast_only, 
                                                    input_labels_patient = docetaxel.labels$patient, 
                                                    test_amount = temp.test_amount)

partition$IC50_breast <- generate_random_partition(input_labels_cell_lines = docetaxel.labels$IC50_breast_only, 
                                                   input_labels_patient = docetaxel.labels$patient, 
                                                   test_amount = temp.test_amount)

partition$AUC_breast <- generate_random_partition(input_labels_cell_lines = docetaxel.labels$AUC_breast_only, 
                                                  input_labels_patient = docetaxel.labels$patient, 
                                                  test_amount = temp.test_amount)

### get l1000 features
Landmark_Genes_n978 <- read.csv("../../Common/Landmark_Genes_n978.csv")
stopifnot(colnames(docetaxel$cgp_IC50) == colnames(docetaxel$cgp_AUC))
stopifnot(colnames(docetaxel$cgp_slope) == colnames(docetaxel$cgp_AUC))
feature.l1000 <- list()
feature.l1000$cc <- which(colnames(docetaxel$cgp_slope) %in% Landmark_Genes_n978$Gene.Symbol)
feature.l1000$cp <- which(colnames(docetaxel$combined_slope.ComBat) %in% Landmark_Genes_n978$Gene.Symbol)
feature.l1000$pp <- which(colnames(docetaxel$patient) %in% Landmark_Genes_n978$Gene.Symbol)
length(feature.l1000$cc)
length(feature.l1000$cp)
length(feature.l1000$pp)
getwd()

####
partition_var <- list()

temp.cp <- foreach (training_amount.cp = seq(from = 10, to = 600, by = 30)) %do% {    
  
  temp.slope <- generate_random_partition.cp_var(input_labels_cell_lines = docetaxel.labels$slope, 
                                                 input_labels_patient = docetaxel.labels$patient, 
                                                 training_amount = training_amount.cp)
  
  temp.IC50 <- generate_random_partition.cp_var(input_labels_cell_lines = docetaxel.labels$IC50, 
                                                input_labels_patient = docetaxel.labels$patient, 
                                                training_amount = training_amount.cp)
  
  temp.AUC <- generate_random_partition.cp_var(input_labels_cell_lines = docetaxel.labels$AUC, 
                                               input_labels_patient = docetaxel.labels$patient, 
                                               training_amount = training_amount.cp)
  list(slope = temp.slope, IC50 = temp.IC50, AUC = temp.AUC)
}

# double check
for (temp.ind in 1:20) {
  stopifnot(length(temp.cp[[temp.ind]]$slope[[100]]$test_index) == length(docetaxel.labels$patient))
  stopifnot(length(temp.cp[[temp.ind]]$slope[[100]]$training_index.single) == 10 + (temp.ind - 1) * 30)
}

partition_var$cp <- temp.cp

############ breast only cp
temp.cp <- foreach (training_amount.cp = seq(from = 10, to = 30, by = 5)) %do% {    
  
  temp.slope <- generate_random_partition.cp_var(input_labels_cell_lines = docetaxel.labels$slope_breast_only, 
                                                 input_labels_patient = docetaxel.labels$patient, 
                                                 training_amount = training_amount.cp)
  
  temp.IC50 <- generate_random_partition.cp_var(input_labels_cell_lines = docetaxel.labels$IC50_breast_only, 
                                                input_labels_patient = docetaxel.labels$patient, 
                                                training_amount = training_amount.cp)
  
  temp.AUC <- generate_random_partition.cp_var(input_labels_cell_lines = docetaxel.labels$AUC_breast_only, 
                                               input_labels_patient = docetaxel.labels$patient, 
                                               training_amount = training_amount.cp)
  list(slope = temp.slope, IC50 = temp.IC50, AUC = temp.AUC)
}

# double check
for (temp.ind in 1:5) {
  stopifnot(length(temp.cp[[temp.ind]]$slope[[100]]$test_index) == length(docetaxel.labels$patient))
  stopifnot(length(temp.cp[[temp.ind]]$slope[[100]]$training_index.single) == 10 + (temp.ind - 1) * 5)
}

partition_var$cp_breast_only <- temp.cp
###################
temp.cpp <- foreach (training_amount.cpp = seq(from = 1, to = 22, by = 1)) %do% {    
  
  temp.slope <- generate_random_partition.cpp_var(input_labels_cell_lines = docetaxel.labels$slope, 
                                                  input_labels_patient = docetaxel.labels$patient, 
                                                  training_amount = training_amount.cpp)
  
  temp.IC50 <- generate_random_partition.cpp_var(input_labels_cell_lines = docetaxel.labels$IC50, 
                                                 input_labels_patient = docetaxel.labels$patient, 
                                                 training_amount = training_amount.cpp)
  
  temp.AUC <- generate_random_partition.cpp_var(input_labels_cell_lines = docetaxel.labels$AUC, 
                                                input_labels_patient = docetaxel.labels$patient, 
                                                training_amount = training_amount.cpp)
  list(slope = temp.slope, IC50 = temp.IC50, AUC = temp.AUC)
}

# double check
for (temp.ind in 1:15) {
  stopifnot(length(temp.cpp[[temp.ind]]$slope[[100]]$test_index) == length(docetaxel.labels$patient) - temp.ind)
  stopifnot(length(temp.cpp[[temp.ind]]$slope[[100]]$training_index.single) == temp.ind + length(docetaxel.labels$slope))
}

partition_var$cpp <- temp.cpp
####################

temp.cpp <- foreach (training_amount.cpp = seq(from = 1, to = 22, by = 1)) %do% {    
  
  temp.slope <- generate_random_partition.cpp_var(input_labels_cell_lines = docetaxel.labels$slope_breast_only, 
                                                  input_labels_patient = docetaxel.labels$patient, 
                                                  training_amount = training_amount.cpp)
  
  temp.IC50 <- generate_random_partition.cpp_var(input_labels_cell_lines = docetaxel.labels$IC50_breast_only, 
                                                 input_labels_patient = docetaxel.labels$patient, 
                                                 training_amount = training_amount.cpp)
  
  temp.AUC <- generate_random_partition.cpp_var(input_labels_cell_lines = docetaxel.labels$AUC_breast_only, 
                                                input_labels_patient = docetaxel.labels$patient, 
                                                training_amount = training_amount.cpp)
  list(slope = temp.slope, IC50 = temp.IC50, AUC = temp.AUC)
}

# double check
for (temp.ind in 1:22) {
  stopifnot(length(temp.cpp[[temp.ind]]$slope[[100]]$test_index) == length(docetaxel.labels$patient) - temp.ind)
  stopifnot(length(temp.cpp[[temp.ind]]$slope[[100]]$training_index.single) == temp.ind + length(docetaxel.labels$slope_breast_only))
}

partition_var$cpp_breast_only <- temp.cpp

## patients
temp.pp <- foreach (training_amount.cpp = seq(from = 2, to = 22, by = 1)) %do% {      
  generate_random_partition.pp_var(input_labels_patient = docetaxel.labels$patient, 
                                   training_amount = training_amount.cpp)
}

# double check
for (temp.ind in 1:22) {
  stopifnot(length(temp.pp[[temp.ind]][[100]]$test_index) == length(docetaxel.labels$patient) - temp.ind - 1)
  stopifnot(length(temp.pp[[temp.ind]][[100]]$training_index.single) == temp.ind)
}

partition_var$pp <- temp.pp


### UNCOMMENT if you want to save the homogenized dataset
save(docetaxel, docetaxel.labels, sampleinfo.cgp, partition, feature.l1000, partition_var, file = "../WS/docetaxel_data.RData")


