library("doParallel")
library("sva")

source('Bortezomib/Script/generate_random_partition.R')
source('Common/preparing_data_helper.R')
source('Common/drug_cut/callingWaterfall.R')
source('Common/drug_cut/distancePointLine.R')
source('Common/drug_cut/distancePointSegment.R')
source('Common/comGENE.R')
source("Common/generate_random_partition.R")

###

load("CGP/cdrug2_cgp_ccle_all.RData")

bortezomib.labels <- list()

### Get bortezomib response labels from Slope Summary Statistics ###
cgp_sensitivity_1 <- read.csv("CGP/cgp_sensitivity_1.csv")

temp.bortezomib_ind <- which(cgp_sensitivity_1$drug.name == "Bortezomib")
cgp_sensitivity_1 <- cgp_sensitivity_1[temp.bortezomib_ind, ]
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

bortezomib.labels$slope <- temp.labels_slope[temp.slope_ind]
bortezomib.labels$slope_ind <- which(rownames(data.ge.cgp) %in% names(temp.labels_slope))
length(bortezomib.labels$slope_ind)

stopifnot(names(bortezomib.labels$slope) == rownames(data.ge.cgp)[bortezomib.labels$slope_ind])
### END ###

### Get bortezomib response labels from IC50 and AUC Summary Statistics ###
temp.drug_ind <- which(druginfo.cgp$drug.name == "BORTEZOMIB")

temp.ic50 <- drugpheno.cgp$IC50[, temp.drug_ind]
temp.ic50_ind <- which(!is.na(temp.ic50))
temp.ic50 <- temp.ic50[temp.ic50_ind]
length(temp.ic50)

temp.auc <- drugpheno.cgp$AUC[, temp.drug_ind]
temp.auc_ind <- which(!is.na(temp.auc))
temp.auc <- temp.auc[temp.auc_ind]
length(temp.auc)

# binarize the responses
temp.response <- callingWaterfall(temp.ic50, type="IC50")

temp.response_auc <- callingWaterfall(temp.auc, type="AUC")

temp.label_IC50 <- temp.response != "resistant"
table(temp.label_IC50)
names(temp.label_IC50) <- names(temp.response)

temp.label_auc <- temp.response_auc != "resistant"
table(temp.label_auc)
names(temp.label_auc) <- names(temp.response_auc)

bortezomib.labels$IC50_ind <- which(rownames(data.ge.cgp) %in% names(temp.label_IC50))
bortezomib.labels$AUC_ind <- which(rownames(data.ge.cgp) %in% names(temp.label_auc))

stopifnot(names(temp.label_IC50) == rownames(data.ge.cgp)[bortezomib.labels$IC50_ind])
bortezomib.labels$IC50 <- temp.label_IC50
stopifnot(names(temp.label_auc) == rownames(data.ge.cgp)[bortezomib.labels$AUC_ind])
bortezomib.labels$AUC <- temp.label_auc
### END ###

### Getting the expression values
## get the gene expression data
bortezomib <- list()
bortezomib$cgp_slope <- data.ge.cgp[bortezomib.labels$slope_ind,  ]
bortezomib$cgp_IC50 <- data.ge.cgp[bortezomib.labels$IC50_ind,  ]
bortezomib$cgp_AUC <- data.ge.cgp[bortezomib.labels$AUC_ind,  ]

temp <- remove_features(input_data = bortezomib$cgp_slope, input_ind = bortezomib.labels$slope_ind,
                        input_label = bortezomib.labels$slope)
bortezomib$cgp_slope <- temp$data
bortezomib.labels$slope_ind <- temp$ind
bortezomib.labels$slope <- temp$label

rm(temp)
temp <- remove_features(input_data = bortezomib$cgp_IC50, input_ind = bortezomib.labels$IC50_ind,
                        input_label = bortezomib.labels$IC50)
bortezomib$cgp_IC50 <- temp$data
bortezomib.labels$IC50_ind <- temp$ind
bortezomib.labels$IC50 <- temp$label

rm(temp)
temp <- remove_features(input_data = bortezomib$cgp_AUC, input_ind = bortezomib.labels$AUC_ind,
                        input_label = bortezomib.labels$AUC)
bortezomib$cgp_AUC <- temp$data
bortezomib.labels$AUC_ind <- temp$ind
bortezomib.labels$AUC <- temp$label

## check agreement among the 3 different summary statistics
temp.ind <- match(names(bortezomib.labels$slope), names(bortezomib.labels$AUC))
stopifnot(names(bortezomib.labels$AUC)[temp.ind] == names(bortezomib.labels$slope))
sum(bortezomib.labels$AUC[temp.ind] == bortezomib.labels$slope)
sum(!is.na(temp.ind))
# 266 / 311 for slope vs. AUC
temp.ind <- match(names(bortezomib.labels$slope), names(bortezomib.labels$IC50))
sum(bortezomib.labels$IC50[temp.ind] == bortezomib.labels$slope)
sum(!is.na(temp.ind))
# 271 / 311 for slope vs. IC50
temp.ind <- match(names(bortezomib.labels$AUC), names(bortezomib.labels$IC50))
sum(bortezomib.labels$IC50[temp.ind] == bortezomib.labels$AUC)
sum(!is.na(temp.ind))
# 286 / 313 for AUC vs. IC50 

#save(bortezomib, bortezomib.labels, sampleinfo.cgp, file = "Bortezomib/WS/bortezomib_data.RData")
### END

# Using sva to harmonize across different tissue types
# show_pca(input_data = bortezomib$cgp_slope, label = bortezomi b.labels$slope)
# show_pca(input_data = bortezomib$cgp_AUC, label = bortezomib.labels$AUC)
# show_pca(input_data = bortezomib$cgp_IC50, label = bortezomib.labels$IC50)
# slope
bortezomib$cgp_slope.sva <- sva_combine(batch = sampleinfo.cgp$tissue.type[bortezomib.labels$slope_ind],  
                                        label = bortezomib.labels$slope, input_data = scale(bortezomib$cgp_slope))

#show_pca(input_data = bortezomib$cgp_slope.sva, label = bortezomib.labels$slope)
#show_pca(input_data = bortezomib$cgp_slope.sva, label = sampleinfo.cgp$tissue.type[bortezomib.labels$slope_ind])

# IC50
bortezomib$cgp_IC50.sva <- sva_combine(batch = sampleinfo.cgp$tissue.type[bortezomib.labels$IC50_ind],  
                                       label = bortezomib.labels$IC50, input_data = scale(bortezomib$cgp_IC50))
# show_pca(input_data = bortezomib$cgp_IC50.sva, label = bortezomib.labels$IC50)
# show_pca(input_data = bortezomib$cgp_IC50.sva, label = sampleinfo.cgp$tissue.type[bortezomib.labels$IC50_ind])

# AUC
bortezomib$cgp_AUC.sva <- sva_combine(batch = sampleinfo.cgp$tissue.type[bortezomib.labels$AUC_ind],  
                                      label = bortezomib.labels$AUC, input_data = scale(bortezomib$cgp_AUC))
# show_pca(input_data = bortezomib$cgp_AUC.sva, label = bortezomib.labels$AUC)
# show_pca(input_data = bortezomib$cgp_AUC.sva, label = sampleinfo.cgp$tissue.type[bortezomib.labels$AUC_ind])

# Using sva to harmonize across different tissue types
#save("Bortezomib/WS/bortezomib_data.RData")

### Get Patient expression data ###
##load("Bortezomib/WS/bortezomib_data.RData")

# get patient data 
load("Bortezomib/WS/bortezomib.patient.RData")

bortezomib$patient.combat <- scale(t(bortezomib.patient_ComBat))
bortezomib.labels$patient <- binaryResponse == 1
table(bortezomib.labels$patient)

### Slopes
# Using sva to harmonize patients and cell lines
bortezomib.labels$slope_combined <- c(bortezomib.labels$slope, bortezomib.labels$patient)
stopifnot(substring(colnames(bortezomib$cgp_slope.sva), 8) == annot.ge.cgp$EntrezGene.ID)
colnames(bortezomib$cgp_slope.sva) <- annot.ge.cgp$symbol
colnames(bortezomib$cgp_slope) <- annot.ge.cgp$symbol

temp.data <- comGENE(bortezomib$cgp_slope.sva, bortezomib$patient.combat)
#temp.data <- comGENE(scale(bortezomib$cgp_slope), bortezomib$patient.combat)
temp.data <- comGENE(bortezomib$cgp_slope.sva, data$patient.combat)

mean(temp.data[[1]])
mean(temp.data[[2]])
dim(temp.data[[1]])
dim(temp.data[[2]])

bortezomib.labels$slope_combined.source <- c(rep("CGP", dim(temp.data[[1]])[1]), rep("Patient", dim(temp.data[[2]])[1]))

bortezomib$combined_slope.sva <- sva_combine(batch = bortezomib.labels$slope_combined.source,
                                             label = bortezomib.labels$slope_combined,
                                             input_data = rbind(temp.data[[1]], scale(temp.data[[2]])), n.sv=2)
mean(bortezomib$combined_slope.sva)

show_pca(input_data = bortezomib$combined_slope.sva, label = bortezomib.labels$slope_combined)
show_pca(input_data = bortezomib$combined_slope.sva, label = bortezomib.labels$slope_combined.source)
rm(temp.data)
#save(bortezomib, bortezomib.labels, sampleinfo.cgp, file = "Bortezomib/WS/bortezomib_data.RData")

## IC50
# Using sva to harmonize patients and cell lines
bortezomib.labels$IC50_combined <- c(bortezomib.labels$IC50, bortezomib.labels$patient)
stopifnot(substring(colnames(bortezomib.labels$IC50), 8) == annot.ge.cgp$EntrezGene.ID)
colnames(bortezomib$cgp_IC50.sva) <- annot.ge.cgp$symbol

temp.data <- comGENE(scale(bortezomib$cgp_IC50.sva), bortezomib$patient.combat)
mean(temp.data[[1]])
mean(temp.data[[2]])
dim(temp.data[[1]])
dim(temp.data[[2]])

bortezomib.labels$IC50_combined.source <- c(rep("cgp", dim(temp.data[[1]])[1]), rep("patient", dim(temp.data[[2]])[1]))
bortezomib$combined_IC50.sva <- sva_combine(batch = bortezomib.labels$IC50_combined.source,
                                            label = bortezomib.labels$IC50_combined, 
                                            input_data = rbind(temp.data[[1]], temp.data[[2]]), n.sv = 2)

show_pca(input_data = bortezomib$combined_IC50.sva, label = bortezomib.labels$IC50_combined)
show_pca(input_data = bortezomib$combined_IC50.sva, label = bortezomib.labels$IC50_combined.source)
rm(temp.data)
#save(bortezomib, bortezomib.labels, sampleinfo.cgp, file = "Bortezomib/WS/bortezomib_data.RData")

## AUC
# Using sva to harmonize patients and cell lines
bortezomib.labels$AUC_combined <- c(bortezomib.labels$AUC, bortezomib.labels$patient)
stopifnot(substring(colnames(bortezomib.labels$AUC), 8) == annot.ge.cgp$EntrezGene.ID)
colnames(bortezomib$cgp_AUC.sva) <- annot.ge.cgp$symbol

temp.data <- comGENE(bortezomib$cgp_AUC.sva, bortezomib$patient.combat)
mean(temp.data[[1]])
mean(temp.data[[2]])

bortezomib.labels$AUC_combined.source <- c(rep("cgp", dim(temp.data[[1]])[1]), rep("patient", dim(temp.data[[2]])[1]))
bortezomib$combined_AUC.sva <- sva_combine(batch = bortezomib.labels$AUC_combined.source,
                                           label = bortezomib.labels$AUC_combined, 
                                           input_data = rbind(temp.data[[1]], temp.data[[2]]), n.sv = 2)

show_pca(input_data = bortezomib$combined_AUC.sva, label = bortezomib.labels$AUC_combined)
show_pca(input_data = bortezomib$combined_AUC.sva, label = bortezomib.labels$AUC_combined.source)
rm(temp.data)

#save(bortezomib, bortezomib.labels, sampleinfo.cgp, file = "Bortezomib/WS/bortezomib_data.RData")

### get the partitioning
partition <- list()

temp.test_amount <- list(cpp = round(0.4 * length(bortezomib.labels$patient)), cc = round(0.2 * length(bortezomib.labels$slope)))
partition$slope <- generate_random_partition(input_labels_cell_lines = bortezomib.labels$slope, 
                                             input_labels_patient = bortezomib.labels$patient, 
                                             test_amount = temp.test_amount)

temp.test_amount <- list(cpp = round(0.4 * length(bortezomib.labels$patient)), cc = round(0.2 * length(bortezomib.labels$IC50)))
partition$IC50 <- generate_random_partition(input_labels_cell_lines = bortezomib.labels$IC50, 
                                            input_labels_patient = bortezomib.labels$patient, 
                                            test_amount = temp.test_amount)

temp.test_amount <- list(cpp = round(0.4 * length(bortezomib.labels$patient)), cc = round(0.2 * length(bortezomib.labels$AUC)))
partition$AUC <- generate_random_partition(input_labels_cell_lines = bortezomib.labels$AUC, 
                                           input_labels_patient = bortezomib.labels$patient, 
                                           test_amount = temp.test_amount)

#save(bortezomib, bortezomib.labels, sampleinfo.cgp, partition, file = "Bortezomib/WS/bortezomib_data.RData")

### get l1000 features
stopifnot(colnames(bortezomib$cgp_IC50.sva) == colnames(bortezomib$cgp_AUC.sva))
stopifnot(colnames(bortezomib$cgp_slope.sva) == colnames(bortezomib$cgp_AUC.sva))

Landmark_Genes_n978 <- read.csv("Common/Landmark_Genes_n978.csv")

feature.l1000 <- list()
feature.l1000$cc <- which(colnames(bortezomib$cgp_slope.sva) %in% Landmark_Genes_n978$Gene.Symbol) # not used
feature.l1000$cp <- which(colnames(bortezomib$combined_slope.sva) %in% Landmark_Genes_n978$Gene.Symbol)
feature.l1000$pp <- which(colnames(bortezomib$patient.combat) %in% Landmark_Genes_n978$Gene.Symbol)

#save(bortezomib, bortezomib.labels, sampleinfo.cgp, partition, feature.l1000, file = "Bortezomib/WS/bortezomib_data.RData")

### create partitions for varying number of patients
partition_var <- list()

temp.cp <- foreach (training_amount.cp = seq(from = 10, to = 300, by = 10)) %do% {    
  
  temp.slope <- generate_random_partition.cp_var(input_labels_cell_lines = bortezomib.labels$slope, 
                                                 input_labels_patient = bortezomib.labels$patient, 
                                                 training_amount = training_amount.cp)
  
  temp.IC50 <- generate_random_partition.cp_var(input_labels_cell_lines = bortezomib.labels$IC50, 
                                                input_labels_patient = bortezomib.labels$patient, 
                                                training_amount = training_amount.cp)
  
  temp.AUC <- generate_random_partition.cp_var(input_labels_cell_lines = bortezomib.labels$AUC, 
                                               input_labels_patient = bortezomib.labels$patient, 
                                               training_amount = training_amount.cp)
  list(slope = temp.slope, IC50 = temp.IC50, AUC = temp.AUC)
}

# double check
for (temp.ind in 1:30) {
  stopifnot(length(temp.cp[[temp.ind]]$slope[[100]]$test_index) == length(bortezomib.labels$patient))
  stopifnot(length(temp.cp[[temp.ind]]$slope[[100]]$training_index.single) == temp.ind * 10)
}

partition_var$cp <- temp.cp

temp.cpp <- foreach (training_amount.cpp = seq(from = 10, to = 150, by = 10)) %do% {    
  
  temp.slope <- generate_random_partition.cpp_var(input_labels_cell_lines = bortezomib.labels$slope, 
                                                  input_labels_patient = bortezomib.labels$patient, 
                                                  training_amount = training_amount.cpp)
  
  temp.IC50 <- generate_random_partition.cpp_var(input_labels_cell_lines = bortezomib.labels$IC50, 
                                                 input_labels_patient = bortezomib.labels$patient, 
                                                 training_amount = training_amount.cpp)
  
  temp.AUC <- generate_random_partition.cpp_var(input_labels_cell_lines = bortezomib.labels$AUC, 
                                                input_labels_patient = bortezomib.labels$patient, 
                                                training_amount = training_amount.cpp)
  list(slope = temp.slope, IC50 = temp.IC50, AUC = temp.AUC)
}

partition_var$cpp <- temp.cpp

# double check
for (temp.ind in 1:15) {
  stopifnot(length(temp.cpp[[temp.ind]]$slope[[100]]$test_index) == length(bortezomib.labels$patient) - temp.ind * 10 )
  stopifnot(length(temp.cpp[[temp.ind]]$slope[[100]]$training_index.single) == temp.ind * 10 + length(bortezomib.labels$slope))
}

## patients
temp.pp <- foreach (training_amount.cpp = seq(from = 10, to = 150, by = 10)) %do% {      
  generate_random_partition.pp_var(input_labels_patient = bortezomib.labels$patient, 
                                                  training_amount = training_amount.cpp)
}

# double check
for (temp.ind in 1:15) {
  stopifnot(length(temp.pp[[temp.ind]][[100]]$test_index) == length(bortezomib.labels$patient) - temp.ind * 10 )
  stopifnot(length(temp.pp[[temp.ind]][[100]]$training_index.single) == temp.ind * 10)
}

partition_var$pp <- temp.pp

#### get partitions for 100 patients and varying number of cell lines in training set
temp.cpp <- foreach (cell_line_training_amount = seq(from = 10, to = 310, by = 10)) %do% {    
  
  temp.slope <- generate_random_partition.cpp_var2(input_labels_cell_lines = bortezomib.labels$slope, 
                                                  input_labels_patient = bortezomib.labels$patient, 
                                                  cell_line_training_amount = cell_line_training_amount,
                                                  patient_training_amount = 100)
  
  temp.IC50 <- generate_random_partition.cpp_var2(input_labels_cell_lines = bortezomib.labels$IC50, 
                                                 input_labels_patient = bortezomib.labels$patient, 
                                                 cell_line_training_amount = cell_line_training_amount,
                                                 patient_training_amount = 100)
  
  temp.AUC <- generate_random_partition.cpp_var2(input_labels_cell_lines = bortezomib.labels$AUC, 
                                                input_labels_patient = bortezomib.labels$patient, 
                                                cell_line_training_amount = cell_line_training_amount,
                                                patient_training_amount = 100)
  list(slope = temp.slope, IC50 = temp.IC50, AUC = temp.AUC)
}

partition_var$cVar_p100_p <- temp.cpp

#### get partitions for 100 patients and varying number of cell lines in training set
temp.cpp <- foreach (cell_line_training_amount = seq(from = 10, to = 310, by = 10)) %do% {    
  
  temp.slope <- generate_random_partition.cpp_var2(input_labels_cell_lines = bortezomib.labels$slope, 
                                                   input_labels_patient = bortezomib.labels$patient, 
                                                   cell_line_training_amount = cell_line_training_amount,
                                                   patient_training_amount = 50)
  
  temp.IC50 <- generate_random_partition.cpp_var2(input_labels_cell_lines = bortezomib.labels$IC50, 
                                                  input_labels_patient = bortezomib.labels$patient, 
                                                  cell_line_training_amount = cell_line_training_amount,
                                                  patient_training_amount = 50)
  
  temp.AUC <- generate_random_partition.cpp_var2(input_labels_cell_lines = bortezomib.labels$AUC, 
                                                 input_labels_patient = bortezomib.labels$patient, 
                                                 cell_line_training_amount = cell_line_training_amount,
                                                 patient_training_amount = 50)
  list(slope = temp.slope, IC50 = temp.IC50, AUC = temp.AUC)
}

partition_var$cVar_p50_p <- temp.cpp



#save(bortezomib, bortezomib.labels, sampleinfo.cgp, partition, feature.l1000, partition_var, file = "Bortezomib/WS/bortezomib_data.RData")
