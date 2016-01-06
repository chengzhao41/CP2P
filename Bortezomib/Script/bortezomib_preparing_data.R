# Load Libraries ----------------------------------------------------------
library("doParallel")
library("sva")

source('Common/preparing_data_helper.R')
source('Common/drug_cut/callingWaterfall.R')
source('Common/drug_cut/distancePointLine.R')
source('Common/drug_cut/distancePointSegment.R')
source('Common/comGENE.R')
source("Common/generate_random_partition.R")
source("Common/ordering_by_similarity.R")

# Load CGP Cell line data --------------------------------------------
load("CGP/cdrug2_cgp_ccle_all.RData")

# Get Slope Summary Statistics --------------------------------------------
bortezomib.labels <- list()

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

# Get IC50 and AUC Summary Statistics -------------------------------------
temp.drug_ind <- which(druginfo.cgp$drug.name == "BORTEZOMIB")

temp.ic50 <- drugpheno.cgp$IC50[, temp.drug_ind]
temp.ic50_ind <- which(!is.na(temp.ic50))
temp.ic50 <- temp.ic50[temp.ic50_ind]
length(temp.ic50)

temp.auc <- drugpheno.cgp$AUC[, temp.drug_ind]
temp.auc_ind <- which(!is.na(temp.auc))
temp.auc <- temp.auc[temp.auc_ind]
length(temp.auc)

# Binarize Responses ------------------------------------------------------
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

# Get the gene expression data --------------------------------------------
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

# check agreement among the 3 different summary statistics ----------------
temp.ind <- match(names(bortezomib.labels$slope), names(bortezomib.labels$AUC))
stopifnot(names(bortezomib.labels$AUC)[temp.ind] == names(bortezomib.labels$slope))
print(paste(sum(bortezomib.labels$AUC[temp.ind] == bortezomib.labels$slope), "/", sum(!is.na(temp.ind))))
# 266 / 311 for slope vs. AUC

temp.ind <- match(names(bortezomib.labels$slope), names(bortezomib.labels$IC50))
stopifnot(names(bortezomib.labels$IC50)[temp.ind] == names(bortezomib.labels$slope))
print(paste(sum(bortezomib.labels$IC50[temp.ind] == bortezomib.labels$slope), "/", sum(!is.na(temp.ind))))
# 271 / 311 for slope vs. IC50

temp.ind <- match(names(bortezomib.labels$AUC), names(bortezomib.labels$IC50))
stopifnot(names(bortezomib.labels$IC50)[temp.ind] == names(bortezomib.labels$AUC))
print(paste(sum(bortezomib.labels$IC50[temp.ind] == bortezomib.labels$AUC), "/", sum(!is.na(temp.ind))))
# 286 / 313 for AUC vs. IC50 

# Using sva to harmonize across different tissue types --------------------
#show_pca(input_data = bortezomib$cgp_slope, label = bortezomib.labels$slope)
#show_pca(input_data = bortezomib$cgp_AUC, label = bortezomib.labels$AUC)
#show_pca(input_data = bortezomib$cgp_IC50, label = bortezomib.labels$IC50)
mean(scale(bortezomib$cgp_slope)) #1.066058e-17
mean(scale(bortezomib$cgp_AUC)) #-2.85792e-18
mean(scale(bortezomib$cgp_IC50)) #-2.85792e-18

# slope
bortezomib$cgp_slope.sva <- sva_combine(batch = sampleinfo.cgp$tissue.type[bortezomib.labels$slope_ind],  
                                        label = bortezomib.labels$slope, input_data = scale(bortezomib$cgp_slope))
mean(bortezomib$cgp_slope.sva) #1.131253e-17
#show_pca(input_data = bortezomib$cgp_slope.sva, label = bortezomib.labels$slope)
#show_pca(input_data = bortezomib$cgp_slope.sva, label = sampleinfo.cgp$tissue.type[bortezomib.labels$slope_ind])

# IC50
bortezomib$cgp_IC50.sva <- sva_combine(batch = sampleinfo.cgp$tissue.type[bortezomib.labels$IC50_ind],  
                                       label = bortezomib.labels$IC50, input_data = scale(bortezomib$cgp_IC50))
mean(bortezomib$cgp_IC50.sva) #-4.927404e-18
# show_pca(input_data = bortezomib$cgp_IC50.sva, label = bortezomib.labels$IC50)
# show_pca(input_data = bortezomib$cgp_IC50.sva, label = sampleinfo.cgp$tissue.type[bortezomib.labels$IC50_ind])

# AUC
bortezomib$cgp_AUC.sva <- sva_combine(batch = sampleinfo.cgp$tissue.type[bortezomib.labels$AUC_ind],  
                                      label = bortezomib.labels$AUC, input_data = scale(bortezomib$cgp_AUC))
mean(bortezomib$cgp_AUC.sva) #-5.764707e-18
# show_pca(input_data = bortezomib$cgp_AUC.sva, label = bortezomib.labels$AUC)
# show_pca(input_data = bortezomib$cgp_AUC.sva, label = sampleinfo.cgp$tissue.type[bortezomib.labels$AUC_ind])

# Get Patient expression data ---------------------------------------------
load("Bortezomib/WS/bortezomib.patient.RData")

bortezomib$patient.combat <- scale(t(bortezomib.patient_ComBat))
mean(bortezomib$patient.combat) #-2.371065e-18
bortezomib.labels$patient <- binaryResponse == 1
table(bortezomib.labels$patient)
#FALSE  TRUE 
#84    85 

# Using slope labels for SVA  ---------------------------------
bortezomib.labels$slope_combined <- c(bortezomib.labels$patient, bortezomib.labels$slope)
stopifnot(substring(colnames(bortezomib$cgp_slope.sva), 8) == annot.ge.cgp$EntrezGene.ID)
colnames(bortezomib$cgp_slope.sva) <- annot.ge.cgp$symbol
colnames(bortezomib$cgp_slope) <- annot.ge.cgp$symbol

temp.data <- comGENE(bortezomib$patient.combat, bortezomib$cgp_slope.sva)

mean(temp.data[[1]]) #-2.442217e-18
mean(temp.data[[2]]) #1.236016e-17
dim(temp.data[[1]])
dim(temp.data[[2]])

bortezomib.labels$slope_combined.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("cgp", dim(temp.data[[2]])[1]))
stopifnot(names(bortezomib.labels$slope_combined) == rownames(rbind(temp.data[[1]], scale(temp.data[[2]]))))
bortezomib$combined_slope.sva <- sva_combine(batch = bortezomib.labels$slope_combined.source,
                                             label = bortezomib.labels$slope_combined,
                                             input_data = rbind(temp.data[[1]], temp.data[[2]]), n.sv=2)
mean(bortezomib$combined_slope.sva) #8.868008e-18

show_pca(input_data = bortezomib$combined_slope.sva, label = bortezomib.labels$slope_combined)
show_pca(input_data = bortezomib$combined_slope.sva, label = bortezomib.labels$slope_combined.source)
rm(temp.data)

# Using IC50 labels for SVA  ---------------------------------
bortezomib.labels$IC50_combined <- c(bortezomib.labels$patient, bortezomib.labels$IC50)
stopifnot(substring(colnames(bortezomib.labels$IC50), 8) == annot.ge.cgp$EntrezGene.ID)
colnames(bortezomib$cgp_IC50.sva) <- annot.ge.cgp$symbol

temp.data <- comGENE(bortezomib$patient.combat, bortezomib$cgp_IC50.sva)
mean(temp.data[[1]]) #-2.442187e-18
mean(temp.data[[2]]) #-6.538526e-18
dim(temp.data[[1]])
dim(temp.data[[2]])

bortezomib.labels$IC50_combined.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("cgp", dim(temp.data[[2]])[1]))

stopifnot(names(bortezomib.labels$IC50_combined.source) == rownames(rbind(temp.data[[1]], temp.data[[2]])))

bortezomib$combined_IC50.sva <- sva_combine(batch = bortezomib.labels$IC50_combined.source,
                                            label = bortezomib.labels$IC50_combined, 
                                            input_data = rbind(temp.data[[1]], temp.data[[2]]), n.sv = 2)
mean(bortezomib$combined_IC50.sva) #-1.05882e-17
show_pca(input_data = bortezomib$combined_IC50.sva, label = bortezomib.labels$IC50_combined)
show_pca(input_data = bortezomib$combined_IC50.sva, label = bortezomib.labels$IC50_combined.source)
rm(temp.data)

# Using AUC labels for SVA  ---------------------------------
bortezomib.labels$AUC_combined <- c(bortezomib.labels$patient, bortezomib.labels$AUC)
stopifnot(substring(colnames(bortezomib.labels$AUC), 8) == annot.ge.cgp$EntrezGene.ID)
colnames(bortezomib$cgp_AUC.sva) <- annot.ge.cgp$symbol

temp.data <- comGENE(bortezomib$patient.combat, bortezomib$cgp_AUC.sva)
mean(temp.data[[1]]) #-2.442187e-18
mean(temp.data[[2]]) #-7.117706e-18

bortezomib.labels$AUC_combined.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("cgp", dim(temp.data[[2]])[1]))
bortezomib$combined_AUC.sva <- sva_combine(batch = bortezomib.labels$AUC_combined.source,
                                           label = bortezomib.labels$AUC_combined, 
                                           input_data = rbind(temp.data[[1]], temp.data[[2]]), n.sv = 2)
mean(bortezomib$combined_AUC.sva) #-1.10284e-17 

show_pca(input_data = bortezomib$combined_AUC.sva, label = bortezomib.labels$AUC_combined)
show_pca(input_data = bortezomib$combined_AUC.sva, label = bortezomib.labels$AUC_combined.source)
rm(temp.data)

# find the common cell lines that exists for all 3 response labels --------
bortezomib.labels$IC50_ind_common = intersect(
  which(names(bortezomib.labels$IC50) %in% names(bortezomib.labels$AUC)),
  which(names(bortezomib.labels$IC50) %in% names(bortezomib.labels$slope)))
bortezomib.labels$AUC_ind_common = intersect(
  which(names(bortezomib.labels$AUC) %in% names(bortezomib.labels$IC50)),
  which(names(bortezomib.labels$AUC) %in% names(bortezomib.labels$slope)))
bortezomib.labels$slope_ind_common = intersect(
  which(names(bortezomib.labels$slope) %in% names(bortezomib.labels$IC50)),
  which(names(bortezomib.labels$slope) %in% names(bortezomib.labels$AUC)))

stopifnot(names(bortezomib.labels$IC50[bortezomib.labels$IC50_ind_common]) 
          == names(bortezomib.labels$AUC[bortezomib.labels$AUC_ind_common]))
stopifnot(names(bortezomib.labels$IC50[bortezomib.labels$IC50_ind_common]) 
          == names(bortezomib.labels$slope[bortezomib.labels$slope_ind_common]))

# get l1000 features ------------------------------------------------------
stopifnot(colnames(bortezomib$cgp_IC50.sva) == colnames(bortezomib$cgp_AUC.sva))
stopifnot(colnames(bortezomib$cgp_slope.sva) == colnames(bortezomib$cgp_AUC.sva))

Landmark_Genes_n978 <- read.csv("Common/Landmark_Genes_n978.csv")

feature.l1000 <- list()
feature.l1000$cp <- which(colnames(bortezomib$combined_slope.sva) %in% Landmark_Genes_n978$Gene.Symbol)
feature.l1000$pp <- which(colnames(bortezomib$patient.combat) %in% Landmark_Genes_n978$Gene.Symbol)

# create partitions for 10 to 150 patients using all cell lines ------------------------
partition <- list()

input.labels_cell_lines = list()
input.labels_cell_lines$slope = bortezomib.labels$slope[bortezomib.labels$slope_ind_common]
input.labels_cell_lines$IC50 = bortezomib.labels$IC50[bortezomib.labels$IC50_ind_common]
input.labels_cell_lines$AUC = bortezomib.labels$AUC[bortezomib.labels$AUC_ind_common]

input.cell_line_order = list()
input.cell_line_order$slope = 1:length(input.labels_cell_lines$slope)
input.cell_line_order$IC50 = 1:length(input.labels_cell_lines$IC50)
input.cell_line_order$AUC = 1:length(input.labels_cell_lines$AUC)

stopifnot(input.cell_line_order$slope == input.cell_line_order$IC50)  
stopifnot(input.cell_line_order$AUC == input.cell_line_order$IC50)  

cell_lines_all <- foreach(input.training_amount.p = seq(from = 10, to = 150, by = 10)
                          , .errorhandling = "stop") %dopar% {
  
  generate_random_partition(labels_cell_lines = input.labels_cell_lines, 
                            labels_patient = bortezomib.labels$patient,
                            training_amount.p = input.training_amount.p,
                            leave_one_out = FALSE,
                            cell_line_order = input.cell_line_order,
                            training_amount.c = length(input.cell_line_order$AUC)
                            )
}

partition$cell_lines_all <- cell_lines_all

# order cell lines by similarity using 100 and 50 patients ------------------------------------------
input.labels_cell_lines = list()
input.labels_cell_lines$slope = bortezomib.labels$slope
input.labels_cell_lines$IC50 = bortezomib.labels$IC50
input.labels_cell_lines$AUC = bortezomib.labels$AUC

input.cell_line_order = list()

temp.slope_combined_ind_common = c(1:length(bortezomib.labels$patient), 
                                   length(bortezomib.labels$patient) + bortezomib.labels$slope_ind_common)
stopifnot(max(table(temp.slope_combined_ind_common)) == 1)

temp.auc_combined_ind_common = c(1:length(bortezomib.labels$patient), 
                                   length(bortezomib.labels$patient) + bortezomib.labels$AUC_ind_common)
stopifnot(max(table(temp.auc_combined_ind_common)) == 1)

temp.ic50_combined_ind_common = c(1:length(bortezomib.labels$patient), 
                                   length(bortezomib.labels$patient) + bortezomib.labels$IC50_ind_common)
stopifnot(max(table(temp.ic50_combined_ind_common)) == 1)

input.cell_line_order$slope =
  order_by_similarity(data = bortezomib$combined_slope.sva[temp.slope_combined_ind_common, ],
                      labels = bortezomib.labels$slope_combined[temp.slope_combined_ind_common],
                      source_labels = bortezomib.labels$slope_combined.source[temp.slope_combined_ind_common])

input.cell_line_order$AUC =
  order_by_similarity(data = bortezomib$combined_AUC.sva[temp.auc_combined_ind_common, ],
                      labels = bortezomib.labels$AUC_combined[temp.auc_combined_ind_common],
                      source_labels = bortezomib.labels$AUC_combined.source[temp.auc_combined_ind_common])

input.cell_line_order$IC50 =
  order_by_similarity(data = bortezomib$combined_IC50.sva[temp.ic50_combined_ind_common, ],
                      labels = bortezomib.labels$IC50_combined[temp.ic50_combined_ind_common],
                      source_labels = bortezomib.labels$IC50_combined.source[temp.ic50_combined_ind_common])

stopifnot(length(input.cell_line_order$slope) == length(input.cell_line_order$IC50))  
stopifnot(length(input.cell_line_order$AUC) == length(input.cell_line_order$IC50))  

patient_100 <- foreach(input.training_amount.c = seq(from = 20, to = 300, by = 20), .errorhandling = "stop") %dopar% {
                            generate_random_partition(labels_cell_lines = input.labels_cell_lines, 
                                                      labels_patient = bortezomib.labels$patient,
                                                      training_amount.p = 100,
                                                      leave_one_out = FALSE,
                                                      cell_line_order = input.cell_line_order,
                                                      training_amount.c = input.training_amount.c,
                                                      input_partition = partition$cell_lines_all[[10]])
  }
partition$patient_100 <- patient_100

patient_50 <- foreach(input.training_amount.c = seq(from = 20, to = 300, by = 20), .errorhandling = "stop") %dopar% {
  generate_random_partition(labels_cell_lines = input.labels_cell_lines, 
                            labels_patient = bortezomib.labels$patient,
                            training_amount.p = 50,
                            leave_one_out = FALSE,
                            cell_line_order = input.cell_line_order,
                            training_amount.c = input.training_amount.c,
                            input_partition = partition$cell_lines_all[[5]])
}
partition$patient_50 <- patient_50

# save worksapce ----------------------------------------------------------
save(sampleinfo.cgp, bortezomib, bortezomib.labels, feature.l1000, partition,
     file = "Bortezomib/WS/bortezomib_data.RData")
