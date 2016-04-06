# Load Libraries ----------------------------------------------------------
library("doParallel")
library("sva")

source('Common/preparing_data_helper.R')
source('Common/comGENE.R')
source("Common/generate_random_partition.R")
source("Common/ordering_by_similarity.R")

# Load CGP Cell line data --------------------------------------------
load("Bortezomib/WS/bortezomib_gdsc.RData")

# sanity check to make sure we have the right version of data
mean(bortezomib$gdsc_slope) # 6.293544
mean(bortezomib$gdsc_AUC) # 6.293544
mean(bortezomib$gdsc_IC50) # 6.293544
dim(bortezomib$gdsc_slope) # 313 11833
dim(bortezomib$gdsc_AUC) # 313 11833
dim(bortezomib$gdsc_IC50) # 313 11833
table(bortezomib.labels$slope)
# FALSE  TRUE 
# 123   190 
table(bortezomib.labels$AUC)
# FALSE  TRUE 
# 123   190 
table(bortezomib.labels$IC50)
# FALSE  TRUE 
# 73   240
stopifnot(sum(bortezomib.labels$slope != bortezomib.labels$AUC) == 6)

# Using sva to harmonize across different tissue types --------------------

# slope
# before sva
# show_pca(input_data = bortezomib$gdsc_slope, label = bortezomib.labels$slope)
show_pca(input_data = scale(bortezomib$gdsc_slope), label = sampleinfo.gdsc$tissue.type[bortezomib.labels$slope_ind])
# sva
bortezomib$gdsc_slope.sva <- sva_combine(batch = sampleinfo.gdsc$tissue.type[bortezomib.labels$slope_ind],  
                                        label = bortezomib.labels$slope, 
                                        input_data = scale(bortezomib$gdsc_slope),
                                        n.sv = 3)
mean(bortezomib$gdsc_slope.sva) # 5.587487e-18
# after sva
# show_pca(input_data = bortezomib$gdsc_slope.sva, label = bortezomib.labels$slope)
show_pca(input_data = bortezomib$gdsc_slope.sva, label = sampleinfo.gdsc$tissue.type[bortezomib.labels$slope_ind])

# IC50
# before sva
# show_pca(input_data = bortezomib$gdsc_IC50, label = bortezomib.labels$IC50)
show_pca(input_data = scale(bortezomib$gdsc_IC50), label = sampleinfo.gdsc$tissue.type[bortezomib.labels$IC50_ind])
# sva
bortezomib$gdsc_IC50.sva <- sva_combine(batch = sampleinfo.gdsc$tissue.type[bortezomib.labels$IC50_ind],  
                                       label = bortezomib.labels$IC50, 
                                       input_data = scale(bortezomib$gdsc_IC50),
                                       n.sv = 3)
mean(bortezomib$gdsc_IC50.sva) # 2.748525e-18
# after sva
# show_pca(input_data = bortezomib$gdsc_IC50.sva, label = bortezomib.labels$IC50)
show_pca(input_data = bortezomib$gdsc_IC50.sva, label = sampleinfo.gdsc$tissue.type[bortezomib.labels$IC50_ind])

# AUC
# before sva
# show_pca(input_data = bortezomib$gdsc_AUC, label = bortezomib.labels$AUC)
show_pca(input_data = scale(bortezomib$gdsc_AUC), label = sampleinfo.gdsc$tissue.type[bortezomib.labels$AUC_ind])
# sva
bortezomib$gdsc_AUC.sva <- sva_combine(batch = sampleinfo.gdsc$tissue.type[bortezomib.labels$AUC_ind],  
                                      label = bortezomib.labels$AUC, 
                                      input_data = scale(bortezomib$gdsc_AUC),
                                      n.sv = 3)
mean(bortezomib$gdsc_AUC.sva) # 1.632235e-17
# show_pca(input_data = bortezomib$gdsc_AUC.sva, label = bortezomib.labels$AUC)
show_pca(input_data = bortezomib$gdsc_AUC.sva, label = sampleinfo.gdsc$tissue.type[bortezomib.labels$AUC_ind])

# Get Patient expression data ---------------------------------------------
load("Bortezomib/WS/bortezomib.patient.RData")

bortezomib$patient.combat <- scale(bortezomib.patient_ComBat)
mean(bortezomib$patient.combat) # -5.982353e-19
bortezomib.labels$patient <- binaryResponse == 1
names(bortezomib.labels$patient) <- rownames(bortezomib$patient.combat)
table(bortezomib.labels$patient)
#FALSE  TRUE 
#84    85 

# Using slope labels for SVA  ---------------------------------
temp.data <- comGENE(bortezomib$patient.combat, bortezomib$gdsc_slope.sva)
mean(temp.data[[1]]) #-6.211426e-19
mean(temp.data[[2]]) #5.614e-18
dim(temp.data[[1]]) #169 11791
dim(temp.data[[2]]) #313 11791

# before sva
bortezomib$slope_combined <- rbind(temp.data[[1]], temp.data[[2]])
bortezomib.labels$slope_combined <- c(bortezomib.labels$patient, bortezomib.labels$slope)
bortezomib.labels$slope_combined.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("gdsc", dim(temp.data[[2]])[1]))
#show_pca(input_data = bortezomib$slope_combined, label = bortezomib.labels$slope_combined)
show_pca(input_data = bortezomib$slope_combined, label = bortezomib.labels$slope_combined.source)

# sva
bortezomib$slope_combined.sva <- sva_combine(batch = bortezomib.labels$slope_combined.source,
                                             label = bortezomib.labels$slope_combined,
                                             input_data = bortezomib$slope_combined,
                                             n.sv = 2)
mean(bortezomib$slope_combined.sva) 
# -3.482237e-18

# after sva
#show_pca(input_data = bortezomib$slope_combined.sva, label = bortezomib.labels$slope_combined)
show_pca(input_data = bortezomib$slope_combined.sva, label = bortezomib.labels$slope_combined.source)
rm(temp.data)

# Using IC50 labels for SVA  ---------------------------------
temp.data <- comGENE(bortezomib$patient.combat, bortezomib$gdsc_IC50.sva)
mean(temp.data[[1]]) # -6.211426e-19
mean(temp.data[[2]]) # 2.758975e-18
dim(temp.data[[1]]) # 169 11791
dim(temp.data[[2]]) # 313 11791

# before sva
bortezomib.labels$IC50_combined <- c(bortezomib.labels$patient, bortezomib.labels$IC50)
bortezomib.labels$IC50_combined.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("gdsc", dim(temp.data[[2]])[1]))
bortezomib$IC50_combined <- rbind(temp.data[[1]], temp.data[[2]])
#show_pca(input_data = bortezomib$IC50_combined, label = bortezomib.labels$IC50_combined)
show_pca(input_data = bortezomib$IC50_combined, label = bortezomib.labels$IC50_combined.source)

# sva
bortezomib$IC50_combined.sva <- sva_combine(batch = bortezomib.labels$IC50_combined.source,
                                            label = bortezomib.labels$IC50_combined, 
                                            input_data = bortezomib$IC50_combined, 
                                            n.sv = 2)
mean(bortezomib$IC50_combined.sva) 
# 7.379986e-18

#show_pca(input_data = bortezomib$IC50_combined.sva, label = bortezomib.labels$IC50_combined)
show_pca(input_data = bortezomib$IC50_combined.sva, label = bortezomib.labels$IC50_combined.source)
rm(temp.data)

# Using AUC labels for SVA  ---------------------------------
temp.data <- comGENE(bortezomib$patient.combat, bortezomib$gdsc_AUC.sva)
mean(temp.data[[1]]) #-6.211426e-19
mean(temp.data[[2]]) #1.632534e-17

# before sva
bortezomib.labels$AUC_combined <- c(bortezomib.labels$patient, bortezomib.labels$AUC)
bortezomib.labels$AUC_combined.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("gdsc", dim(temp.data[[2]])[1]))
bortezomib$AUC_combined <- rbind(temp.data[[1]], temp.data[[2]])
#show_pca(input_data = bortezomib$AUC_combined, label = bortezomib.labels$AUC_combined)
show_pca(input_data = bortezomib$AUC_combined, label = bortezomib.labels$AUC_combined.source)

bortezomib$AUC_combined.sva <- sva_combine(batch = bortezomib.labels$AUC_combined.source,
                                           label = bortezomib.labels$AUC_combined, 
                                           input_data = bortezomib$AUC_combined, 
                                           n.sv = 2)
mean(bortezomib$AUC_combined.sva) 
# 3.370679e-18

#show_pca(input_data = bortezomib$AUC_combined.sva, label = bortezomib.labels$AUC_combined)
show_pca(input_data = bortezomib$AUC_combined.sva, label = bortezomib.labels$AUC_combined.source)
rm(temp.data)

# get l1000 features ------------------------------------------------------
stopifnot(colnames(bortezomib$gdsc_IC50.sva) == colnames(bortezomib$gdsc_AUC.sva))
stopifnot(colnames(bortezomib$gdsc_slope.sva) == colnames(bortezomib$gdsc_AUC.sva))

Landmark_Genes_n978 <- read.csv("Common/Landmark_Genes_n978.csv")

feature.l1000 <- list()
feature.l1000$cp <- which(colnames(bortezomib$AUC_combined.sva) %in% Landmark_Genes_n978$Ensembl)
feature.l1000$pp <- which(colnames(bortezomib$patient.combat) %in% Landmark_Genes_n978$Ensembl)

stopifnot(length(feature.l1000$cp) > 0)
stopifnot(length(feature.l1000$pp) > 0)

# create partitions for 10 to 120 patients using all cell lines ------------------------
partition <- list()

input.labels_cell_lines = list()
input.labels_cell_lines$slope = bortezomib.labels$slope
input.labels_cell_lines$IC50 = bortezomib.labels$IC50
input.labels_cell_lines$AUC = bortezomib.labels$AUC

input.cell_line_order = list()
input.cell_line_order$slope = 1:length(input.labels_cell_lines$slope)
input.cell_line_order$IC50 = 1:length(input.labels_cell_lines$IC50)
input.cell_line_order$AUC = 1:length(input.labels_cell_lines$AUC)

stopifnot(input.cell_line_order$slope == input.cell_line_order$IC50)  
stopifnot(input.cell_line_order$AUC == input.cell_line_order$IC50)  
stopifnot(length(input.cell_line_order$AUC) > 0)

cell_lines_all <- foreach(parInd = c(1:11), .errorhandling = "stop") %dopar% {
  
  generate_random_partition(labels.cell_lines = input.labels_cell_lines, 
                            labels.patient = bortezomib.labels$patient,
                            num.training.p = seq(from = 20, to = 120, by = 10)[parInd],
                            cell_line_order = input.cell_line_order,
                            num.training.c = length(input.cell_line_order$AUC),
                            num.test_size = 49
                            )
}

partition$cell_lines_all <- cell_lines_all

# order cell lines by similarity using 100 and 50 patients ------------------------------------------
input.labels_cell_lines = list()
input.labels_cell_lines$slope = bortezomib.labels$slope
input.labels_cell_lines$IC50 = bortezomib.labels$IC50
input.labels_cell_lines$AUC = bortezomib.labels$AUC

input.cell_line_order = list()

input.cell_line_order$slope =
  order_by_similarity(data = bortezomib$slope_combined.sva,
                      labels = bortezomib.labels$slope_combined,
                      source_labels = bortezomib.labels$slope_combined.source)

input.cell_line_order$AUC =
  order_by_similarity(data = bortezomib$AUC_combined.sva,
                      labels = bortezomib.labels$AUC_combined,
                      source_labels = bortezomib.labels$AUC_combined.source)

input.cell_line_order$IC50 =
  order_by_similarity(data = bortezomib$IC50_combined.sva,
                      labels = bortezomib.labels$IC50_combined,
                      source_labels = bortezomib.labels$IC50_combined.source)

stopifnot(length(input.cell_line_order$slope) == length(input.cell_line_order$IC50))
stopifnot(length(input.cell_line_order$AUC) == length(input.cell_line_order$IC50))
stopifnot(length(input.cell_line_order$AUC) > 0)

stopifnot(length(partition$cell_lines_all[[1]]$p2p[[1]]$training_index) == 20)
input.training.c = c(seq(from = 10, to = 100, by = 10), seq(from = 120, to = 313, by = 20))
stopifnot(length(input.training.c) == 20)
patient_20 <- foreach(parInd = c(1:20), .errorhandling = "stop") %dopar% {
  if (parInd == 1) {
    num.min_labels.training = 5
  } else {
    num.min_labels.training = 8
  }
  generate_random_partition(labels.cell_lines = input.labels_cell_lines, 
                            labels.patient = bortezomib.labels$patient,
                            num.training.p = 20,
                            cell_line_order = input.cell_line_order,
                            num.training.c = input.training.c[parInd],
                            input_partition = partition$cell_lines_all[[1]],
                            num.min_labels.training = num.min_labels.training)
}
partition$patient_20 <- patient_20

stopifnot(length(partition$cell_lines_all[[3]]$p2p[[1]]$training_index) == 40)
patient_40 <- foreach(parInd = c(1:20), .errorhandling = "stop") %dopar% {
  if (parInd == 1) {
    num.min_labels.training = 5
  } else {
    num.min_labels.training = 8
  }
  generate_random_partition(labels.cell_lines = input.labels_cell_lines, 
                            labels.patient = bortezomib.labels$patient,
                            num.training.p = 40,
                            cell_line_order = input.cell_line_order,
                            num.training.c = input.training.c[parInd],
                            input_partition = partition$cell_lines_all[[3]],
                            num.min_labels.training = num.min_labels.training)
}
partition$patient_40 <- patient_40

# save worksapce ----------------------------------------------------------
names(bortezomib)
names(bortezomib.labels)
names(partition)
save(sampleinfo.gdsc, bortezomib, bortezomib.labels, feature.l1000, partition,
     file = "Bortezomib/WS/bortezomib_data.RData")
