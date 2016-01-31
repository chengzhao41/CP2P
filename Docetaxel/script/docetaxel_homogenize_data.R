# Load Libraries ----------------------------------------------------------
library("doParallel")
library("sva")

source('Common/preparing_data_helper.R')
source('Common/comGENE.R')
source("Common/generate_random_partition.R")
source("Common/ordering_by_similarity.R")

# Load GDSC Cell line data --------------------------------------------
load("Docetaxel/WS/docetaxel_gdsc.RData")

mean(docetaxel$gdsc_slope) # 6.293791
mean(docetaxel$gdsc_AUC) # 6.293791
mean(docetaxel$gdsc_IC50) # 6.293791
dim(docetaxel$gdsc_slope) #  618 11833
stopifnot(dim(docetaxel$gdsc_slope) == dim(docetaxel$gdsc_AUC))
stopifnot(dim(docetaxel$gdsc_slope) == dim(docetaxel$gdsc_IC50))
table(docetaxel.labels$slope)
# FALSE  TRUE 
# 231   387 
table(docetaxel.labels$AUC)
# FALSE  TRUE 
# 253   365 
table(docetaxel.labels$IC50)
# FALSE  TRUE 
# 219   399

# get patient data --------------------------------------------------------
load("Docetaxel/WS/docetaxel.patient.RData")

docetaxel$patient <- scale(docetaxel.patient)
mean(docetaxel$patient) #-1.359742e-17
dim(docetaxel$patient) #24 8147
docetaxel.labels$patient <- binaryResponse
names(docetaxel.labels$patient) <- rownames(docetaxel.patient)
table(docetaxel.labels$patient)
#FALSE  TRUE 
#14    10 

# Using ComBat to harmonize patients and cell lines - Slopes -----------------
temp.data <- comGENE(docetaxel$patient, scale(docetaxel$gdsc_slope))
mean(temp.data[[1]]) #-1.290866e-17
mean(temp.data[[2]]) #-1.006278e-17
dim(temp.data[[1]]) #24 7969
dim(temp.data[[2]]) #618 7969

# before ComBat
docetaxel.labels$slope_combined <- c(docetaxel.labels$patient, docetaxel.labels$slope)
docetaxel.labels$slope_combined.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("gdsc", dim(temp.data[[2]])[1]))
docetaxel$slope_combined <- rbind(temp.data[[1]], temp.data[[2]])
#show_pca(input_data = docetaxel$slope_combined, label = docetaxel.labels$slope_combined)
show_pca(input_data = docetaxel$slope_combined, label = docetaxel.labels$slope_combined.source)

# ComBat
docetaxel$slope_combined.ComBat = t(ComBat(dat=t(docetaxel$slope_combined), batch = docetaxel.labels$slope_combined.source, 
                                     mod=NULL, par.prior=TRUE, prior.plots=FALSE))
mean(docetaxel$slope_combined.ComBat) #-1.016799e-17

# after ComBat
#show_pca(input_data = docetaxel$slope_combined.ComBat, label = docetaxel.labels$slope_combined)
show_pca(input_data = docetaxel$slope_combined.ComBat, label = docetaxel.labels$slope_combined.source)
rm(temp.data)

# Using ComBat to harmonize patients and cell lines - IC50 -----------------
temp.data <- comGENE(docetaxel$patient, scale(docetaxel$gdsc_IC50))
mean(temp.data[[1]]) #-1.290866e-17
mean(temp.data[[2]]) #-1.006278e-17
dim(temp.data[[1]]) #24 7969
dim(temp.data[[2]]) #618 7969

# before ComBat
docetaxel.labels$IC50_combined <- c(docetaxel.labels$patient, docetaxel.labels$IC50)
docetaxel.labels$IC50_combined.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("gdsc", dim(temp.data[[2]])[1]))
docetaxel$IC50_combined <- rbind(temp.data[[1]], temp.data[[2]])
show_pca(input_data = docetaxel$IC50_combined, label = docetaxel.labels$IC50_combined.source)

# ComBat
docetaxel$IC50_combined.ComBat = t(ComBat(dat=t(docetaxel$IC50_combined), batch = docetaxel.labels$IC50_combined.source, 
                                           mod=NULL, par.prior=TRUE, prior.plots=FALSE))
mean(docetaxel$IC50_combined.ComBat) #-1.016799e-17

# after ComBat
#show_pca(input_data = docetaxel$IC50_combined.ComBat, label = docetaxel.labels$IC50_combined)
show_pca(input_data = docetaxel$IC50_combined.ComBat, label = docetaxel.labels$IC50_combined.source)
rm(temp.data)

# Using sva to harmonize patients and cell lines - AUC -----------------
temp.data <- comGENE(docetaxel$patient, scale(docetaxel$gdsc_AUC))
mean(temp.data[[1]]) #-1.290866e-17
mean(temp.data[[2]]) #-1.006278e-17
dim(temp.data[[1]]) #24 7969
dim(temp.data[[2]]) #618 7969

# before ComBat
docetaxel.labels$AUC_combined <- c(docetaxel.labels$patient, docetaxel.labels$AUC)
docetaxel.labels$AUC_combined.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("gdsc", dim(temp.data[[2]])[1]))
docetaxel$AUC_combined <- rbind(temp.data[[1]], temp.data[[2]])
show_pca(input_data = docetaxel$AUC_combined, label = docetaxel.labels$AUC_combined.source)

# ComBat
docetaxel$AUC_combined.ComBat = t(ComBat(dat=t(docetaxel$AUC_combined), batch = docetaxel.labels$AUC_combined.source, 
                                          mod=NULL, par.prior=TRUE, prior.plots=FALSE))
mean(docetaxel$AUC_combined.ComBat) #-1.016799e-17

# after ComBat
#show_pca(input_data = docetaxel$AUC_combined.ComBat, label = docetaxel.labels$AUC_combined)
show_pca(input_data = docetaxel$AUC_combined.ComBat, label = docetaxel.labels$AUC_combined.source)
rm(temp.data)

# Breast Cell lines only --------------------------------------------------
## IC50
temp.ind <- which(sampleinfo.gdsc$tissue.type[docetaxel.labels$IC50_ind] == "breast")
temp.data <- comGENE(docetaxel$patient, scale(docetaxel$gdsc_IC50[temp.ind, ]))
mean(temp.data[[1]]) #-1.290866e-17
mean(temp.data[[2]]) #-9.652955e-18
dim(temp.data[[2]]) #35 7969

# before ComBat
docetaxel.labels$IC50_breast.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("gdsc", dim(temp.data[[2]])[1]))
docetaxel.labels$IC50_breast_only <- docetaxel.labels$IC50[temp.ind]
docetaxel.labels$IC50_breast <- c(docetaxel.labels$patient, docetaxel.labels$IC50[temp.ind])
docetaxel$IC50_breast <- rbind(temp.data[[1]], temp.data[[2]])
show_pca(input_data = docetaxel$IC50_breast, label = docetaxel.labels$IC50_breast.source)

# ComBat
docetaxel$IC50_breast.ComBat <- t(ComBat(dat=t(docetaxel$IC50_breast), batch = docetaxel.labels$IC50_breast.source, 
         mod=NULL, par.prior=TRUE, prior.plots=FALSE))
mean(docetaxel$IC50_breast.ComBat) #-1.091286e-17

# after ComBat
#show_pca(input_data = docetaxel$IC50_breast.ComBat, label = docetaxel.labels$IC50_breast)
show_pca(input_data = docetaxel$IC50_breast.ComBat, label = docetaxel.labels$IC50_breast.source)
rm(temp.data)

# AUC
temp.ind <- which(sampleinfo.gdsc$tissue.type[docetaxel.labels$AUC_ind] == "breast")
temp.data <- comGENE(docetaxel$patient, scale(docetaxel$gdsc_AUC[temp.ind, ]))
mean(temp.data[[1]]) #-1.290866e-17
mean(temp.data[[2]]) #-9.652955e-18
dim(temp.data[[2]]) #35 7969

# before ComBat
docetaxel.labels$AUC_breast.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("gdsc", dim(temp.data[[2]])[1]))
docetaxel.labels$AUC_breast_only <- docetaxel.labels$AUC[temp.ind]
docetaxel.labels$AUC_breast <- c(docetaxel.labels$patient, docetaxel.labels$AUC[temp.ind])
docetaxel$AUC_breast <- rbind(temp.data[[1]], temp.data[[2]])
show_pca(input_data = docetaxel$AUC_breast, label = docetaxel.labels$AUC_breast.source)

# ComBat
docetaxel$AUC_breast.ComBat <- t(ComBat(dat=t(docetaxel$AUC_breast), batch = docetaxel.labels$AUC_breast.source, 
                                         mod=NULL, par.prior=TRUE, prior.plots=FALSE))
mean(docetaxel$AUC_breast.ComBat) #-1.091286e-17

# after ComBat
#show_pca(input_data = docetaxel$AUC_breast.ComBat, label = docetaxel.labels$AUC_breast)
show_pca(input_data = docetaxel$AUC_breast.ComBat, label = docetaxel.labels$AUC_breast.source)
rm(temp.data)

# Slope
temp.ind <- which(sampleinfo.gdsc$tissue.type[docetaxel.labels$slope_ind] == "breast")
temp.data <- comGENE(docetaxel$patient, scale(docetaxel$gdsc_slope[temp.ind, ]))
mean(temp.data[[1]]) #-1.290866e-17
mean(temp.data[[2]]) #-9.652955e-18
dim(temp.data[[2]]) #35 7969

# before ComBat
docetaxel.labels$slope_breast.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("gdsc", dim(temp.data[[2]])[1]))
docetaxel.labels$slope_breast_only <- docetaxel.labels$slope[temp.ind]
docetaxel.labels$slope_breast <- c(docetaxel.labels$patient, docetaxel.labels$slope[temp.ind])
docetaxel$slope_breast <- rbind(temp.data[[1]], temp.data[[2]])
show_pca(input_data = docetaxel$slope_breast, label = docetaxel.labels$slope_breast.source)

# ComBat
docetaxel$slope_breast.ComBat <- t(ComBat(dat=t(docetaxel$slope_breast), batch = docetaxel.labels$slope_breast.source, 
                                        mod=NULL, par.prior=TRUE, prior.plots=FALSE))
mean(docetaxel$slope_breast.ComBat) #-1.091286e-17

# after ComBat
#show_pca(input_data = docetaxel$slope_breast.ComBat, label = docetaxel.labels$slope_breast)
show_pca(input_data = docetaxel$slope_breast.ComBat, label = docetaxel.labels$slope_breast.source)
rm(temp.data)

# get l1000 features ------------------------------------------------------
stopifnot(colnames(docetaxel$gdsc_IC50) == colnames(docetaxel$gdsc_AUC))
stopifnot(colnames(docetaxel$gdsc_slope) == colnames(docetaxel$gdsc_AUC))

Landmark_Genes_n978 <- read.csv("Common/Landmark_Genes_n978.csv")
feature.l1000 <- list()
feature.l1000$cp <- which(colnames(docetaxel$slope_combined.ComBat) %in% Landmark_Genes_n978$Ensembl)
feature.l1000$pp <- which(colnames(docetaxel$patient) %in% Landmark_Genes_n978$Ensembl)

stopifnot(length(feature.l1000$cp) > 0)
stopifnot(length(feature.l1000$pp) > 0)

# create partitions -------------------------------------------------------
partition <- list()

input.labels_cell_lines = list()
input.labels_cell_lines$slope = docetaxel.labels$slope
input.labels_cell_lines$IC50 = docetaxel.labels$IC50
input.labels_cell_lines$AUC = docetaxel.labels$AUC

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
input.labels_cell_lines$slope = docetaxel.labels$slope_breast
input.labels_cell_lines$IC50 = docetaxel.labels$IC50_breast
input.labels_cell_lines$AUC = docetaxel.labels$AUC_breast

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
input.labels_cell_lines$slope = docetaxel.labels$slope
input.labels_cell_lines$IC50 = docetaxel.labels$IC50
input.labels_cell_lines$AUC = docetaxel.labels$AUC

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
warnings() # should only warn about "using leave-one"
partition$patient_23 <- patient_23

# save worksapce ----------------------------------------------------------
save(docetaxel, docetaxel.labels, sampleinfo.gdsc, partition, feature.l1000, partition, 
     file = "Docetaxel/WS/docetaxel_data.RData")
