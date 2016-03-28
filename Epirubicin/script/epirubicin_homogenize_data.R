# Load Libraries ----------------------------------------------------------
library("doParallel")
library("sva")

source('Common/preparing_data_helper.R')
source('Common/comGENE.R')
source("Common/generate_random_partition.R")
source("Common/ordering_by_similarity.R")

# Load CGP Cell line data --------------------------------------------
load("Epirubicin/WS/epirubicin_gray.RData")

# sanity check to make sure we have the right version of data
mean(epirubicin$gray_slope) # 0.9209632
mean(epirubicin$gray_AUC) # 0.9209632
mean(epirubicin$gray_IC50) # 0.9209632
dim(epirubicin$gray_slope) # 313 61958
dim(epirubicin$gray_AUC) # 313 61958
dim(epirubicin$gray_IC50) # 313 61958

# remove genes with 0 variance
temp.var <- apply(epirubicin$gray_slope, 2, var)
epirubicin$gray_slope <- epirubicin$gray_slope[, -which(temp.var == 0)]
dim(epirubicin$gray_slope)
#38 45940
temp.var <- apply(epirubicin$gray_AUC, 2, var)
epirubicin$gray_AUC <- epirubicin$gray_AUC[, -which(temp.var == 0)]
dim(epirubicin$gray_AUC)
#38 45940
temp.var <- apply(epirubicin$gray_IC50, 2, var)
epirubicin$gray_IC50 <- epirubicin$gray_IC50[, -which(temp.var == 0)]
dim(epirubicin$gray_IC50)
#38 45940

table(epirubicin.labels$slope)
# FALSE  TRUE 
# 5    33 
table(epirubicin.labels$AUC)
# FALSE  TRUE 
# 7    31 
table(epirubicin.labels$IC50)
# FALSE  TRUE 
# 2    36 

# Get Patient expression data ---------------------------------------------
load("Epirubicin/WS/epirubicin.patient.RData")

epirubicin$patient <- scale(t(epirubicin.patient))
mean(epirubicin$patient) # 3.138273e-18
epirubicin.labels$patient <- binaryResponse
names(epirubicin.labels$patient) <- rownames(epirubicin$patient)
table(epirubicin.labels$patient)
# FALSE  TRUE 
# 101    17 

# Using slope labels for SVA  ---------------------------------
temp.data <- comGENE(epirubicin$patient, scale(epirubicin$gray_slope))
mean(temp.data[[1]]) #1.191286e-18
mean(temp.data[[2]]) #1.56495e-19
dim(temp.data[[1]]) #118 19100
dim(temp.data[[2]]) #38 19343

# before sva
epirubicin$slope_combined <- rbind(temp.data[[1]], temp.data[[2]])
epirubicin.labels$slope_combined <- c(epirubicin.labels$patient, epirubicin.labels$slope)
epirubicin.labels$slope_combined.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("gray", dim(temp.data[[2]])[1]))
#show_pca(input_data = epirubicin$slope_combined, label = epirubicin.labels$slope_combined)
show_pca(input_data = epirubicin$slope_combined, label = epirubicin.labels$slope_combined.source)
# not the outlier greatly skews the red eclipse 

# sva
epirubicin$slope_combined.sva <- sva_combine(batch = epirubicin.labels$slope_combined.source,
                                             label = epirubicin.labels$slope_combined,
                                             input_data = epirubicin$slope_combined,
                                             n.sv = 2)
mean(epirubicin$slope_combined.sva) 
# 1.994808e-17

# after sva
#show_pca(input_data = epirubicin$slope_combined.sva, label = epirubicin.labels$slope_combined)
show_pca(input_data = epirubicin$slope_combined.sva, label = epirubicin.labels$slope_combined.source)
rm(temp.data)

# Using IC50 labels for SVA  ---------------------------------
temp.data <- comGENE(epirubicin$patient, scale(epirubicin$gray_IC50))
mean(temp.data[[1]]) # 1.191286e-18
mean(temp.data[[2]]) # 1.56495e-19
dim(temp.data[[1]]) # 118 19100
dim(temp.data[[2]]) # 38 19100

# before sva
epirubicin.labels$IC50_combined <- c(epirubicin.labels$patient, epirubicin.labels$IC50)
epirubicin.labels$IC50_combined.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("gray", dim(temp.data[[2]])[1]))
epirubicin$IC50_combined <- rbind(temp.data[[1]], temp.data[[2]])
#show_pca(input_data = epirubicin$IC50_combined, label = epirubicin.labels$IC50_combined)
show_pca(input_data = epirubicin$IC50_combined, label = epirubicin.labels$IC50_combined.source)

# sva
epirubicin$IC50_combined.sva <- sva_combine(batch = epirubicin.labels$IC50_combined.source,
                                            label = epirubicin.labels$IC50_combined, 
                                            input_data = epirubicin$IC50_combined, 
                                            n.sv = 2)
mean(epirubicin$IC50_combined.sva) 
# -1.031379e-17

#show_pca(input_data = epirubicin$IC50_combined.sva, label = epirubicin.labels$IC50_combined)
show_pca(input_data = epirubicin$IC50_combined.sva, label = epirubicin.labels$IC50_combined.source)
rm(temp.data)

# Using AUC labels for SVA  ---------------------------------
temp.data <- comGENE(epirubicin$patient, scale(epirubicin$gray_AUC))
mean(temp.data[[1]]) #1.191286e-18
mean(temp.data[[2]]) #1.56495e-19

# before sva
epirubicin.labels$AUC_combined <- c(epirubicin.labels$patient, epirubicin.labels$AUC)
epirubicin.labels$AUC_combined.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("gray", dim(temp.data[[2]])[1]))
epirubicin$AUC_combined <- rbind(temp.data[[1]], temp.data[[2]])
#show_pca(input_data = epirubicin$AUC_combined, label = epirubicin.labels$AUC_combined)
show_pca(input_data = epirubicin$AUC_combined, label = epirubicin.labels$AUC_combined.source)

epirubicin$AUC_combined.sva <- sva_combine(batch = epirubicin.labels$AUC_combined.source,
                                           label = epirubicin.labels$AUC_combined, 
                                           input_data = epirubicin$AUC_combined, 
                                           n.sv = 2)
mean(epirubicin$AUC_combined.sva) 
# 3.530146e-17

#show_pca(input_data = epirubicin$AUC_combined.sva, label = epirubicin.labels$AUC_combined)
show_pca(input_data = epirubicin$AUC_combined.sva, label = epirubicin.labels$AUC_combined.source)
rm(temp.data)

# get l1000 features ------------------------------------------------------
stopifnot(colnames(epirubicin$gray_IC50.sva) == colnames(epirubicin$gray_AUC.sva))
stopifnot(colnames(epirubicin$gray_slope.sva) == colnames(epirubicin$gray_AUC.sva))

Landmark_Genes_n978 <- read.csv("Common/Landmark_Genes_n978.csv")

feature.l1000 <- list()
feature.l1000$cp <- which(colnames(epirubicin$AUC_combined.sva) %in% Landmark_Genes_n978$Ensembl)
feature.l1000$pp <- which(colnames(epirubicin$patient) %in% Landmark_Genes_n978$Ensembl)

stopifnot(length(feature.l1000$cp) > 0)
stopifnot(length(feature.l1000$pp) > 0)

# create partitions for 10 to 150 patients using all cell lines ------------------------
partition <- list()

input.labels_cell_lines = list()
input.labels_cell_lines$slope = epirubicin.labels$slope
input.labels_cell_lines$IC50 = epirubicin.labels$IC50
input.labels_cell_lines$AUC = epirubicin.labels$AUC

input.cell_line_order = list()
input.cell_line_order$slope = 1:length(input.labels_cell_lines$slope)
input.cell_line_order$IC50 = 1:length(input.labels_cell_lines$IC50)
input.cell_line_order$AUC = 1:length(input.labels_cell_lines$AUC)

stopifnot(input.cell_line_order$slope == input.cell_line_order$IC50)  
stopifnot(input.cell_line_order$AUC == input.cell_line_order$IC50)  
stopifnot(length(input.cell_line_order$AUC) > 0)

cell_lines_all <- foreach(input.training_amount.p = seq(from = 20, to = 100, by = 10)
                          , .errorhandling = "stop") %dopar% {
                            
                            generate_random_partition(labels_cell_lines = input.labels_cell_lines, 
                                                      labels_patient = epirubicin.labels$patient,
                                                      training_amount.p = input.training_amount.p,
                                                      leave_one_out = FALSE,
                                                      cell_line_order = input.cell_line_order,
                                                      training_amount.c = length(input.cell_line_order$AUC),
                                                      acc_training = FALSE
                            )
                          }

partition$cell_lines_all <- cell_lines_all

# order cell lines by similarity using 80 and 40 patients ------------------------------------------
load("CGP/cosmic.tcga.RData")
input.labels_cell_lines = list()
input.labels_cell_lines$slope = epirubicin.labels$slope
input.labels_cell_lines$IC50 = epirubicin.labels$IC50
input.labels_cell_lines$AUC = epirubicin.labels$AUC

input.cell_line_order = list()
input.cell_line_order$slope <- order(match(names(input.labels_cell_lines$slope), brca_ordered$cell.line))
input.cell_line_order$IC50 <- order(match(names(input.labels_cell_lines$IC50), brca_ordered$cell.line))
input.cell_line_order$AUC <- order(match(names(input.labels_cell_lines$AUC), brca_ordered$cell.line))

stopifnot(!is.unsorted(match(brca_ordered$cell.line, 
                             names(input.labels_cell_lines$slope)[input.cell_line_order$slope] ), na.rm = TRUE))
stopifnot( length(table(input.labels_cell_lines$slope[input.cell_line_order$slope][1:30])) == 2 )
stopifnot( length(table(input.labels_cell_lines$AUC[input.cell_line_order$AUC][1:30])) == 2 )
stopifnot( length(table(input.labels_cell_lines$IC50[input.cell_line_order$IC50][1:30])) == 2 )

stopifnot(length(input.cell_line_order$slope) == length(input.cell_line_order$IC50))
stopifnot(length(input.cell_line_order$AUC) == length(input.cell_line_order$IC50))
stopifnot(length(input.cell_line_order$AUC) > 0)

stopifnot(length(partition$cell_lines_all[[7]]$p2p[[1]]$training_index) == 80)
patient_80 <- foreach(input.training_amount.c = seq(from = 5, to = 38, by = 5), .errorhandling = "stop") %dopar% {
  generate_random_partition(labels_cell_lines = input.labels_cell_lines, 
                            labels_patient = epirubicin.labels$patient,
                            training_amount.p = 80,
                            leave_one_out = FALSE,
                            cell_line_order = input.cell_line_order,
                            training_amount.c = input.training_amount.c,
                            input_partition = partition$cell_lines_all[[7]],
                            acc_training = FALSE)
}
partition$patient_80 <- patient_80

stopifnot(length(partition$cell_lines_all[[3]]$p2p[[3]]$training_index) == 40)
patient_40 <- foreach(input.training_amount.c = seq(from = 5, to = 38, by = 5), .errorhandling = "stop") %dopar% {
  generate_random_partition(labels_cell_lines = input.labels_cell_lines, 
                            labels_patient = epirubicin.labels$patient,
                            training_amount.p = 40,
                            leave_one_out = FALSE,
                            cell_line_order = input.cell_line_order,
                            training_amount.c = input.training_amount.c,
                            input_partition = partition$cell_lines_all[[3]],
                            acc_training = FALSE)
}
partition$patient_40 <- patient_40

# save worksapce ----------------------------------------------------------
save(sampleinfo.gray, epirubicin, epirubicin.labels, feature.l1000, partition,
     file = "Epirubicin/WS/epirubicin_data.RData")
