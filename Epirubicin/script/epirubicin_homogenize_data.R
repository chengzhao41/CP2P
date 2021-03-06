# Load Libraries ----------------------------------------------------------
library("doParallel")

source('Common/preparing_data_helper.R')
source('Common/comGENE.R')
source("Common/generate_random_partition.R")

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

# not using IC50 labels, because it is too imbalanced
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

# Using slope labels ---------------------------------
temp.data <- comGENE(epirubicin$patient, scale(epirubicin$gray_slope))
mean(temp.data[[1]]) #1.191286e-18
mean(temp.data[[2]]) #1.56495e-19
dim(temp.data[[1]]) #118 19100
dim(temp.data[[2]]) #38 19343

epirubicin$slope_combined <- rbind(temp.data[[1]], temp.data[[2]])
epirubicin.labels$slope_combined <- c(epirubicin.labels$patient, epirubicin.labels$slope)
epirubicin.labels$slope_combined.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("gray", dim(temp.data[[2]])[1]))
# show_pca(input_data = epirubicin$slope_combined, label = epirubicin.labels$slope_combined)
# show_pca(input_data = epirubicin$slope_combined, label = epirubicin.labels$slope_combined.source)

# Using AUC labels ---------------------------------
rm(temp.data)
temp.data <- comGENE(epirubicin$patient, scale(epirubicin$gray_AUC))
mean(temp.data[[1]]) #1.191286e-18
mean(temp.data[[2]]) #1.56495e-19

epirubicin.labels$AUC_combined <- c(epirubicin.labels$patient, epirubicin.labels$AUC)
epirubicin.labels$AUC_combined.source <- c(rep("patient", dim(temp.data[[1]])[1]), rep("gray", dim(temp.data[[2]])[1]))
epirubicin$AUC_combined <- rbind(temp.data[[1]], temp.data[[2]])
# show_pca(input_data = epirubicin$AUC_combined, label = epirubicin.labels$AUC_combined)
# show_pca(input_data = epirubicin$AUC_combined, label = epirubicin.labels$AUC_combined.source)

# get l1000 features ------------------------------------------------------
stopifnot(colnames(epirubicin$gray_IC50) == colnames(epirubicin$gray_AUC))
stopifnot(colnames(epirubicin$gray_slope) == colnames(epirubicin$gray_AUC))

Landmark_Genes_n978 <- read.csv("Common/Landmark_Genes_n978.csv")

feature.l1000 <- list()
feature.l1000$cp <- which(colnames(epirubicin$AUC_combined) %in% Landmark_Genes_n978$Ensembl)
feature.l1000$pp <- which(colnames(epirubicin$patient) %in% Landmark_Genes_n978$Ensembl)

stopifnot(length(feature.l1000$cp) > 0)
stopifnot(length(feature.l1000$pp) > 0)

# create partitions for 10 to 150 patients using all cell lines ------------------------
partition <- list()

input.labels_cell_lines = list()
input.labels_cell_lines$slope = epirubicin.labels$slope
input.labels_cell_lines$AUC = epirubicin.labels$AUC

input.cell_line_order = list()
input.cell_line_order$slope = 1:length(input.labels_cell_lines$slope)
input.cell_line_order$AUC = 1:length(input.labels_cell_lines$AUC)

stopifnot(input.cell_line_order$slope == input.cell_line_order$AUC)  
stopifnot(length(input.cell_line_order$AUC) > 0)

cell_lines_all <- foreach(parInd = 1:11, .errorhandling = "stop") %dopar% {
                            
    generate_random_partition(labels.cell_lines = input.labels_cell_lines, 
                              labels.patient = epirubicin.labels$patient,
                              num.training.p = seq(from = 20, to = 70, by = 5)[parInd],
                              cell_line_order = input.cell_line_order,
                              num.training.c = length(input.cell_line_order$AUC),
                              num.test_size = 48,
                              num.min_labels.training.c = 5
                              )
}

partition$cell_lines_all <- cell_lines_all

# order cell lines by brca mutations ------------------------------------------
load("CGP/cosmic.tcga.RData")
input.labels_cell_lines = list()
input.labels_cell_lines$slope = epirubicin.labels$slope
input.labels_cell_lines$AUC = epirubicin.labels$AUC

input.cell_line_order = list()
input.cell_line_order$slope <- order(match(names(input.labels_cell_lines$slope), brca_ordered$cell.line))
input.cell_line_order$AUC <- order(match(names(input.labels_cell_lines$AUC), brca_ordered$cell.line))

stopifnot(!is.unsorted(match(brca_ordered$cell.line, 
                             names(input.labels_cell_lines$slope)[input.cell_line_order$slope] ), na.rm = TRUE))
stopifnot( length(table(input.labels_cell_lines$slope[input.cell_line_order$slope][1:30])) == 2 )
stopifnot( length(table(input.labels_cell_lines$AUC[input.cell_line_order$AUC][1:30])) == 2 )

stopifnot(length(input.cell_line_order$slope) == length(input.cell_line_order$AUC))
stopifnot(length(input.cell_line_order$AUC) > 0)

input.training.c = seq(from = 10, to = 38, by = 5)
patient_20 <- foreach(parInd = 1:6, .errorhandling = "stop") %dopar% {
  
  stopifnot(length(partition$cell_lines_all[[parInd + 2]]$p2p[[1]]$training_index) == 20 + input.training.c[parInd])
  
  generate_random_partition(labels.cell_lines = input.labels_cell_lines, 
                            labels.patient = epirubicin.labels$patient,
                            num.training.p = 20,
                            cell_line_order = input.cell_line_order,
                            num.training.c = input.training.c[parInd],
                            input_partition = partition$cell_lines_all[[parInd + 2]],
                            num.test_size = 48,
                            num.min_labels.training.c = 5)
}
partition$patient_20 <- patient_20


# save worksapce ----------------------------------------------------------
save(sampleinfo.gray, epirubicin, epirubicin.labels, feature.l1000, partition,
     file = "Epirubicin/WS/epirubicin_data.RData")
