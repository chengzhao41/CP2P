# Load Libraries ----------------------------------------------------------
library("doParallel")
library("sva")

source('Common/preparing_data_helper.R')
source('Common/comGENE.R')
source("Common/generate_random_partition.R")
source("Common/ordering_by_similarity.R")

# Load CGP Cell line data --------------------------------------------
load("Erlotinib/WS/erlotinib_ccle.RData")

# sanity check to make sure we have the right version of data
mean(erlotinib$ccle_slope) # 6.045301
mean(erlotinib$ccle_AUC) # 6.045301
mean(erlotinib$ccle_IC50) # 6.045301
dim(erlotinib$ccle_slope) # 489 20049
dim(erlotinib$ccle_AUC) # 489 20049
dim(erlotinib$ccle_IC50) # 489 20049
table(erlotinib.labels$slope)
# FALSE  TRUE 
# 396    93 
table(erlotinib.labels$AUC)
# FALSE  TRUE 
# 385   104 
table(erlotinib.labels$IC50)
# FALSE  TRUE 
# 402    87 

# Get Patient expression data ---------------------------------------------
load("Erlotinib/WS/erlotinib.patient.RData")

erlotinib$patient <- scale(erlotinib.patient)
mean(erlotinib$patient) # -1.622537e-17
erlotinib.labels$patient <- binaryResponse
stopifnot(names(erlotinib.labels$patient) == rownames(erlotinib$patient))
table(erlotinib.labels$patient)
# FALSE  TRUE 
# 14    11 

# Get the Battle expression data ------------------------------------------
load("Erlotinib/WS/erlotinib_battle.RData")

erlotinib$battle <- scale(battle$IC50)
mean(erlotinib$battle) # -3.273899e-18
erlotinib.labels$battle_IC50 <- battle.labels$IC50
stopifnot(names(erlotinib.labels$battle_IC50) == rownames(erlotinib$battle))
table(erlotinib.labels$battle_IC50)
# FALSE  TRUE 
# 20    24

# Finding common set of genes  ---------------------------------
temp.data1 <- comGENE(erlotinib$patient, erlotinib$battle)
mean(temp.data1[[1]]) #-1.85417e-17
mean(temp.data1[[2]]) #-2.855084e-18
dim(temp.data1[[1]]) #25 16502
dim(temp.data1[[2]]) #44 16502

temp.data2 <- comGENE(temp.data1[[1]], scale(erlotinib$ccle_slope))
mean(temp.data2[[1]]) #-2.189319e-17
mean(temp.data2[[2]]) #6.591745e-18
dim(temp.data2[[1]]) #25 15237
dim(temp.data2[[2]]) #489 15237

temp.data3 <- comGENE(temp.data1[[2]], temp.data2[[2]])
mean(temp.data3[[1]]) #-6.766092e-18
mean(temp.data3[[2]]) #6.591745e-18
dim(temp.data3[[1]]) #44 15237
dim(temp.data3[[2]]) #489 15237

temp.patient <- temp.data2[[1]]
erlotinib$battle <- temp.data3[[1]]
erlotinib$ccle_slope <- temp.data3[[2]]
stopifnot(colnames(temp.patient) == colnames(erlotinib$battle))
stopifnot(colnames(erlotinib$ccle_slope) == colnames(erlotinib$battle))

# remove intersecting cell lines, use CCLE
temp.dup_ind <- which(names(erlotinib.labels$battle_IC50) %in% names(erlotinib.labels$slope))
erlotinib$battle <- erlotinib$battle[-temp.dup_ind, ]
erlotinib.labels$battle_IC50 <- erlotinib.labels$battle_IC50[-temp.dup_ind]

# Homogenizing with ccle slope labels using ComBat  ---------------------------------
mean(temp.patient) #-2.189319e-17
mean(erlotinib$battle) #-6.766092e-18
mean(erlotinib$ccle_slope) # 6.591745e-18

# before
erlotinib$slope_combined <- rbind(temp.patient, erlotinib$battle, erlotinib$ccle_slope)
erlotinib.labels$slope_combined <- c(erlotinib.labels$patient, erlotinib.labels$battle_IC50, erlotinib.labels$slope)
erlotinib.labels$slope_combined.source <- c(rep("patient", dim(temp.patient)[1]), rep("battle", dim(erlotinib$battle)[1]), rep("ccle", dim(erlotinib$ccle_slope)[1]))

#show_pca(input_data = erlotinib$slope_combined, label = erlotinib.labels$slope_combined)
show_pca(input_data = erlotinib$slope_combined, label = erlotinib.labels$slope_combined.source)

# ComBat
erlotinib$slope_combined.ComBat <- ComBat_combine(batch = erlotinib.labels$slope_combined.source,
                                                  label = erlotinib.labels$slope_combined,
                                                  input_data = erlotinib$slope_combined)
mean(erlotinib$slope_combined.ComBat) 
# 0.0001254162

# after
#show_pca(input_data = erlotinib$slope_combined.ComBat, label = erlotinib.labels$slope_combined)
show_pca(input_data = erlotinib$slope_combined.ComBat, label = erlotinib.labels$slope_combined.source)

# Homogenizing with ccle IC50 labels using ComBat  ---------------------------------
temp.data <- comGENE(temp.patient, scale(erlotinib$ccle_IC50))
mean(temp.data[[1]]) # -2.189319e-17
mean(temp.data[[2]]) # 6.591745e-18
dim(temp.data[[1]]) # 25 15237
dim(temp.data[[2]]) # 489 15237

erlotinib$ccle_IC50 <- temp.data[[2]]

# before
erlotinib$IC50_combined <- rbind(temp.patient, erlotinib$battle, erlotinib$ccle_IC50)
erlotinib.labels$IC50_combined <- c(erlotinib.labels$patient, erlotinib.labels$battle_IC50, erlotinib.labels$IC50)
erlotinib.labels$IC50_combined.source <- c(rep("patient", dim(temp.patient)[1]), 
                                           rep("battle", dim(erlotinib$battle)[1]), 
                                           rep("ccle", dim(erlotinib$ccle_IC50)[1]))

#show_pca(input_data = erlotinib$IC50_combined, label = erlotinib.labels$IC50_combined)
show_pca(input_data = erlotinib$IC50_combined, label = erlotinib.labels$IC50_combined.source)

# ComBat
erlotinib$IC50_combined.ComBat <- ComBat_combine(batch = erlotinib.labels$IC50_combined.source,
                                                 label = erlotinib.labels$IC50_combined,
                                                 input_data = erlotinib$IC50_combined)
mean(erlotinib$IC50_combined.ComBat) 
# 0.000121424

#show_pca(input_data = erlotinib$IC50_combined.ComBat, label = erlotinib.labels$IC50_combined)
show_pca(input_data = erlotinib$IC50_combined.ComBat, label = erlotinib.labels$IC50_combined.source)

# Homogenizing with ccle AUC labels using ComBat  ---------------------------------
temp.data <- comGENE(temp.patient, scale(erlotinib$ccle_AUC))
mean(temp.data[[1]]) # -2.189319e-17
mean(temp.data[[2]]) # 6.591745e-18
dim(temp.data[[1]]) # 25 15237
dim(temp.data[[2]]) # 489 15237

erlotinib$ccle_AUC <- temp.data[[2]]

# before
erlotinib$AUC_combined <- rbind(temp.patient, erlotinib$battle, erlotinib$ccle_AUC)
erlotinib.labels$AUC_combined <- c(erlotinib.labels$patient, erlotinib.labels$battle_IC50, erlotinib.labels$AUC)
erlotinib.labels$AUC_combined.source <- c(rep("patient", dim(temp.patient)[1]), 
                                          rep("battle", dim(erlotinib$battle)[1]), 
                                          rep("ccle", dim(erlotinib$ccle_AUC)[1]))

#show_pca(input_data = erlotinib$AUC_combined, label = erlotinib.labels$AUC_combined)
show_pca(input_data = erlotinib$AUC_combined, label = erlotinib.labels$AUC_combined.source)

# ComBat
erlotinib$AUC_combined.ComBat <- ComBat_combine(batch = erlotinib.labels$AUC_combined.source,
                                                label = erlotinib.labels$AUC_combined,
                                                input_data = erlotinib$AUC_combined)
mean(erlotinib$AUC_combined.ComBat) 
# 0.0001268236

#show_pca(input_data = erlotinib$IC50_combined.ComBat, label = erlotinib.labels$IC50_combined)
show_pca(input_data = erlotinib$IC50_combined.ComBat, label = erlotinib.labels$IC50_combined.source)


# homogenizing with only lung cell lines ------------------------------------------------------
# Homogenizing with ccle slope labels using ComBat  ---------------------------------
mean(temp.patient) #-2.189319e-17
mean(erlotinib$battle) #0.002000544
mean(erlotinib$ccle_slope) #6.591745e-18

# before
temp.ind <- which(sampleinfo.ccle$tissueid[erlotinib.labels$slope_ind] == "lung")
erlotinib$slope_lung <- rbind(temp.patient, erlotinib$battle, erlotinib$ccle_slope[temp.ind, ])
erlotinib.labels$slope_lung <- c(erlotinib.labels$patient, erlotinib.labels$battle_IC50, erlotinib.labels$slope[temp.ind])
erlotinib.labels$slope_lung.source <- c(rep("patient", dim(temp.patient)[1]), rep("battle", dim(erlotinib$battle)[1]), rep("ccle", dim(erlotinib$ccle_slope[temp.ind, ])[1]))

#show_pca(input_data = erlotinib$slope_lung, label = erlotinib.labels$slope_lung)
show_pca(input_data = erlotinib$slope_lung, label = erlotinib.labels$slope_lung.source)

# ComBat
erlotinib$slope_lung.ComBat <- ComBat_combine(batch = erlotinib.labels$slope_lung.source,
                                              label = erlotinib.labels$slope_lung,
                                              input_data = erlotinib$slope_lung)
mean(erlotinib$slope_lung.ComBat) 
# -0.006142389

# after
#show_pca(input_data = erlotinib$slope_lung.ComBat, label = erlotinib.labels$slope_lung)
show_pca(input_data = erlotinib$slope_lung.ComBat, label = erlotinib.labels$slope_lung.source)

# Homogenizing with ccle IC50 labels using ComBat  ---------------------------------
# before
temp.ind <- which(sampleinfo.ccle$tissueid[erlotinib.labels$IC50_ind] == "lung")
erlotinib$IC50_lung <- rbind(temp.patient, erlotinib$battle, erlotinib$ccle_IC50[temp.ind, ])
erlotinib.labels$IC50_lung <- c(erlotinib.labels$patient, erlotinib.labels$battle_IC50, erlotinib.labels$IC50[temp.ind])
erlotinib.labels$IC50_lung.source <- c(rep("patient", dim(temp.patient)[1]), 
                                       rep("battle", dim(erlotinib$battle)[1]), 
                                       rep("ccle", dim(erlotinib$ccle_IC50[temp.ind, ])[1]))

#show_pca(input_data = erlotinib$IC50_lung, label = erlotinib.labels$IC50_lung)
show_pca(input_data = erlotinib$IC50_lung, label = erlotinib.labels$IC50_lung.source)

# ComBat
erlotinib$IC50_lung.ComBat <- ComBat_combine(batch = erlotinib.labels$IC50_lung.source,
                                             label = erlotinib.labels$IC50_lung,
                                             input_data = erlotinib$IC50_lung)
mean(erlotinib$IC50_lung.ComBat) 
# -0.006148108

#show_pca(input_data = erlotinib$IC50_lung.ComBat, label = erlotinib.labels$IC50_lung)
show_pca(input_data = erlotinib$IC50_lung.ComBat, label = erlotinib.labels$IC50_lung.source)

# Homogenizing with ccle AUC labels using ComBat  ---------------------------------
# before
temp.ind <- which(sampleinfo.ccle$tissueid[erlotinib.labels$AUC_ind] == "lung")
erlotinib$AUC_lung <- rbind(temp.patient, erlotinib$battle, erlotinib$ccle_AUC[temp.ind, ])
erlotinib.labels$AUC_lung <- c(erlotinib.labels$patient, erlotinib.labels$battle_IC50, erlotinib.labels$AUC[temp.ind])
erlotinib.labels$AUC_lung.source <- c(rep("patient", dim(temp.patient)[1]), 
                                      rep("battle", dim(erlotinib$battle)[1]), 
                                      rep("ccle", dim(erlotinib$ccle_AUC[temp.ind, ])[1]))

#show_pca(input_data = erlotinib$AUC_lung, label = erlotinib.labels$AUC_lung)
show_pca(input_data = erlotinib$AUC_lung, label = erlotinib.labels$AUC_lung.source)

# ComBat
erlotinib$AUC_lung.ComBat <- ComBat_combine(batch = erlotinib.labels$AUC_lung.source,
                                            label = erlotinib.labels$AUC_lung,
                                            input_data = erlotinib$AUC_lung)
mean(erlotinib$AUC_lung.ComBat) 
# -0.006164507

#show_pca(input_data = erlotinib$IC50_lung.ComBat, label = erlotinib.labels$IC50_lung)
show_pca(input_data = erlotinib$IC50_lung.ComBat, label = erlotinib.labels$IC50_lung.source)

# get l1000 features ------------------------------------------------------
stopifnot(colnames(erlotinib$IC50_combined.ComBat) == colnames(erlotinib$AUC_combined.ComBat))
stopifnot(colnames(erlotinib$slope_combined.ComBat) == colnames(erlotinib$AUC_combined.ComBat))

Landmark_Genes_n978 <- read.csv("Common/Landmark_Genes_n978.csv")

feature.l1000 <- list()
feature.l1000$cp <- which(colnames(erlotinib$AUC_combined.ComBat) %in% Landmark_Genes_n978$Ensembl)
feature.l1000$pp <- which(colnames(erlotinib$patient) %in% Landmark_Genes_n978$Ensembl)

stopifnot(length(feature.l1000$cp) > 0)
stopifnot(length(feature.l1000$pp) > 0)
length(feature.l1000$cp)
length(feature.l1000$pp)

# create partitions for 13 to 24 patients using all cell lines ------------------------
partition <- list()

input.labels_cell_lines = list()
input.labels_cell_lines$slope = erlotinib.labels$slope_combined[which(erlotinib.labels$slope_combined.source != "patient")]
input.labels_cell_lines$IC50 = erlotinib.labels$IC50_combined[which(erlotinib.labels$IC50_combined.source != "patient")]
input.labels_cell_lines$AUC = erlotinib.labels$AUC_combined[which(erlotinib.labels$AUC_combined.source != "patient")]

input.cell_line_order = list()
input.cell_line_order$slope = 1:length(input.labels_cell_lines$slope)
input.cell_line_order$IC50 = 1:length(input.labels_cell_lines$IC50)
input.cell_line_order$AUC = 1:length(input.labels_cell_lines$AUC)

stopifnot(input.cell_line_order$slope == input.cell_line_order$IC50)  
stopifnot(input.cell_line_order$AUC == input.cell_line_order$IC50)  
stopifnot(length(input.cell_line_order$AUC) > 0)

cell_lines_all <- foreach(input.training_amount.p = seq(from = 13, to = 24, by = 1)
                          , .errorhandling = "stop") %dopar% {
                            
                            generate_random_partition(labels_cell_lines = input.labels_cell_lines, 
                                                      labels_patient = erlotinib.labels$patient,
                                                      training_amount.p = input.training_amount.p,
                                                      cell_line_order = input.cell_line_order,
                                                      training_amount.c = length(input.cell_line_order$AUC),
                                                      acc_training = TRUE
                            )
                          }

partition$cell_lines_all <- cell_lines_all

# using cell lines based on lusc.comsic similarities -----------------------------
load("CGP/cosmic.tcga.RData")

input.labels_cell_lines = list()
input.labels_cell_lines$slope = erlotinib.labels$slope_combined[which(erlotinib.labels$slope_combined.source != "patient")]
input.labels_cell_lines$IC50 = erlotinib.labels$IC50_combined[which(erlotinib.labels$IC50_combined.source != "patient")]
input.labels_cell_lines$AUC = erlotinib.labels$AUC_combined[which(erlotinib.labels$AUC_combined.source != "patient")]

input.cell_line_order = list()
input.cell_line_order$slope <- order(match(names(input.labels_cell_lines$slope),lusc_ordered$cell.line))
input.cell_line_order$IC50 <- order(match(names(input.labels_cell_lines$IC50), lusc_ordered$cell.line))
input.cell_line_order$AUC <- order(match(names(input.labels_cell_lines$AUC), lusc_ordered$cell.line))

stopifnot(length(input.cell_line_order$slope) == length(input.cell_line_order$IC50))
stopifnot(length(input.cell_line_order$AUC) == length(input.cell_line_order$IC50))
stopifnot(length(input.cell_line_order$AUC) > 0)

stopifnot( length(table(input.labels_cell_lines$slope[input.cell_line_order$slope][1:30])) == 2 )
stopifnot( length(table(input.labels_cell_lines$AUC[input.cell_line_order$AUC][1:30])) == 2 )
stopifnot( length(table(input.labels_cell_lines$IC50[input.cell_line_order$IC50][1:30])) == 2 )

stopifnot(length(input.labels_cell_lines$slope) == 525)

patient_24 <- foreach(input.training_amount.c = seq(from = 50, to = 525, by = 30), .errorhandling = "stop") %dopar% {
  generate_random_partition(labels_cell_lines = input.labels_cell_lines, 
                            labels_patient = erlotinib.labels$patient,
                            training_amount.p = 24,
                            cell_line_order = input.cell_line_order,
                            training_amount.c = input.training_amount.c,
                            acc_training = TRUE)
}
# expect all leave-one-out warnigns
partition$patient_24 <- patient_24

# create partitions for 13 to 24 patients using lung cell lines ------------------------
input.labels_cell_lines = list()
input.labels_cell_lines$slope = erlotinib.labels$slope_lung[which(erlotinib.labels$slope_lung.source != "patient")]
input.labels_cell_lines$IC50 = erlotinib.labels$IC50_lung[which(erlotinib.labels$IC50_lung.source != "patient")]
input.labels_cell_lines$AUC = erlotinib.labels$AUC_lung[which(erlotinib.labels$AUC_lung.source != "patient")]

input.cell_line_order = list()
input.cell_line_order$slope = 1:length(input.labels_cell_lines$slope)
input.cell_line_order$IC50 = 1:length(input.labels_cell_lines$IC50)
input.cell_line_order$AUC = 1:length(input.labels_cell_lines$AUC)

stopifnot(input.cell_line_order$slope == input.cell_line_order$IC50)  
stopifnot(input.cell_line_order$AUC == input.cell_line_order$IC50)  
stopifnot(length(input.cell_line_order$AUC) > 0)

cell_lines_lung <- foreach(input.training_amount.p = seq(from = 13, to = 24, by = 1)
                           , .errorhandling = "stop") %dopar% {
                             
                             generate_random_partition(labels_cell_lines = input.labels_cell_lines, 
                                                       labels_patient = erlotinib.labels$patient,
                                                       training_amount.p = input.training_amount.p,
                                                       cell_line_order = input.cell_line_order,
                                                       training_amount.c = length(input.cell_line_order$AUC),
                                                       acc_training = TRUE,
                                                       input_partition = cell_lines_all[[input.training_amount.p - 12]]
                             )
                           }

partition$cell_lines_lung <- cell_lines_lung

names(partition)
#[1] "cell_lines_all"  "patient_24"      "cell_lines_lung"
names(erlotinib)
names(erlotinib.labels)
# save worksapce ----------------------------------------------------------
save(sampleinfo.ccle, erlotinib, erlotinib.labels, feature.l1000, partition,
     file = "Erlotinib/WS/erlotinib_homogenized_data_ccle.RData")
