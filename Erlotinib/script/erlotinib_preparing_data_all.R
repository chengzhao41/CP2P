# Load Libraries ----------------------------------------------------------
library("doParallel")

source('Common/preparing_data_helper.R')
source('Common/comGENE.R')
source("Common/generate_random_partition.R")

# Load GDSC Cell line data --------------------------------------------
load("Erlotinib/WS/erlotinib_gdsc.RData")
erlotinib.gdsc <- erlotinib
erlotinib.labels.gdsc <- erlotinib.labels

# Load CCLE Cell line data --------------------------------------------
load("Erlotinib/WS/erlotinib_ccle.RData")

erlotinib$gdsc_slope <- erlotinib.gdsc$gdsc_slope
erlotinib$gdsc_AUC <- erlotinib.gdsc$gdsc_AUC
erlotinib$gdsc_IC50 <- erlotinib.gdsc$gdsc_IC50

erlotinib.labels$gdsc_slope <- erlotinib.labels.gdsc$slope
erlotinib.labels$gdsc_AUC <- erlotinib.labels.gdsc$AUC
erlotinib.labels$gdsc_IC50 <- erlotinib.labels.gdsc$IC50

erlotinib.labels$gdsc_slope_ind <- erlotinib.labels.gdsc$slope_ind  
erlotinib.labels$gdsc_AUC_ind <- erlotinib.labels.gdsc$AUC_ind
erlotinib.labels$gdsc_IC50_ind <- erlotinib.labels.gdsc$IC50_ind

rm(erlotinib.labels.gdsc, erlotinib.gdsc)

# sanity check to make sure we have the right version of data
mean(erlotinib$gdsc_slope) # 6.293544
mean(erlotinib$gdsc_AUC) # 6.293544
mean(erlotinib$gdsc_IC50) # 6.293544
dim(erlotinib$gdsc_slope) # 287 11833
dim(erlotinib$gdsc_AUC) # 287 11833
dim(erlotinib$gdsc_IC50) # 287 11833

table(erlotinib.labels$gdsc_slope)
# FALSE  TRUE 
# 218    69 
table(erlotinib.labels$gdsc_AUC)
# FALSE  TRUE 
# 215    72 
table(erlotinib.labels$gdsc_IC50)
# FALSE  TRUE 
# 277    10

mean(erlotinib$ccle_slope) # 6.045301
mean(erlotinib$ccle_AUC) # 6.045301
mean(erlotinib$ccle_IC50) # 6.045301
dim(erlotinib$ccle_slope) # 489 20049
dim(erlotinib$ccle_AUC) # 489 20049
dim(erlotinib$ccle_IC50) # 489 20049

erlotinib.labels$ccle_slope <- erlotinib.labels$slope
erlotinib.labels$ccle_AUC <- erlotinib.labels$AUC
erlotinib.labels$ccle_IC50 <- erlotinib.labels$IC50

erlotinib.labels$ccle_slope_ind <- erlotinib.labels$slope_ind
erlotinib.labels$ccle_AUC_ind <- erlotinib.labels$AUC_ind
erlotinib.labels$ccle_IC50_ind <- erlotinib.labels$IC50_ind

erlotinib.labels$slope <- NULL
erlotinib.labels$AUC <- NULL
erlotinib.labels$IC50 <- NULL
erlotinib.labels$AUC.new <- NULL
erlotinib.labels$slope.new <- NULL
erlotinib.labels$IC50.new <- NULL
erlotinib.labels$IC50.cont <- NULL
erlotinib.labels$AUC.new <- NULL
erlotinib.labels$AUC.cont <- NULL
erlotinib.labels$slope.new <- NULL
erlotinib.labels$slope.cont <- NULL

erlotinib.labels$AUC_ind <- NULL
erlotinib.labels$slope_ind <- NULL
erlotinib.labels$IC50_ind <- NULL

table(erlotinib.labels$ccle_slope)
# FALSE  TRUE 
# 396    93 
table(erlotinib.labels$ccle_AUC)
# FALSE  TRUE 
# 385   104 
table(erlotinib.labels$ccle_IC50)
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

temp.data2 <- comGENE(temp.data1[[1]], scale(erlotinib$gdsc_slope))
mean(temp.data2[[1]]) #-1.158721e-17
mean(temp.data2[[2]]) #1.103038e-17
dim(temp.data2[[1]]) #25 10752
dim(temp.data2[[2]]) #287 10752

temp.data3 <- comGENE(temp.data1[[2]], temp.data2[[2]])
mean(temp.data3[[1]]) #-2.654878e-18
mean(temp.data3[[2]]) #1.103038e-17
dim(temp.data3[[1]]) #44 10752
dim(temp.data3[[2]]) #287 10752

temp.patient <- temp.data2[[1]]
erlotinib$battle <- temp.data3[[1]]
erlotinib$gdsc_slope <- temp.data3[[2]]
stopifnot(colnames(temp.patient) == colnames(erlotinib$battle))
stopifnot(colnames(erlotinib$gdsc_slope) == colnames(erlotinib$battle))

temp.data4 <- comGENE(erlotinib$gdsc_slope, scale(erlotinib$ccle_slope))
mean(temp.data4[[1]]) #-2.654878e-18
mean(temp.data4[[2]]) #1.103038e-17
dim(temp.data4[[1]]) #287 10752
dim(temp.data4[[2]]) #489 10752

erlotinib$ccle_slope <- temp.data4[[2]]

dim(erlotinib$ccle_slope) == dim(erlotinib$gdsc_slope)

# find common set of samples, and take them out from battle ----------------------------------------------
temp.common <- intersect(names(erlotinib.labels$ccle_AUC), names(erlotinib.labels$battle_IC50))
erlotinib$battle <- erlotinib$battle[-which(names(erlotinib.labels$battle_IC50) %in% temp.common), ]
erlotinib.labels$battle_IC50 <- erlotinib.labels$battle_IC50[-which(names(erlotinib.labels$battle_IC50) %in% temp.common)]

rm(temp.common)
temp.common <- intersect(names(erlotinib.labels$ccle_AUC), names(erlotinib.labels$gdsc_slope))
erlotinib$ccle_AUC <- erlotinib$ccle_AUC[-which(names(erlotinib.labels$ccle_AUC) %in% temp.common), ]
erlotinib$ccle_IC50 <- erlotinib$ccle_IC50[-which(names(erlotinib.labels$ccle_IC50) %in% temp.common), ]
erlotinib$ccle_slope <- erlotinib$ccle_slope[-which(names(erlotinib.labels$ccle_slope) %in% temp.common), ]

erlotinib.labels$ccle_AUC <- erlotinib.labels$ccle_AUC[-which(names(erlotinib.labels$ccle_AUC) %in% temp.common)]
erlotinib.labels$ccle_IC50 <- erlotinib.labels$ccle_IC50[-which(names(erlotinib.labels$ccle_IC50) %in% temp.common)]
erlotinib.labels$ccle_slope <- erlotinib.labels$ccle_slope[-which(names(erlotinib.labels$ccle_slope) %in% temp.common)]

erlotinib.labels$ccle_AUC_ind <- erlotinib.labels$ccle_AUC_ind[-which(names(erlotinib.labels$ccle_AUC) %in% temp.common)]
erlotinib.labels$ccle_IC50_ind <- erlotinib.labels$ccle_IC50_ind[-which(names(erlotinib.labels$ccle_IC50) %in% temp.common)]
erlotinib.labels$ccle_slope_ind <- erlotinib.labels$ccle_slope_ind[-which(names(erlotinib.labels$ccle_slope) %in% temp.common)]

# Combining with slope labels ---------------------------------
mean(temp.patient) #-1.158721e-17
mean(erlotinib$battle) #0.0009030523
mean(erlotinib$gdsc_slope) #1.103038e-17
mean(erlotinib$ccle_slope) #-0.001864834

erlotinib$slope_combined <- rbind(temp.patient, erlotinib$battle, 
                                  erlotinib$gdsc_slope, erlotinib$ccle_slope)
erlotinib.labels$slope_combined <- c(erlotinib.labels$patient, erlotinib.labels$battle_IC50, 
                                     erlotinib.labels$gdsc_slope, erlotinib.labels$ccle_slope)
erlotinib.labels$slope_combined.source <- c(rep("patient", dim(temp.patient)[1]), 
                                            rep("battle", dim(erlotinib$battle)[1]), 
                                            rep("gdsc", dim(erlotinib$gdsc_slope)[1]),
                                            rep("ccle", dim(erlotinib$ccle_slope)[1]))

mean(erlotinib$slope_combined) #-0.0009763436
#show_pca(input_data = erlotinib$slope_combined, label = erlotinib.labels$slope_combined)
#show_pca(input_data = erlotinib$slope_combined, label = erlotinib.labels$slope_combined.source)

# Combining with gdsc IC50 labels   ---------------------------------
erlotinib$IC50_combined <- erlotinib$slope_combined
erlotinib.labels$IC50_combined <- c(erlotinib.labels$patient, erlotinib.labels$battle_IC50, 
                                    erlotinib.labels$gdsc_IC50, erlotinib.labels$ccle_IC50)
erlotinib.labels$IC50_combined.source <- c(rep("patient", dim(temp.patient)[1]), 
                                           rep("battle", dim(erlotinib$battle)[1]), 
                                           rep("gdsc", dim(erlotinib$gdsc_IC50)[1]),
                                           rep("ccle", dim(erlotinib$ccle_IC50)[1]))
#show_pca(input_data = erlotinib$IC50_combined, label = erlotinib.labels$IC50_combined)
#show_pca(input_data = erlotinib$IC50_combined, label = erlotinib.labels$IC50_combined.source)

# Combining with gdsc AUC labels   ---------------------------------
erlotinib$AUC_combined <- erlotinib$slope_combined
erlotinib.labels$AUC_combined <- c(erlotinib.labels$patient, erlotinib.labels$battle_IC50, 
                                   erlotinib.labels$gdsc_AUC, erlotinib.labels$ccle_AUC)
erlotinib.labels$AUC_combined.source <- c(rep("patient", dim(temp.patient)[1]), 
                                          rep("battle", dim(erlotinib$battle)[1]), 
                                          rep("gdsc", dim(erlotinib$gdsc_AUC)[1]),
                                          rep("ccle", dim(erlotinib$ccle_AUC)[1]))
#show_pca(input_data = erlotinib$AUC_combined, label = erlotinib.labels$AUC_combined)
#show_pca(input_data = erlotinib$AUC_combined, label = erlotinib.labels$AUC_combined.source)

# Combining with only lung cell lines ------------------------------------------------------
mean(temp.patient) #-1.158721e-17
mean(erlotinib$battle) #0.0009030523
mean(erlotinib$gdsc_slope) #1.103038e-17

gdsc.lung.ind <- which(sampleinfo.gdsc$tissue.type[erlotinib.labels$gdsc_slope_ind] == "lung")
ccle.lung.ind <- which(sampleinfo.ccle$tissueid[erlotinib.labels$ccle_slope_ind] == "lung")

erlotinib$slope_lung <- rbind(temp.patient, 
                              erlotinib$battle, 
                              erlotinib$gdsc_slope[gdsc.lung.ind, ],
                              erlotinib$ccle_slope[ccle.lung.ind, ])
erlotinib.labels$slope_lung <- c(erlotinib.labels$patient, 
                                 erlotinib.labels$battle_IC50, 
                                 erlotinib.labels$gdsc_slope[gdsc.lung.ind],
                                 erlotinib.labels$ccle_slope[ccle.lung.ind]
                                 )
erlotinib.labels$slope_lung.source <- c(rep("patient", dim(temp.patient)[1]), 
                                        rep("battle", dim(erlotinib$battle)[1]), 
                                        rep("gdsc", dim(erlotinib$gdsc_slope[gdsc.lung.ind, ])[1]),
                                        rep("ccle", dim(erlotinib$ccle_slope[ccle.lung.ind, ])[1])
                                        )
mean(erlotinib$slope_lung) #-0.001877834
#show_pca(input_data = erlotinib$slope_lung, label = erlotinib.labels$slope_lung)
#show_pca(input_data = erlotinib$slope_lung, label = erlotinib.labels$slope_lung.source)

# Combining with gdsc IC50 labels   ---------------------------------
gdsc.lung.ind <- which(sampleinfo.gdsc$tissue.type[erlotinib.labels$gdsc_IC50_ind] == "lung")
ccle.lung.ind <- which(sampleinfo.ccle$tissueid[erlotinib.labels$ccle_IC50_ind] == "lung")

erlotinib$IC50_lung <- erlotinib$slope_lung
erlotinib.labels$IC50_lung <- c(erlotinib.labels$patient, 
                                erlotinib.labels$battle_IC50, 
                                erlotinib.labels$gdsc_IC50[gdsc.lung.ind],
                                erlotinib.labels$ccle_IC50[ccle.lung.ind])
erlotinib.labels$IC50_lung.source <- c(rep("patient", dim(temp.patient)[1]), 
                                       rep("battle", dim(erlotinib$battle)[1]), 
                                       rep("gdsc", dim(erlotinib$gdsc_IC50[gdsc.lung.ind, ])[1]),
                                       rep("ccle", dim(erlotinib$ccle_IC50[ccle.lung.ind, ])[1]))
mean(erlotinib$IC50_lung) #-0.001877834
#show_pca(input_data = erlotinib$IC50_lung, label = erlotinib.labels$IC50_lung)
#show_pca(input_data = erlotinib$IC50_lung, label = erlotinib.labels$IC50_lung.source)

# Combining with gdsc AUC labels   ---------------------------------
gdsc.lung.ind <- which(sampleinfo.gdsc$tissue.type[erlotinib.labels$gdsc_AUC_ind] == "lung")
ccle.lung.ind <- which(sampleinfo.ccle$tissueid[erlotinib.labels$ccle_AUC_ind] == "lung")

erlotinib$AUC_lung <- erlotinib$slope_lung
erlotinib.labels$AUC_lung <- c(erlotinib.labels$patient, 
                               erlotinib.labels$battle_IC50, 
                               erlotinib.labels$gdsc_AUC[gdsc.lung.ind],
                               erlotinib.labels$ccle_AUC[ccle.lung.ind])
erlotinib.labels$AUC_lung.source <- c(rep("patient", dim(temp.patient)[1]), 
                                      rep("battle", dim(erlotinib$battle)[1]), 
                                      rep("gdsc", dim(erlotinib$gdsc_AUC[gdsc.lung.ind, ])[1]),
                                      rep("ccle", dim(erlotinib$ccle_AUC[ccle.lung.ind, ])[1]))
mean(erlotinib$AUC_lung) #-0.001877834
#show_pca(input_data = erlotinib$AUC_lung, label = erlotinib.labels$AUC_lung)
#show_pca(input_data = erlotinib$AUC_lung, label = erlotinib.labels$AUC_lung.source)

# get l1000 features ------------------------------------------------------
stopifnot(colnames(erlotinib$IC50_combined) == colnames(erlotinib$AUC_combined))
stopifnot(colnames(erlotinib$slope_combined) == colnames(erlotinib$AUC_combined))

Landmark_Genes_n978 <- read.csv("Common/Landmark_Genes_n978.csv")

feature.l1000 <- list()
feature.l1000$cp <- which(colnames(erlotinib$AUC_combined) %in% Landmark_Genes_n978$Ensembl)
feature.l1000$pp <- which(colnames(erlotinib$patient) %in% Landmark_Genes_n978$Ensembl)

stopifnot(length(feature.l1000$cp) > 0)
stopifnot(length(feature.l1000$pp) > 0)
length(feature.l1000$cp) #879
length(feature.l1000$pp) #939

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

cell_lines_all <- foreach(parInd = 1:12, .errorhandling = "stop") %dopar% {
  
  generate_random_partition(labels.cell_lines = input.labels_cell_lines, 
                            labels.patient = erlotinib.labels$patient,
                            num.training.p = seq(13, 24, 1)[parInd],
                            cell_line_order = input.cell_line_order,
                            num.training.c = length(input.cell_line_order$AUC),
                            metric = 'acc',
                            num.min_labels.training.p = 6,
                            num.min_labels.training.c = 6,
                            num.min_labels.test = 1
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

stopifnot(length(input.labels_cell_lines$slope) == 742)

patient_24 <- foreach(parInd = 1:15, .errorhandling = "stop") %dopar% {
                            generate_random_partition(labels.cell_lines = input.labels_cell_lines, 
                                                      labels.patient = erlotinib.labels$patient,
                                                      num.training.p = 24,
                                                      cell_line_order = input.cell_line_order,
                                                      num.training.c = seq(from = 30, to = 742, by = 50)[parInd],
                                                      metric = 'acc',
                                                      num.min_labels.training.p = 6,
                                                      num.min_labels.training.c = 6,
                                                      num.min_labels.test = 1)
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

cell_lines_lung <- foreach(parInd = 1:12, .errorhandling = "stop") %dopar% {
  
  generate_random_partition(labels.cell_lines = input.labels_cell_lines, 
                            labels.patient = erlotinib.labels$patient,
                            num.training.p = seq(13, 24, 1)[parInd],
                            cell_line_order = input.cell_line_order,
                            num.training.c = length(input.cell_line_order$AUC),
                            metric = 'acc',
                            num.min_labels.training.p = 6,
                            num.min_labels.training.c = 6,
                            num.min_labels.test = 1,
                            input_partition = cell_lines_all[[parInd]])
                             
}

partition$cell_lines_lung <- cell_lines_lung

names(partition)
#[1] "cell_lines_all"  "patient_24"      "cell_lines_lung"
names(erlotinib)
# save worksapce ----------------------------------------------------------
save(sampleinfo.gdsc, erlotinib, erlotinib.labels, feature.l1000, partition,
     file = "Erlotinib/WS/erlotinib_homogenized_data_all.RData")
