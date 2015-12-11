library("doParallel")
library("sva")

source('Epirubicin/Script/generate_random_partition.R')
source('Common/preparing_data_helper.R')
source('Common/drug_cut/callingWaterfall.R')
source('Common/drug_cut/distancePointLine.R')
source('Common/drug_cut/distancePointSegment.R')
source('Common/comGENE.R')

###

load("~/Dropbox/SNF_DRUG_PROJECT/CELL_LINE/su2c_frma.RData")

temp.drug_id = rownames(druginfo_su2c)[match("Epirubicin", druginfo_su2c$name)]
stopifnot(temp.drug_id == "drug_su2c.20")

temp.median <- median(pheno_su2c$drug_su2c.20, na.rm = TRUE)
temp.include_ind <- which(!is.na(pheno_su2c$drug_su2c.20))

epirubicin <- list()

epirubicin$heiser <- data_su2c
stopifnot(length(colnames(epirubicin$heiser)) == length(annot_su2c$hgnc_symbol))
colnames(epirubicin$heiser) <- annot_su2c$hgnc_symbol

epirubicin$heiser <- epirubicin$heiser[temp.include_ind, ]

dim(epirubicin$heiser)

#### Get the new Sensitivity Data #### 
epirubicin.labels <- list()
gray_sensitivity_2 <- read.csv("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/New_sensitivity_measurments/gray_sensitivity_2.csv")
temp.epirubicin_ind <- which(gray_sensitivity_2$drug == "Epirubicin")
gray_sensitivity_2 <- gray_sensitivity_2[temp.epirubicin_ind, ]
dim(gray_sensitivity_2)

temp.na_ind <- which(!is.na(gray_sensitivity_2$slope0.sensitivity.call))
length(temp.na_ind)
temp.response_slope <- gray_sensitivity_2$slope0.sensitivity.call[temp.na_ind]
length(temp.response_slope) #48
stopifnot(sum(is.na(temp.response_slope)) == 0)

### get the data
temp.ind <- which(rownames(epirubicin$heiser) %in% gray_sensitivity_2$cellline)
epirubicin$heiser <- epirubicin$heiser[temp.ind, ]
temp.ind <- match(rownames(epirubicin$heiser), gray_sensitivity_2$cellline)
length(temp.ind)

# slope
epirubicin.labels$slope <- gray_sensitivity_2$slope0.sensitivity.call[temp.ind]
names(epirubicin.labels$slope) <- gray_sensitivity_2$cellline[temp.ind]
stopifnot(sum(is.na(epirubicin.labels$slope)) == 0)
epirubicin.labels$slope = epirubicin.labels$slope == 1
table(epirubicin.labels$slope)

### get AUC labels
temp.response <- callingWaterfall(gray_sensitivity_2$AUC[temp.ind], type="AUC")
table(temp.response)
temp.response <- temp.response != "resistant"
names(temp.response) <- gray_sensitivity_2$cellline[temp.ind]
table(epirubicin.labels$AUC)
epirubicin.labels$AUC <- temp.response
length(epirubicin.labels$AUC)
stopifnot(sum(is.na(epirubicin.labels$AUC)) == 0)
sum(epirubicin.labels$AUC == epirubicin.labels$slope)
# 35 / 38 for AUC and slope

### get the IC50 labels
temp.response <- callingWaterfall(gray_sensitivity_2$IC50_Published[temp.ind], type="IC50")
table(temp.response)
temp.response <- temp.response != "resistant"
names(temp.response) <- gray_sensitivity_2$cellline[temp.ind]

epirubicin.labels$IC50 <- temp.response
length(epirubicin.labels$IC50)
epirubicin$IC50 <- epirubicin$heiser[-which(is.na(epirubicin.labels$IC50)), ]
epirubicin.labels$IC50 <- epirubicin.labels$IC50[-which(is.na(epirubicin.labels$IC50))]
stopifnot(sum(is.na(epirubicin.labels$IC50)) == 0)

temp.ind <- match(names(epirubicin.labels$IC50), names(epirubicin.labels$slope))
sum(epirubicin.labels$IC50 == epirubicin.labels$slope[temp.ind])
# 26 / 36 for IC50 and slope
temp.ind <- match(names(epirubicin.labels$IC50), names(epirubicin.labels$AUC))
sum(epirubicin.labels$IC50 == epirubicin.labels$AUC[temp.ind])
# 25 / 36 for IC50 and slope


table(epirubicin.labels$IC50)

###
load("Epirubicin/WS/desmedt2009.RData")
dim(demo)

### get the labels for epirubicin
table(demo$pCR)
temp.include_ind = which(!is.na(demo$pCR))

stopifnot(rownames(data) == demo$samplename)
epirubicin$patient <- data[temp.include_ind, ]
epirubicin.labels$patient <- demo$pCR[temp.include_ind]
epirubicin.labels$patient[which(epirubicin.labels$patient == "YES")] = TRUE
epirubicin.labels$patient[which(epirubicin.labels$patient == "NO")] = FALSE

epirubicin.labels$patient <- as.logical(epirubicin.labels$patient)

names(epirubicin.labels$patient) <- demo$samplename[temp.include_ind]
table(epirubicin.labels$patient)

stopifnot(colnames(epirubicin$patient) == annot$probe)
colnames(epirubicin$patient) <- annot$Gene.symbol

### combine slope and patient
temp.data <- comGENE(scale(epirubicin$heiser), scale(epirubicin$patient))
dim(temp.data[[1]])
dim(temp.data[[2]])
mean(temp.data[[1]])
mean(temp.data[[2]])

epirubicin.labels$slope_combined.source <- c(rep("heiser", dim(temp.data[[1]])[1]), rep("patient", dim(temp.data[[2]])[1]))
epirubicin.labels$slope_combined <- c(epirubicin.labels$slope, epirubicin.labels$patient)
table(epirubicin.labels$slope_combined)
cheng <- rbind(scale(temp.data[[1]]), scale(temp.data[[2]]))
show_pca(input_data = epirubicin$combined_slope.sva, label = epirubicin.labels$slope_combined)
show_pca(input_data = cheng, label = epirubicin.labels$slope_combined.source)

epirubicin$combined_slope.sva <- sva_combine(batch = epirubicin.labels$slope_combined.source,
                                            label = epirubicin.labels$slope_combined, 
                                              input_data = rbind(scale(temp.data[[1]]), scale(temp.data[[2]])), n.sv = 3) 
show_pca(input_data = epirubicin$combined_slope.sva, label = epirubicin.labels$slope_combined)
show_pca(input_data = epirubicin$combined_slope.sva, label = epirubicin.labels$slope_combined.source)

epirubicin$combined_slope.combat <- ComBat_combine(batch = epirubicin.labels$slope_combined.source,
                         label = epirubicin.labels$slope_combined, 
                         input_data = rbind(scale(temp.data[[1]]), scale(temp.data[[2]])),  par.prior = FALSE)
show_pca(input_data = epirubicin$combined_slope.combat, label = epirubicin.labels$slope_combined)
show_pca(input_data = epirubicin$combined_slope.combat, label = epirubicin.labels$slope_combined.source)
rm(temp.data)
# used ComBat

### combine AUC and patient
temp.data <- comGENE(scale(epirubicin$heiser), scale(epirubicin$patient))
dim(temp.data[[1]])
dim(temp.data[[2]])
mean(temp.data[[1]])
mean(temp.data[[2]])

epirubicin.labels$AUC_combined.source <- c(rep("heiser", dim(temp.data[[1]])[1]), rep("patient", dim(temp.data[[2]])[1]))
epirubicin.labels$AUC_combined <- c(epirubicin.labels$AUC, epirubicin.labels$patient)
table(epirubicin.labels$AUC_combined)
epirubicin$combined_AUC.sva <- sva_combine(batch = epirubicin.labels$AUC_combined.source,
                                             label = epirubicin.labels$AUC_combined, 
                                             input_data = rbind(scale(temp.data[[1]]), scale(temp.data[[2]]))
                                             , n.sv = 3)
show_pca(input_data = epirubicin$combined_AUC.sva, label = epirubicin.labels$AUC_combined)
show_pca(input_data = epirubicin$combined_AUC.sva, label = epirubicin.labels$AUC_combined.source)

epirubicin$combined_AUC.combat <- ComBat_combine(batch = epirubicin.labels$AUC_combined.source,
                                                label = epirubicin.labels$AUC_combined, 
                                                input_data = rbind(scale(temp.data[[1]]), scale(temp.data[[2]])))
show_pca(input_data = epirubicin$combined_AUC.combat, label = epirubicin.labels$AUC_combined)
show_pca(input_data = epirubicin$combined_AUC.combat, label = epirubicin.labels$AUC_combined.source)
rm(temp.data)
# combat

### combine IC50 and patient
source('~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Script/preparing_data_helper.R')
source('~/Dropbox/SNF_DRUG_PROJECT/Script/comGENE.R')
temp.data <- comGENE(scale(epirubicin$IC50), scale(epirubicin$patient))
dim(temp.data[[1]])
dim(temp.data[[2]])
mean(temp.data[[1]])
mean(temp.data[[2]])

epirubicin.labels$IC50_combined.source <- c(rep("heiser", dim(temp.data[[1]])[1]), rep("patient", dim(temp.data[[2]])[1]))
epirubicin.labels$IC50_combined <- c(epirubicin.labels$IC50, epirubicin.labels$patient)
table(epirubicin.labels$IC50_combined)
epirubicin$combined_IC50.sva <- sva_combine(batch = epirubicin.labels$IC50_combined.source,
                                              label = epirubicin.labels$IC50_combined, 
                                              input_data = rbind(scale(temp.data[[1]]), scale(temp.data[[2]]))
                                              , n.sv = 3)
show_pca(input_data = epirubicin$combined_IC50.sva, label = epirubicin.labels$IC50_combined)
show_pca(input_data = epirubicin$combined_IC50.sva, label = epirubicin.labels$IC50_combined.source)
epirubicin$combined_IC50.combat <- ComBat_combine(batch = epirubicin.labels$IC50_combined.source,
                                                 label = epirubicin.labels$IC50_combined, 
                                                 input_data = rbind(scale(temp.data[[1]]), scale(temp.data[[2]])))
show_pca(input_data = epirubicin$combined_IC50.combat, label = epirubicin.labels$IC50_combined)
show_pca(input_data = epirubicin$combined_IC50.combat, label = epirubicin.labels$IC50_combined.source)
rm(temp.data)
# combat

### get the partitioning
temp.test_amount <- list(cpp = round(0.2 * length(epirubicin.labels$patient)))

partition <- list()
partition$slope <- generate_random_partition(input_labels_cell_lines = epirubicin.labels$slope, 
                                             input_labels_patient = epirubicin.labels$patient, 
                                             test_amount = temp.test_amount)

partition$IC50 <- generate_random_partition(input_labels_cell_lines = epirubicin.labels$IC50, 
                                            input_labels_patient = epirubicin.labels$patient, 
                                            test_amount = temp.test_amount)

partition$AUC <- generate_random_partition(input_labels_cell_lines = epirubicin.labels$AUC, 
                                           input_labels_patient = epirubicin.labels$patient, 
                                           test_amount = temp.test_amount)


## get the L1000 genes
Landmark_Genes_n978 <- read.csv("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Landmark_Genes_n978.csv")
feature.l1000 <- list()
feature.l1000$cc <- which(colnames(epirubicin$heiser) %in% Landmark_Genes_n978$Gene.Symbol)
length(feature.l1000$cc)
feature.l1000$cp <- which(colnames(epirubicin$combined_slope.combat) %in% Landmark_Genes_n978$Gene.Symbol)
length(feature.l1000$cp)

source('~/Dropbox/SNF_DRUG_PROJECT/Script/comGENE.R')
temp.data <- comGENE(epirubicin$patient, epirubicin$patient)
epirubicin$patient <- temp.data[[1]]

temp.variances <- apply(epirubicin$patient, 2, var)
quantile(temp.variances, 0.1)
temp.low_var_genes <- which(temp.variances < quantile(temp.variances, 0.1))
length(temp.low_var_genes)
epirubicin$patient <- epirubicin$patient[, -temp.low_var_genes]
dim(epirubicin$patient)
feature.l1000$pp <- which(colnames(epirubicin$patient) %in% Landmark_Genes_n978$Gene.Symbol)
length(feature.l1000$pp)

### create partitions for varying number of patients
partition_var <- list()

temp.cp <- foreach (training_amount.cp = seq(from = 20, to = 38, by = 5)) %do% {    
  
  temp.slope <- generate_random_partition.cp_var(input_labels_cell_lines = epirubicin.labels$slope, 
                                                 input_labels_patient = epirubicin.labels$patient, 
                                                 training_amount = training_amount.cp)  
  
  temp.AUC <- generate_random_partition.cp_var(input_labels_cell_lines = epirubicin.labels$AUC, 
                                               input_labels_patient = epirubicin.labels$patient, 
                                               training_amount = training_amount.cp)
  list(slope = temp.slope, AUC = temp.AUC)
}
length(temp.cp) #6
# double check
for (temp.ind in 1:length(temp.cp)) {
  stopifnot(length(temp.cp[[temp.ind]]$slope[[100]]$test_index) == length(epirubicin.labels$patient))
  stopifnot(length(temp.cp[[temp.ind]]$slope[[100]]$training_index.single) == temp.ind * 5 + 15)
}

partition_var$cp <- temp.cp

temp.cpp <- foreach (training_amount.cpp = seq(from = 10, to = 100, by = 10)) %do% {    
  
  temp.slope <- generate_random_partition.cpp_var(input_labels_cell_lines = epirubicin.labels$slope, 
                                                  input_labels_patient = epirubicin.labels$patient, 
                                                  training_amount = training_amount.cpp)  
  
  temp.AUC <- generate_random_partition.cpp_var(input_labels_cell_lines = epirubicin.labels$AUC, 
                                                input_labels_patient = epirubicin.labels$patient, 
                                                training_amount = training_amount.cpp)
  list(slope = temp.slope, AUC = temp.AUC)
}

partition_var$cpp <- temp.cpp

# double check
for (temp.ind in 1:10) {
  stopifnot(length(temp.cpp[[temp.ind]]$slope[[100]]$test_index) == length(epirubicin.labels$patient) - temp.ind * 10 )
  stopifnot(length(temp.cpp[[temp.ind]]$slope[[100]]$training_index.single) == temp.ind * 10 + length(epirubicin.labels$slope))
}

## patients
temp.pp <- foreach (training_amount.pp = seq(from = 20, to = 100, by = 10)) %do% {      
  generate_random_partition.pp_var(input_labels_patient = epirubicin.labels$patient, 
                                   training_amount = training_amount.pp)
}

# double check
for (temp.ind in 1:9) {
  stopifnot(length(temp.pp[[temp.ind]][[100]]$test_index) == length(epirubicin.labels$patient) - temp.ind * 10 - 10)
  stopifnot(length(temp.pp[[temp.ind]][[100]]$training_index.single) == temp.ind * 10 + 10)
}

partition_var$pp <- temp.pp

###
temp.cp100p <- foreach (cell_line_training_amount = seq(from = 5, to = 35, by = 5)) %do% {    
  
  temp.slope <- generate_random_partition.cpp_var2(input_labels_cell_lines = epirubicin.labels$slope, 
                                                  input_labels_patient = epirubicin.labels$patient, 
                                                  patient_training_amount = 100,
                                                  cell_line_training_amount = cell_line_training_amount
                                                  )  
  
  temp.AUC <- generate_random_partition.cpp_var2(input_labels_cell_lines = epirubicin.labels$AUC, 
                                                input_labels_patient = epirubicin.labels$patient,
                                                cell_line_training_amount = cell_line_training_amount,
                                                patient_training_amount = 100)
  list(slope = temp.slope, AUC = temp.AUC)
}

partition_var$cVar_p100_p <- temp.cp100p

###
temp.cp80p <- foreach (cell_line_training_amount = seq(from = 5, to = 35, by = 5)) %do% {    
  
  temp.slope <- generate_random_partition.cpp_var2(input_labels_cell_lines = epirubicin.labels$slope, 
                                                   input_labels_patient = epirubicin.labels$patient, 
                                                   patient_training_amount = 80,
                                                   cell_line_training_amount = cell_line_training_amount
  )  
  
  temp.AUC <- generate_random_partition.cpp_var2(input_labels_cell_lines = epirubicin.labels$AUC, 
                                                 input_labels_patient = epirubicin.labels$patient,
                                                 cell_line_training_amount = cell_line_training_amount,
                                                 patient_training_amount = 80)
  list(slope = temp.slope, AUC = temp.AUC)
}

partition_var$cVar_p80_p <- temp.cp80p


###

temp.cp50p <- foreach (cell_line_training_amount = seq(from = 5, to = 35, by = 5)) %do% {    
  
  temp.slope <- generate_random_partition.cpp_var2(input_labels_cell_lines = epirubicin.labels$slope, 
                                                   input_labels_patient = epirubicin.labels$patient, 
                                                   patient_training_amount = 50,
                                                   cell_line_training_amount = cell_line_training_amount
  )  
  
  temp.AUC <- generate_random_partition.cpp_var2(input_labels_cell_lines = epirubicin.labels$AUC, 
                                                 input_labels_patient = epirubicin.labels$patient,
                                                 cell_line_training_amount = cell_line_training_amount,
                                                 patient_training_amount = 50)
  list(slope = temp.slope, AUC = temp.AUC)
}

partition_var$cVar_p50_p <- temp.cp50p

####
load("CGP/cosmic.tcga.RData")

cell_line_order <- list()
cell_line_order$AUC <- order(match(names(epirubicin.labels$AUC), brca_ordered$cell.line))
cell_line_order$slope <- order(match(names(epirubicin.labels$slope), brca_ordered$cell.line))

stopifnot(!is.unsorted(match(brca_ordered$cell.line, names(epirubicin.labels$slope)[cell_line_order$AUC] ), na.rm = TRUE))

stopifnot( length(table(epirubicin.labels$AUC[cell_line_order$AUC][1:5])) == 2 )
stopifnot( length(table(epirubicin.labels$slope[cell_line_order$slope][1:5])) == 2 )

# called this with 50 and 80
temp.cp100p <- foreach (cell_line_training_amount = seq(from = 5, to = 35, by = 5)) %do% {    
  
  temp.slope <- generate_random_partition.cpp_var3(input_labels_cell_lines = epirubicin.labels$slope, 
                                                   input_labels_patient = epirubicin.labels$patient, 
                                                   patient_training_amount = 100,
                                                   cell_line_training_amount = cell_line_training_amount,
                                                   cell_line_order = cell_line_order$slope
  )  
  
  temp.AUC <- generate_random_partition.cpp_var3(input_labels_cell_lines = epirubicin.labels$AUC, 
                                                 input_labels_patient = epirubicin.labels$patient,
                                                 cell_line_training_amount = cell_line_training_amount,
                                                 patient_training_amount = 100,
                                                 cell_line_order = cell_line_order$AUC)
  list(slope = temp.slope, AUC = temp.AUC)
}

partition_var$cVar_p50_p_brca <- temp.cp50p
partition_var$cVar_p80_p_brca <- temp.cp80p
partition_var$cVar_p100_p_brca <- temp.cp100p



###########
#save(epirubicin, epirubicin.labels, partition, feature.l1000, partition_var, file = "Epirubicin/WS/epirubicin_data.RData")

