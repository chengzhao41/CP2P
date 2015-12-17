rootdir <- "~/Dropbox/CP2P"
setwd(rootdir)

source("Common/preparing_data_helper.R")
source("Common/comGENE.R")

load("Bortezomib/WS/bortezomib_data.RData"); bortezomib$patient <- bortezomib$patient.combat
# load("Docetaxel/WS/docetaxel_data.RData")
# load("Erlotinib/WS/erlotinib_data.RData")
# load("Epirubicin/WS/epirubicin_data.RData")


generate_single_random_partition.cpp_var <- function(input_labels_cell_lines, input_labels_patient, training_amount) {
  ## Copied and adjusted from Bortezomib/Script/generate_random_partition.R function generate_random_partition.cpp_var
  stopifnot(input_labels_patient < training_amount)  
  
  temp.cell_lines <- 1:length(input_labels_cell_lines)
  temp.patients <- length(input_labels_cell_lines) + 1:length(input_labels_patient)    
  input_combined_labels <- c(input_labels_cell_lines, input_labels_patient)
  
  temp.patients_in_training <- sample(temp.patients, training_amount)
  temp.training_index.single <- c(temp.cell_lines, temp.patients_in_training)
  temp.test_index <- setdiff(temp.patients, temp.patients_in_training)
  
  while ((length(table(input_combined_labels[temp.training_index.single])) != 2) || (length(table(input_combined_labels[temp.test_index])) != 2)) {
    temp.patients_in_training <- sample(temp.patients, training_amount)
    temp.training_index.single <- c(temp.cell_lines, temp.patients_in_training)
    temp.test_index <- setdiff(temp.patients, temp.patients_in_training)      
  }

  # flags to denote: 0=> cell lines in training; 1=> patients in training; 2=> patients in test
  cpp_flags <- rep("gdsc", length(input_combined_labels))
  cpp_flags[temp.patients_in_training] <- "patients (70%)"
  cpp_flags[temp.test_index] <- "patients 30%"
  
  cpp_source <- rep("gdsc", length(input_combined_labels))
  cpp_source[temp.patients] <- "patient"
    
  stopifnot(length(temp.test_index) > 0)
  stopifnot(length(temp.training_index.single) > 0)    
  stopifnot(length(intersect(temp.test_index, temp.training_index.single)) == 0)
  stopifnot(temp.test_index %in% temp.patients)
  
  cpp.partition <- list(
    test_index = temp.test_index
    , training_index = temp.training_index.single
    , cpp_flags = cpp_flags
    , cpp_source = cpp_source
  )
  
  return(cpp.partition)
}


compute_preprocessing <- function(data, label, sampleinfo, n.sv = 0, preprocessCLsva = FALSE){
  ## Split to training set and patient-only test set
  labels_cell_lines <- label$AUC
  cl.len <- length(label$AUC)
  labels_patient <- label$patient
  training_amount <- round(length(labels_patient) * 0.7) #30% test set
  cpp.partition <- generate_single_random_partition.cpp_var(labels_cell_lines, labels_patient, training_amount)
  #print(length(cpp.partition$test_index))
  #print(length(cpp.partition$training_index))
  
  data.raw <- list()
  labels.combined <- list()
  ## Pool AUC, IC50, Slope data and labels together, preprocess cell lines with SVA first
  if (preprocessCLsva) {
    cgp_AUC.sva <- sva_combine(batch = sampleinfo$tissue.type[label$AUC_ind],  
                                 label = label$AUC, input_data = scale(data$cgp_AUC), n.sv=3)
    temp.data <- comGENE(scale(cgp_AUC.sva), scale(data$patient))
  } else {
    temp.data <- comGENE(scale(data$cgp_AUC), scale(data$patient))
  }
  data.raw[[1]] <- rbind(temp.data[[1]], temp.data[[2]])
  labels.combined[[1]] <- c(label$AUC, label$patient)
  #IC50
  if (preprocessCLsva) {
    cgp_IC50.sva <- sva_combine(batch = sampleinfo$tissue.type[label$IC50_ind],  
                                label = label$IC50, input_data = scale(data$cgp_IC50), n.sv=3)
    temp.data <- comGENE(scale(cgp_IC50.sva), scale(data$patient))
  } else {
    temp.data <- comGENE(scale(data$cgp_IC50), scale(data$patient))
  }
  data.raw[[2]] <- rbind(temp.data[[1]], temp.data[[2]])
  labels.combined[[2]] <- c(label$IC50, label$patient)
  #Slope
  if (preprocessCLsva) {
    cgp_slope.sva <- sva_combine(batch = sampleinfo$tissue.type[label$slope_ind],  
                                 label = label$slope, input_data = scale(data$cgp_slope), n.sv=3)
    temp.data <- comGENE(scale(cgp_slope.sva), scale(data$patient))
  } else {
    temp.data <- comGENE(scale(data$cgp_slope), scale(data$patient))
  }
  data.raw[[3]] <- rbind(temp.data[[1]], temp.data[[2]])
  labels.combined[[3]] <- c(label$slope, label$patient)
  
  pID <- 1
  data.all <- list()
  plabel.all <- list()
  ## Store raw data as is
  # temp.data <- comGENE(scale(data$cgp_AUC), scale(data$patient))
  temp.data <- comGENE(scale(data$cgp_slope), scale(data$patient))
  temp.data.raw <- rbind(temp.data[[1]], temp.data[[2]])
  temp.label.combined <- c(label$slope, label$patient)
  
  ## Compensate for different number of cell lines with AUC/IC50/Slope labels
  temp.ddiff <- length(labels.combined[[1]]) - length(temp.label.combined)
  print(paste("adjusting", temp.ddiff))
  temp.plot_label <- cpp.partition$cpp_flags[(temp.ddiff+1) : length(cpp.partition$cpp_flags)]
  
  data.all[[pID]] <- temp.data.raw
  plabel.all[[pID]] <- temp.plot_label
  pID <- pID + 1
  
  for(outcomeID in 1:3) {
    ## Get data and labels
    temp.data.raw <- data.raw[[outcomeID]]
    temp.labels.combined <- labels.combined[[outcomeID]]
    # print(paste("len", length(temp.data.raw), length(temp.labels.combined)))
    
    ## Compensate for different number of cell lines with AUC/IC50/Slope labels
    temp.ddiff <- length(labels.combined[[1]]) - length(labels.combined[[outcomeID]])
    print(paste("adjusting", temp.ddiff))
    temp.plot_label <- cpp.partition$cpp_flags[(temp.ddiff+1) : length(cpp.partition$cpp_flags)]
    temp.cpp_source <- cpp.partition$cpp_source[(temp.ddiff+1) : length(cpp.partition$cpp_source)]
    # print(paste("len2", length(temp.plot_label), length(temp.cpp_source)))
    # print(paste("len22", length(cpp.partition$cpp_flags), length(cpp.partition$cpp_source)))
    temp.train_index <- cpp.partition$training_index[cpp.partition$training_index <= cl.len - temp.ddiff]
    temp.train_index <- c(temp.train_index, cpp.partition$training_index[cpp.partition$training_index > cl.len] - temp.ddiff)
    temp.test_index <- cpp.partition$test_index - temp.ddiff
    # print(paste("len3", length(temp.train_index), length(temp.test_index)))
    
    ## Store Combat processed data
    temp.data.ComBat <- ComBat_combine(batch = temp.cpp_source,
                                       label = temp.labels.combined,
                                       input_data = temp.data.raw)
    data.all[[pID]] <- temp.data.ComBat
    plabel.all[[pID]] <- temp.plot_label
    pID <- pID + 1
    
    ## Store SVA processed data
    temp.data.sva <- sva_combine(batch = temp.cpp_source,
                                 label = temp.labels.combined,
                                 input_data = temp.data.raw,
                                 n.sv = n.sv)
    data.all[[pID]] <- temp.data.sva
    plabel.all[[pID]] <- temp.plot_label
    pID <- pID + 1
    
    ## Store frozen SVA processed data
    temp.data.fsva <- fsva_combine(batch = temp.cpp_source,
                                   label = temp.labels.combined,
                                   input_data = temp.data.raw,
                                   training_ind = temp.train_index,
                                   test_ind = temp.test_index,
                                   n.sv = n.sv)
    #num_sv_method = "leek"
    data.all[[pID]] <- temp.data.fsva
    plabel.all[[pID]] <- temp.plot_label
    pID <- pID + 1
  }
  return(list(data.all, plabel.all))
}


format_PCAplot <- function (plt, plot_label, pOrder) {
  plot_label <<- plot_label
  plt <- plt + theme(legend.direction = 'horizontal', 
                     legend.position = 'top', plot.margin = unit(c(5.1,7,4.5,3.5)/2, "lines"), 
                     text = element_text(size=15), axis.title.x=element_text(vjust=-1.5)) + 
    scale_color_discrete(name = 'Dataset')
  plt <- plt + guides(colour = guide_legend(override.aes = list(size=2.5, linetype=0, shape = c(16, 17, 15))))
  plt <- plt + geom_point(aes(shape = factor(plot_label), colour = factor(plot_label)), show_guide = FALSE)
  plt <- plt + ggtitle(LETTERS[pOrder]) + theme(plot.title=element_text(hjust=-0.12, size = rel(1.75), face='bold'))
  return(eval(ggplotGrob(plt)))
}

#### Plot PCA of all given data
plot_processed_data_all <- function(data.all, plabel.all, name="drugX"){
  library(ggplot2)
  library(grid)
  library(RColorBrewer)
  require(gtable)
  require(gridExtra)
  
  lplots <- list() # list to keep all plots in
  
  for(pID in 1:length(data.all)) {
    ## Plot PCA of the data as is
    p1 <- show_pca(input_data = data.all[[pID]], label = plabel.all[[pID]], pca_line_plot = FALSE, print = FALSE)
    g1 <- format_PCAplot(p1, plabel.all[[pID]], pID)
    lplots <- c(lplots, list(g1))
  }
  
#   ## Add figure letters
#   for(pID in 1:length(lplots)) {
#     lplots[[pID]] <- lplots[[pID]] + ggtitle(LETTERS[pID]) + theme(plot.title=element_text(hjust=-0.12, size = rel(1.75), face='bold'))
#   }
  
  ## Extract Grobs
  # g1 <- ggplotGrob(p1)
  blank <- grid.rect(gp=gpar(col="white"))
  row1 <- list(blank, lplots[[1]], blank)
  gs <- lplots[2:length(lplots)]
  gs <- c(row1, gs)
  
  print("plotting")
  
  pdf(paste(name, "_all_t4.pdf", sep=""), width = 20, height = 25)#, onefile=FALSE)
  #grid.arrange(blank, g1, blank, g2, g3, g4, nrow=2, ncol=3)
  grid.arrange(grobs=gs, nrow=4, ncol=3)
  dev.off()

}

# !!!TODO!!!!: check ComBat function and decide of preprocess CL data with SVA beforehand

#processed.data <- compute_preprocessing(bortezomib, bortezomib.labels, sampleinfo.cgp, n.sv = 2, preprocessCLsva = TRUE)
#plot_processed_data_all(processed.data[[1]], processed.data[[2]], name="bortezomib")

# plot_processed_data_all(docetaxel, docetaxel.labels, name="docetaxel", n.sv = 2)
# plot_processed_data_all(erlotinib, erlotinib.labels, name="erlotinib", n.sv = 2)
# plot_processed_data_all(epirubicin, epirubicin.labels, name="epirubicin", n.sv = 2)
