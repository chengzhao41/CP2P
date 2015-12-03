rootdir <- "~/Dropbox/CP2P"
setwd(rootdir)

source("Common/preparing_data_helper.R")
source("Common/comGENE.R")

load("Bortezomib/WS/bortezomib_data.RData"); bortezomib$patient <- bortezomib$patient.combat
load("Docetaxel/WS/docetaxel_data.RData")
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

plot_processed_data_all <- function(data, label, name="drugX", n.sv = 0){
  library(ggplot2)
  library(grid)
  library(RColorBrewer)
  require(gtable)
  require(gridExtra)
  
  ## Split to training set and patient-only test set
  labels_cell_lines <- label$AUC
  cl.len <- length(label$AUC)
  labels_patient <- label$patient
  training_amount <- round(length(labels_patient) * 0.7) #30% test set
  cpp.partition <- generate_single_random_partition.cpp_var(labels_cell_lines, labels_patient, training_amount)
  #print(length(cpp.partition$test_index))
  #print(length(cpp.partition$training_index))
  
  ## Pool AUC, IC50, Slope data and labels together
  #AUC
  temp.data <- comGENE(scale(data$cgp_AUC), scale(data$patient))
  data.raw <- list(rbind(temp.data[[1]], temp.data[[2]]))
  labels.combined <- list(c(label$AUC, label$patient))
  #IC50
  temp.data <- comGENE(scale(data$cgp_IC50), scale(data$patient))
  data.raw[[2]] <- rbind(temp.data[[1]], temp.data[[2]])
  labels.combined[[2]] <- c(label$IC50, label$patient)
  #Slope
  temp.data <- comGENE(scale(data$cgp_slope), scale(data$patient))
  data.raw[[3]] <- rbind(temp.data[[1]], temp.data[[2]])
  labels.combined[[3]] <- c(label$slope, label$patient)
  
  
  ## Get AUC data and labels
  temp.data.raw <- data.raw[[1]]
  temp.labels.combined <- labels.combined[[1]]
  
  ## Cell lines, train patient and test patient labels for PCA plot
  temp.plot_label <- cpp.partition$cpp_flags
  
  pID <- 1
  ## Plot PCA of the data as is
  p1 <- show_pca(input_data = temp.data.raw, label = temp.plot_label, pca_line_plot = FALSE, print = FALSE)
  g1 <- format_PCAplot(p1, temp.plot_label, pID)
  pID <- pID + 1
  lplots <- list(g1) # list to keep all plots in
  
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
    
    ## Plot PCA of Combat processed data
    temp.data.ComBat <- ComBat_combine(batch = temp.cpp_source,
                                      label = temp.labels.combined,
                                      input_data = temp.data.raw)
    
    p2 <- show_pca(input_data = temp.data.ComBat, label = temp.plot_label, pca_line_plot = FALSE, print = FALSE)
    g2 <- format_PCAplot(p2, temp.plot_label, pID)
    pID <- pID + 1
    lplots <- c(lplots, list(g2))
    
    ## Plot PCA of SVA processed data
    temp.data.sva <- sva_combine(batch = temp.cpp_source,
                                 label = temp.labels.combined,
                                 input_data = temp.data.raw,
                                 n.sv = n.sv)
    
    p3 <- show_pca(input_data = temp.data.sva, label = temp.plot_label, pca_line_plot = FALSE, print = FALSE)
    g3 <- format_PCAplot(p3, temp.plot_label, pID)
    pID <- pID + 1
    lplots <- c(lplots, list(g3))
    
    ## Plot PCA of frozen SVA processed data
    temp.data.fsva <- fsva_combine(batch = temp.cpp_source,
                                 label = temp.labels.combined,
                                 input_data = temp.data.raw,
                                 training_ind = temp.train_index,
                                 test_ind = temp.test_index,
                                 n.sv = n.sv)
                                 #num_sv_method = "leek"
    
    p4 <- show_pca(input_data = temp.data.fsva, label = temp.plot_label, pca_line_plot = FALSE, print = FALSE)
    g4 <- format_PCAplot(p4, temp.plot_label, pID)
    pID <- pID + 1
    lplots <- c(lplots, list(g4))
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
  
  pdf(paste(name, "_all_sva2.pdf", sep=""), width = 20, height = 25)#, onefile=FALSE)
  #grid.arrange(blank, g1, blank, g2, g3, g4, nrow=2, ncol=3)
  grid.arrange(grobs=gs, nrow=4, ncol=3)
  dev.off()

  }

plot_processed_data_all(bortezomib, bortezomib.labels, name="bortezomib", n.sv = 2)
plot_processed_data_all(docetaxel, docetaxel.labels, name="docetaxel", n.sv = 2)
# plot_processed_data_all(erlotinib, erlotinib.labels, name="erlotinib", n.sv = 2)
# plot_processed_data_all(epirubicin, epirubicin.labels, name="epirubicin", n.sv = 2)
