## Main script to compute and plot PCA of various data preprocessing methods

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

  # flags to denote: cell lines in training; patients in training; patients in test
  cell.line.name = "GDSC "
  if (length(input_labels_cell_lines) == 38) { # hack to change name in the legened in case of Epirubicin
    cell.line.name = "GRAY "
  }
  if (length(input_labels_cell_lines) == 489) { # hack to change name in the legened in case of Erlotinib CCLE
    cell.line.name = "CCLE "
  }

  cpp_flags <- rep(cell.line.name, length(input_combined_labels))
  cpp_flags[temp.patients_in_training] <- "Patients  70% " # a hack to make this be before the "patients 30%"
  cpp_flags[temp.test_index] <- "Patients 30%"
  
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


compute_preprocessing <- function(data, label, sampleinfo, n.sv, preprocessCLsva = 0, baselinePlot = c('AUC', 'IC50', 'Slope')){
  ## Split to training set and patient-only test set
  labels_cell_lines <- label$AUC
  cl.len <- length(label$AUC)
  labels_patient <- label$patient
  training_amount <- round(length(labels_patient) * 0.7) #30% test set
  cpp.partition <- generate_single_random_partition.cpp_var(labels_cell_lines, labels_patient, training_amount)
  #print(length(cpp.partition$test_index))
  #print(length(cpp.partition$training_index))
  
  print("preprocessing cell lines only")
  data.raw <- list()
  labels.combined <- list()
  ## Pool AUC, IC50, Slope data and labels together, preprocess cell lines with SVA first
  if (preprocessCLsva != 0) {
    cl_AUC.sva <- sva_combine(batch = sampleinfo$tissue.type[label$AUC_ind],  
                                 label = label$AUC, input_data = scale(data$cl_AUC), n.sv=preprocessCLsva)
    temp.data <- comGENE(scale(cl_AUC.sva), scale(data$patient))
  } else {
    temp.data <- comGENE(scale(data$cl_AUC), scale(data$patient))
  }
  data.raw[[1]] <- rbind(temp.data[[1]], temp.data[[2]])
  labels.combined[[1]] <- c(label$AUC, label$patient)
  #IC50
  if (preprocessCLsva != 0) {
    cl_IC50.sva <- sva_combine(batch = sampleinfo$tissue.type[label$IC50_ind],  
                                label = label$IC50, input_data = scale(data$cl_IC50), n.sv=preprocessCLsva)
    temp.data <- comGENE(scale(cl_IC50.sva), scale(data$patient))
  } else {
    temp.data <- comGENE(scale(data$cl_IC50), scale(data$patient))
  }
  data.raw[[2]] <- rbind(temp.data[[1]], temp.data[[2]])
  labels.combined[[2]] <- c(label$IC50, label$patient)
  #Slope
  if (preprocessCLsva != 0) {
    cl_slope.sva <- sva_combine(batch = sampleinfo$tissue.type[label$slope_ind],  
                                 label = label$slope, input_data = scale(data$cl_slope), n.sv=preprocessCLsva)
    temp.data <- comGENE(scale(cl_slope.sva), scale(data$patient))
  } else {
    temp.data <- comGENE(scale(data$cl_slope), scale(data$patient))
  }
  data.raw[[3]] <- rbind(temp.data[[1]], temp.data[[2]])
  labels.combined[[3]] <- c(label$slope, label$patient)
  
  pID <- 1
  data.all <- list()
  plabel.all <- list()
  ptitle.all <- list()
  ## Store raw data as is
  if (baselinePlot == 'AUC'){
    temp.data <- comGENE(scale(data$cl_AUC), scale(data$patient))
  } else if (baselinePlot == 'IC50'){
    temp.data <- comGENE(scale(data$cl_IC50), scale(data$patient))
  } else if (baselinePlot == 'slope'){
    temp.data <- comGENE(scale(data$cl_slope), scale(data$patient))
  }
  temp.data.raw <- rbind(temp.data[[1]], temp.data[[2]])
  temp.label.combined <- c(label$slope, label$patient)
  
  ## Compensate for different number of cell lines with AUC/IC50/Slope labels
  temp.ddiff <- length(labels.combined[[1]]) - length(temp.label.combined)
  print(paste("adjusting", temp.ddiff))
  temp.plot_label <- cpp.partition$cpp_flags[(temp.ddiff+1) : length(cpp.partition$cpp_flags)]
  
  data.all[[pID]] <- temp.data.raw
  plabel.all[[pID]] <- temp.plot_label
  ptitle.all[[pID]] <- paste(baselinePlot, "labels before cell lines vs patients homogenization")
  pID <- pID + 1
  
  print("processing CL + patients")
  for(outcomeID in 1:1) {
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
    ptitle.all[[pID]] <- paste(c('AUC', 'IC50', 'Slope')[outcomeID], "labels after ComBat homogenization of entire data set")
    pID <- pID + 1
    
    ## Store SVA processed data
    temp.data.sva <- sva_combine(batch = temp.cpp_source,
                                 label = temp.labels.combined,
                                 input_data = temp.data.raw,
                                 n.sv = n.sv)
    data.all[[pID]] <- temp.data.sva
    plabel.all[[pID]] <- temp.plot_label
    ptitle.all[[pID]] <- paste(c('AUC', 'IC50', 'Slope')[outcomeID], "labels after SVA homogenization of entire data set")
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
    ptitle.all[[pID]] <- paste(c('AUC', 'IC50', 'Slope')[outcomeID], "labels after frozen SVA homogenization")
    pID <- pID + 1
  }
  return(list(data.all, plabel.all, ptitle.all))
}

format_PCAplot <- function (plt, plot_label, pOrder) {
  # The color-blind friendly palette with grey
  # cbPalette <- c("#B3B3B3", "#DE9900", "#48A6DB", "#00996F", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  cbPalette <- c("#B3B3B3", "#0D5DE0", "#E6374F")
  
  plot_label <<- plot_label
  plt <- plt + theme(legend.direction = 'horizontal', 
                     legend.position = 'top', plot.margin = unit(c(5.1,7,4.5,3.5)/2, "lines"), 
                     text = element_text(size=15), axis.title.x=element_text(vjust=-1.5))
  plt <- plt + guides(colour = guide_legend(override.aes = list(size=2.5, linetype=0, shape = c(16, 17, 15))))
  plt <- plt + geom_point(aes(shape = factor(plot_label), colour = factor(plot_label)), show.legend = FALSE)
  # plt <- plt + ggtitle(LETTERS[pOrder]) + theme(plot.title=element_text(hjust=-0.12, size = rel(1.75), face='bold'))
  plt <- plt + scale_colour_manual(name = 'Dataset: ', values=cbPalette)
  return(eval(ggplotGrob(plt)))
}

#### Plot PCA of all given data
plot_processed_data_all <- function(data.all, plabel.all, ptitle.all, name="drugX"){
  library(ggplot2)
  library(grid)
  library(RColorBrewer)
  require(gtable)
  require(gridExtra)
  
  print("running PCA of the all data for plotting")
  lplots <- list() # list to keep all plots in
  
  for(pID in 1:length(data.all)) { #length(data.all)
    ## Plot PCA of the data as is
    p1 <- show_pca(input_data = data.all[[pID]], label = plabel.all[[pID]], pca_line_plot = FALSE, print = FALSE, useColorBlindScheme = TRUE)
    p1 <- p1 + ggtitle(ptitle.all[[pID]]) + theme(plot.title=element_text(hjust=0.5, size = rel(1.15)))
    g1 <- format_PCAplot(p1, plabel.all[[pID]], pID)
    g2 <- arrangeGrob(g1, top = textGrob(LETTERS[pID+3], x = unit(0.05, "npc"), y = unit(0.35, "npc"),
                                     just = c("left", "top"),
                                     gp = gpar(col="black", fontsize=22, fontface="bold")))
    lplots <- c(lplots, list(g2))
  }
  
  ## Extract Grobs
  # g1 <- ggplotGrob(p1)
  blank <- grid.rect(gp=gpar(col="white"))
  row1 <- list(blank, lplots[[1]], blank)
  gs <- lplots[2:length(lplots)]
  gs <- c(row1, gs)
  
  print("plotting")
  
  pdf(paste0(name, "_all.pdf"), width = 20, height = 25)#, onefile=FALSE)
  #grid.arrange(blank, g1, blank, g2, g3, g4, nrow=2, ncol=3)
  grid.arrange(grobs=gs, nrow=4, ncol=3)
  dev.off()

}

######################  MAIN  #######################

args <- commandArgs(TRUE)

ENV <- args[1]
stopifnot(ENV %in% c('dp', 'cb'))

DRUG2RUN <- args[2]
stopifnot(DRUG2RUN %in% c('bor', 'doc', 'epir', 'erl-gdsc', 'erl-ccle'))

if (ENV == 'dp'){
  rootdir <- "~/Dropbox/CP2P"
  setwd(rootdir)

  source("Common/preparing_data_helper.R")
  source("Common/comGENE.R")
  rdata_prefix_bor <- "../Bortezomib/WS/"
  rdata_prefix_doc <- "../Docetaxel/WS/"
  rdata_prefix_epir <- "../Epirubicin/WS/"
  rdata_prefix_erl <- "../Erlotinib/WS/"
}
if (ENV == 'cb'){
  rootdir <- "/dupa-filer/laci/CP2P/"
  setwd(rootdir)

  source("SVA_plots/preparing_data_helper.R")
  source("SVA_plots/comGENE.R")
  rdata_prefix_bor <- ""
  rdata_prefix_doc <- ""
  rdata_prefix_epir <- ""
  rdata_prefix_erl <- ""
}

#### Bortezomib
if (DRUG2RUN == 'bor'){
  load(paste0(rdata_prefix_bor, "bortezomib_data.RData")); bortezomib$patient <- bortezomib$patient.combat
  bortezomib$cl_IC50 <- bortezomib$gdsc_IC50
  bortezomib$cl_AUC <- bortezomib$gdsc_AUC
  bortezomib$cl_slope <- bortezomib$gdsc_slope
  bor.processed.data <- compute_preprocessing(bortezomib, bortezomib.labels, sampleinfo.gdsc, n.sv = 2, preprocessCLsva = 3, baselinePlot = 'Slope')
  save(bor.processed.data, file = "bor.processed.data.RData")
  plot_processed_data_all(bor.processed.data[[1]], bor.processed.data[[2]], bor.processed.data[[3]], name="bortezomib_2sva")
  }

#### Docetaxel
if (DRUG2RUN == 'doc'){
  load(paste0(rdata_prefix_doc, "docetaxel_data.RData"))
  docetaxel$cl_IC50 <- docetaxel$gdsc_IC50
  docetaxel$cl_AUC <- docetaxel$gdsc_AUC
  docetaxel$cl_slope <- docetaxel$gdsc_slope
  doc.processed.data <- compute_preprocessing(docetaxel, docetaxel.labels, sampleinfo.gdsc, n.sv = 2, preprocessCLsva = 0, baselinePlot = 'Slope')
  save(doc.processed.data, file = "doc.processed.data.RData")
  plot_processed_data_all(doc.processed.data[[1]], doc.processed.data[[2]], doc.processed.data[[3]], name="docetaxel_2sva")
}

#### Epirubicin
if (DRUG2RUN == 'epir'){
  load(paste0(rdata_prefix_epir, "epirubicin_data.RData"))
  epirubicin$cl_IC50 <- epirubicin$gray_IC50
  epirubicin$cl_AUC <- epirubicin$gray_AUC
  epirubicin$cl_slope <- epirubicin$gray_slope
  epir.processed.data <- compute_preprocessing(epirubicin, epirubicin.labels, sampleinfo.gray, n.sv = 3, preprocessCLsva = 0, baselinePlot = 'Slope')
  save(epir.processed.data, file = "epir.processed.data.RData")
  # load("epir.processed.data.RData")
  plot_processed_data_all(epir.processed.data[[1]], epir.processed.data[[2]], epir.processed.data[[3]], name="epirubicin_3sva")
}

#### Erlotinib
if (DRUG2RUN == 'erl-gdsc'){
  load(paste0(rdata_prefix_erl, "erlotinib_homogenized_data_gdsc.RData"))
  erlotinib$cl_IC50 <- erlotinib$gdsc_IC50
  erlotinib$cl_AUC <- erlotinib$gdsc_AUC
  erlotinib$cl_slope <- erlotinib$gdsc_slope
  colnames(erlotinib$cl_AUC) <- colnames(erlotinib$cl_IC50) # fix gene names
  erl.processed.data <- compute_preprocessing(erlotinib, erlotinib.labels, sampleinfo.gdsc, n.sv = 2, preprocessCLsva = 0, baselinePlot = 'Slope')
  save(erl.processed.data, file = "erl.gdsc.processed.data.RData")
  # load("erl.gdsc.processed.data.RData")
  plot_processed_data_all(erl.processed.data[[1]], erl.processed.data[[2]], erl.processed.data[[3]], name="erlotinibGDSC_2sva")
}
if (DRUG2RUN == 'erl-ccle'){
  load(paste0(rdata_prefix_erl, "erlotinib_homogenized_data_ccle.RData"))
  erlotinib$cl_IC50 <- erlotinib$ccle_IC50
  erlotinib$cl_AUC <- erlotinib$ccle_AUC
  erlotinib$cl_slope <- erlotinib$ccle_slope
  colnames(erlotinib$cl_AUC) <- colnames(erlotinib$cl_IC50) # fix gene names
  erl.processed.data <- compute_preprocessing(erlotinib, erlotinib.labels, sampleinfo.ccle, n.sv = 2, preprocessCLsva = 0, baselinePlot = 'Slope')
  save(erl.processed.data, file = "erl.ccle.processed.data.RData")
  # load("erl.ccle.processed.data.RData")
  plot_processed_data_all(erl.processed.data[[1]], erl.processed.data[[2]], erl.processed.data[[3]], name="erlotinibCCLE_2sva")
}
