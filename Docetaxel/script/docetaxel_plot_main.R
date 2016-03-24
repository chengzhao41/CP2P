# load library ------------------------------------------------------------
library("doParallel")
library("RColorBrewer")
source('Common/plot_data_helper.R')

output_dir <- "Docetaxel/plots/"

# Docetaxel p2p ----------------------------------------------------------
rm(list=ls(pattern="input.*"))
input.seq <- seq(1, 10, 1)
input.num_partitions <- rep(list(1:100), 9)
input.num_partitions[[10]] <- 1:24
input.x_axis <- seq(14, 23, 1)
input.WS_path <- "Docetaxel/computed_WS/Var P/docetaxel_p2p_1to100_parInd"
input.WS_paths <- rep(input.WS_path, 9)
input.WS_paths[10] <- "Docetaxel/computed_WS/Var P/docetaxel_p2p_1to24_parInd"
input.output_file_name_plot <- "doc_p2p.png"
input.output_file_name_WS <- "doc_p2p.RData"
input.xlab <- "# of Patient Samples in Training Data"
# create WS for plot ------------------------------------------------------
create_plot_WS(input.seq, input.num_partitions, input.x_axis, 
               input.WS_paths, input.output_file_name_WS, output_dir, input.type_measure = "acc")
# create plot -------------------------------------------------------------
load(paste0(output_dir, input.output_file_name_WS))
create_plot(input.output_file_name_plot, varying_training_matrix, input.xlab, output_dir)

# Docetaxel cp2p var p ------------------------------------------------------
rm(list=ls(pattern="input.*"))
input.seq <- seq(1, 10, 1)
input.num_partitions <- rep(list(1:100), 9)
input.num_partitions[[10]] <- 1:24
input.x_axis <- seq(14, 23, 1)
input.xlab <- "# of Patient Samples in Training Data"
output_dir <- "Docetaxel/plots/"

input.response_labels <- c("ic50", "auc", "slope", "ic50_breast", "auc_breast", "slope_breast")
for (i in 1:length(input.response_labels)) {
  input.response_label <- input.response_labels[i]
  # based on response_label
  input.WS_path <- paste0("Docetaxel/computed_WS/Var P/docetaxel_cp2p_", input.response_label, "_1to100_parInd")
  input.WS_paths <- rep(input.WS_path, 9)
  input.WS_paths[10] <- paste0("Docetaxel/computed_WS/Var P/docetaxel_cp2p_", input.response_label, "_1to24_parInd")
  input.output_file_name_plot <- paste0("doc_cp2p_", input.response_label, ".png")
  input.output_file_name_WS <- paste0("doc_cp2p_", input.response_label, ".RData")
  # create WS for plot ------------------------------------------------------
  create_plot_WS(input.seq, input.num_partitions, input.x_axis, 
                 input.WS_paths, input.output_file_name_WS, output_dir, input.type_measure = "acc")
  # create plot -------------------------------------------------------------
  load(paste0(output_dir, input.output_file_name_WS))
  create_plot(input.output_file_name_plot, varying_training_matrix, input.xlab, output_dir)
}

# Docetaxel cp2p var c for p23 ---------------------------------------------------
rm(list=ls(pattern="input.*"))
input.seq <- seq(1, 20, 1)
input.num_partitions <- list(1:24)
input.x_axis <- seq(from = 30, to = 605, by = 30)
input.xlab <- "# of Cell Line Samples in Training Data"
output_dir <- "Docetaxel/plots/"

# based on response_label
input.response_labels <- c("cp2p_ic50_p23", "cp2p_auc_p23", "cp2p_slope_p23",
                           "c2p_ic50_p23", "c2p_auc_p23", "c2p_slope_p23")
for (i in 1:length(input.response_labels)) {
  input.response_label <- input.response_labels[i]
  # based on response_label
  input.WS_paths <- paste0("Docetaxel/computed_WS/Var C p23/docetaxel_", input.response_label, "_1to24_parInd")
  input.output_file_name_plot <- paste0("doc_", input.response_label, ".png")
  input.output_file_name_WS <- paste0("doc_", input.response_label, ".RData")
  # create WS for plot ------------------------------------------------------
  create_plot_WS(input.seq, input.num_partitions, input.x_axis, 
                 input.WS_paths, input.output_file_name_WS, output_dir, input.type_measure = "acc")
  # create plot -------------------------------------------------------------
  load(paste0(output_dir, input.output_file_name_WS))
  create_plot(input.output_file_name_plot, varying_training_matrix, input.xlab, output_dir)
}

# Docetaxel best var p plot ----------------------------------------------
load(paste0(output_dir, "doc_p2p.RData"))
p2p.results <- varying_training_matrix
p2p.x_axis <- input.x_axis 
rm("varying_training_matrix")
load(paste0(output_dir, "doc_cp2p_slope.RData"))
cp2p.results <- varying_training_matrix
rm("varying_training_matrix")
load(paste0(output_dir, "doc_c2p_slope_p23.RData"))
c2p.results <- varying_training_matrix
rm("varying_training_matrix")
load(paste0(output_dir, "doc_cp2p_slope_breast.RData"))
cp2p_breast.results <- varying_training_matrix
rm("varying_training_matrix")

cols <- brewer.pal(n = 4, name = 'Dark2')
input.color_manual <- c("P2P" = cols[1], "C2P" = cols[2], "CP2P" = cols[3], "CP2P Breast" = cols[4])
input.shape_manual <- c("P2P" = 16, "C2P" = NA, "CP2P" = 17, "CP2P Breast" = 18)

create_plot_best(input.x_axis = p2p.x_axis, 
                 input_x_labels = "# of Patient Samples in Training",
                 input.results = list(p2p.results, cp2p.results, cp2p_breast.results), 
                 input.results_constant = list(c2p.results),
                 input.output_file_name_plot = "docetaxel_var_p_best.png", 
                 output_dir = output_dir, 
                 input.labels = c("P2P", "CP2P", "CP2P Breast"),
                 input.labels_constant = c("C2P"),
                 include_std = FALSE,
                 input.color_manual,
                 input.shape_manual,
                 input.constant_best_ind = dim(c2p.results$mean)[1])


# Docetaxel best var c p23 plot ----------------------------------------------
rm(cp2p.results)
load(paste0(output_dir, "doc_p2p.RData"))
p2p.results <- varying_training_matrix
rm("varying_training_matrix")
load(paste0(output_dir, "doc_cp2p_slope_p23.RData"))
cp2p.results <- varying_training_matrix
rm("varying_training_matrix")
load(paste0(output_dir, "doc_c2p_slope_p23.RData"))
c2p.results <- varying_training_matrix
c2p.x_axis <- input.x_axis 
rm("varying_training_matrix")

cols <- brewer.pal(n = 4, name = 'Dark2')
input.color_manual <- c("P2P" = cols[1], "C2P" = cols[2], "CP2P" = cols[3])
input.shape_manual <- c("P2P" = NA, "C2P" = 16, "CP2P" = 17)

create_plot_best(input.x_axis = c2p.x_axis, 
                 input_x_labels = "# of Cell Line Samples in Training",
                 input.results = list(c2p.results, cp2p.results), 
                 input.results_constant = list(p2p.results),
                 input.output_file_name_plot = "docetaxel_var_c_p23_best.png", 
                 output_dir = output_dir, 
                 input.labels = c("C2P", "CP2P"),
                 input.labels_constant = c("P2P"),
                 include_std = FALSE,
                 input.color_manual,
                 input.shape_manual,
                 input.constant_best_ind = dim(p2p.results$mean)[1])
