# load library ------------------------------------------------------------
library("doParallel")
library("RColorBrewer")
source('Common/plot_data_helper.R')

output_dir <- "Epirubicin/plots/"

# Epirubicin p2p ----------------------------------------------------------
rm(list=ls(pattern="input.*"))
input.seq <- seq(2, 9, 1)
input.num_partitions <- rep(list(1:100), length(input.seq))
input.x_axis <- seq(30, 100, 10)
stopifnot(length(input.seq) == length(input.x_axis))
input.WS_path <- "Epirubicin/computed_WS/Var P/epirubicin_p2p_1to100_parInd"
input.output_file_name_plot <- "epi_p2p.png"
input.output_file_name_WS <- "epi_p2p.RData"
input.xlab <- "# of Patient Samples in Training Data"
# create WS for plot ------------------------------------------------------
create_plot_WS(input.seq, input.num_partitions, input.x_axis, 
               input.WS_path, input.output_file_name_WS, output_dir, input.type_measure = "auc")
# create plot -------------------------------------------------------------
load(paste0(output_dir, input.output_file_name_WS))
create_plot(input.output_file_name_plot, varying_training_matrix, input.xlab, output_dir)

# Epirubicin cp2p var p ------------------------------------------------------
rm(list=ls(pattern="input.*"))
input.seq <- seq(2, 9, 1)
input.num_partitions <- rep(list(1:100), length(input.seq))
input.x_axis <- seq(30, 100, 10)

stopifnot(length(input.seq) == length(input.x_axis))
input.xlab <- "# of Patient Samples in Training Data"
output_dir <- "Epirubicin/plots/"

input.response_labels <- c("auc", "slope")
for (i in 1:length(input.response_labels)) {
  input.response_label <- input.response_labels[i]
  input.WS_path <- paste0("Epirubicin/computed_WS/Var P/epirubicin_cp2p_", input.response_label, "_1to100_parInd")
  input.output_file_name_plot <- paste0("epi_cp2p_", input.response_label, ".png")
  input.output_file_name_WS <- paste0("epi_cp2p_", input.response_label, ".RData")
  # create WS for plot ------------------------------------------------------
  create_plot_WS(input.seq, input.num_partitions, input.x_axis, 
                 input.WS_path, input.output_file_name_WS, output_dir, input.type_measure = "auc")
  # create plot -------------------------------------------------------------
  load(paste0(output_dir, input.output_file_name_WS))
  create_plot(input.output_file_name_plot, varying_training_matrix, input.xlab, output_dir)
}

# Epirubicin cp2p var c for p50 ---------------------------------------------------
rm(list=ls(pattern="input.*"))
input.seq <- seq(2, 7, 1)
input.num_partitions <- rep(list(1:100), length(input.seq))
input.x_axis <- seq(from = 10, to = 38, by = 5)

stopifnot(length(input.seq) == length(input.x_axis))
input.xlab <- "# of Cell Line Samples in Training Data"
output_dir <- "Epirubicin/plots/"

input.response_labels <- c("cp2p_auc_p50","cp2p_slope_p50")

for (i in 1:length(input.response_labels)) {
  input.response_label <- input.response_labels[i]
  input.WS_path <- paste0("Epirubicin/computed_WS/var c 50/epirubicin_", input.response_label, "_1to100_parInd")
  input.output_file_name_plot <- paste0("epi_", input.response_label, ".png")
  input.output_file_name_WS <- paste0("epi_", input.response_label, ".RData")
  # create WS for plot ------------------------------------------------------
  create_plot_WS(input.seq, input.num_partitions, input.x_axis, 
                 input.WS_path, input.output_file_name_WS, output_dir, input.type_measure = "acc")
  # create plot -------------------------------------------------------------
  load(paste0(output_dir, input.output_file_name_WS))
  create_plot(input.output_file_name_plot, varying_training_matrix, input.xlab, output_dir)
}

# Epirubicin c2p var c for p50 --------------------------------------------
rm(list=ls(pattern="input.*"))
input.seq <- seq(1, 1, 1)
input.num_partitions <- rep(list(1:1), length(input.seq))
input.x_axis <- seq(from = 38, to = 38, by = 1)

stopifnot(length(input.seq) == length(input.x_axis))
input.xlab <- "# of Cell Line Samples in Training Data"
output_dir <- "Epirubicin/plots/"

input.response_labels <- c("c2p_auc_p50", "c2p_slope_p50")

for (i in 1:length(input.response_labels)) {
  input.response_label <- input.response_labels[i]
  input.WS_path <- paste0("Epirubicin/computed_WS/var c 50/epirubicin_", input.response_label, "_1to1_parInd")
  input.output_file_name_plot <- paste0("epi_", input.response_label, ".png")
  input.output_file_name_WS <- paste0("epi_", input.response_label, ".RData")
  # create WS for plot ------------------------------------------------------
  create_plot_WS(input.seq, input.num_partitions, input.x_axis, 
                 input.WS_path, input.output_file_name_WS, output_dir, input.type_measure = "auc")
}
warning("created WS needs some modification in terms of the format.")

# Epirubicin best var p plot ----------------------------------------------
load(paste0(output_dir, "epi_p2p.RData"))
p2p.results <- varying_training_matrix
p2p.x_axis <- input.x_axis 
rm("varying_training_matrix")
load(paste0(output_dir, "epi_cp2p_slope.RData"))
cp2p.results <- varying_training_matrix
rm("varying_training_matrix")
load(paste0(output_dir, "epi_c2p_slope_p50.RData"))
c2p.results <- varying_training_matrix
rm("varying_training_matrix")

cols <- brewer.pal(n = 4, name = 'Dark2')
input.color_manual <- c("P2P" = cols[1], "C2P" = cols[2], "CP2P" = cols[3])
input.shape_manual <- c("P2P" = 16, "C2P" = NA, "CP2P" = 17)

best <- create_plot_best(input.x_axis = p2p.x_axis, 
                 input_x_labels = "# of Patient Samples in Training",
                 input.results = list(p2p.results, cp2p.results), 
                 input.results_constant = list(c2p.results),
                 input.output_file_name_plot = "epirubicin_var_p_best.png", 
                 output_dir = output_dir, 
                 input.labels = c("P2P", "CP2P"),
                 input.labels_constant = ("C2P"),
                 include_std = TRUE,
                 input.color_manual,
                 input.shape_manual,
                 input.constant_best_ind = 1)

temp.best <- c(best$best_model, names(best$best_model_constant))
for (i in 1:3) {
  create_plot_best(input.x_axis = p2p.x_axis, 
                   input_x_labels = "# of Patient Samples in Training",
                   input.results = list(p2p.results, cp2p.results), 
                   input.results_constant = list(c2p.results),
                   input.output_file_name_plot = paste0("epirubicin_var_p_best_", temp.best[i],".png"), 
                   output_dir = output_dir, 
                   input.labels = c("P2P", "CP2P"),
                   input.labels_constant = c("C2P"),
                   include_std = TRUE,
                   input.color_manual,
                   input.shape_manual,
                   input.constant_best_ind = dim(c2p.results$mean)[1],
                   input.best_model_ind = list(which(colnames(p2p.results$mean) == temp.best[i]),
                                               which(colnames(cp2p.results$mean) == temp.best[i])),
                   input.constant_best_model_ind = list(which(colnames(c2p.results$mean) == temp.best[i])))
}

# Epirubicin best var c p50 plot ----------------------------------------------
load(paste0(output_dir, "epi_p2p.RData"))
p2p.results <- varying_training_matrix
p2p.x_axis <- input.x_axis 
rm("varying_training_matrix")
load(paste0(output_dir, "epi_cp2p_slope_p50.RData"))
cp2p.results <- varying_training_matrix
cp2p.x_axis <- input.x_axis
rm("varying_training_matrix")
load(paste0(output_dir, "epi_c2p_slope_p50.RData"))
c2p.results <- varying_training_matrix
rm("varying_training_matrix")
input.output_file_name_plot = "epirubicin_var_c_p50_best.png"

cols <- brewer.pal(n = 4, name = 'Dark2')
input.color_manual <- c("P2P" = cols[1], "CP2P" = cols[3])
input.shape_manual <- c("P2P" = NA, "CP2P" = 17)

create_plot_best(input.x_axis = cp2p.x_axis, 
                 input_x_labels = "# of Cell Line Samples in Training",
                 input.results = list(cp2p.results), 
                 input.results_constant = list(p2p.results),
                 input.output_file_name_plot = "epirubicin_var_c_p50_best.png", 
                 output_dir = output_dir, 
                 input.labels = c("CP2P"),
                 input.labels_constant = c("P2P"),
                 include_std = TRUE,
                 input.color_manual,
                 input.shape_manual,
                 input.constant_best_ind = which(p2p.x_axis == 50))
