# load library ------------------------------------------------------------
library("doParallel")
library("RColorBrewer")
source('Common/plot_data_helper.R')

output_dir <- "Erlotinib/plots/"

# Erlotinib p2p ----------------------------------------------------------
rm(list=ls(pattern="input.*"))
input.seq <- seq(1, 12, 1)
input.num_partitions <- rep(list(1:100), 11)
input.num_partitions[[12]] <- 1:25
input.x_axis <- seq(13, 24, 1)
input.WS_path <- "Erlotinib/computed_WS/p2p/erlotinib_gdsc_p2p_1to100_parInd"
input.WS_paths <- rep(input.WS_path, 11)
input.WS_paths[12] <- "Erlotinib/computed_WS/p2p/erlotinib_gdsc_p2p_1to25_parInd"
input.output_file_name_plot <- "erl_p2p.png"
input.output_file_name_WS <- "erl_p2p.RData"
input.xlab <- "# of Patient Samples in Training Data"
# create WS for plot ------------------------------------------------------
create_plot_WS(input.seq, input.num_partitions, input.x_axis, 
               input.WS_paths, input.output_file_name_WS, output_dir, input.type_measure = "acc")
# create plot -------------------------------------------------------------
load(paste0(output_dir, input.output_file_name_WS))
create_plot(input.output_file_name_plot, varying_training_matrix, input.xlab, output_dir)

# Erlotinib cp2p var p ------------------------------------------------------
rm(list=ls(pattern="input.*"))
#input.seq <- seq(1, 12, 1)
input.seq <- seq(1, 11, 1)
input.num_partitions <- rep(list(1:100), 11)
#input.num_partitions[[12]] <- 1:25
input.x_axis <- seq(13, 24, 1)
input.x_axis <- seq(13, 23, 1)
input.xlab <- "# of Patient Samples in Training Data"
output_dir <- "Erlotinib/plots/"

# input.response_labels <- c("ic50", "auc", "slope", 
#                            "ic50_lung", "auc_lung", "slope_lung")
input.response_labels <- c("ic50", "auc", "slope")
input.response_labels <- c("ic50_lung", "auc_lung", "slope_lung")
for (i in 1:length(input.response_labels)) {
  input.response_label <- input.response_labels[i]
  # based on response_label
  input.WS_path <- paste0("Erlotinib/computed_WS/all/erlotinib_all_cp2p_", input.response_label, "_1to100_parInd")
  input.WS_paths <- rep(input.WS_path, 11)
  #input.WS_paths[12] <- paste0("Erlotinib/computed_WS/gdsc/var p/erlotinib_gdsc_cp2p_", input.response_label, "_1to25_parInd")
  input.output_file_name_plot <- paste0("erl_all_cp2p_", input.response_label, ".png")
  input.output_file_name_WS <- paste0("erl_all_cp2p_", input.response_label, ".RData")
  # create WS for plot ------------------------------------------------------
  create_plot_WS(input.seq, input.num_partitions, input.x_axis, 
                 input.WS_paths, input.output_file_name_WS, output_dir, input.type_measure = "auc")
  # create plot -------------------------------------------------------------
  load(paste0(output_dir, input.output_file_name_WS))
  create_plot(input.output_file_name_plot, varying_training_matrix, input.xlab, output_dir)
}

# Erlotinib cp2p var c for p24 ---------------------------------------------------
rm(list=ls(pattern="input.*"))
input.seq <- seq(1, 13, 1)
input.num_partitions <- list(1:25)
input.x_axis <- seq(from = 30, to = 287, by = 20)
input.xlab <- "# of Cell Line Samples in Training Data"
output_dir <- "Erlotinib/plots/"

# based on response_label
# input.response_labels <- c("cp2p_ic50_p24", "cp2p_auc_p24", "cp2p_slope_p24",
#                            "c2p_ic50_p24", "c2p_auc_p24", "c2p_slope_p24")
input.response_labels <- c("cp2p_auc_p24", "cp2p_slope_p24",
                           "c2p_auc_p24", "c2p_slope_p24")

for (i in 1:length(input.response_labels)) {
  input.response_label <- input.response_labels[i]
  # based on response_label
  input.WS_paths <- paste0("Erlotinib/computed_WS/gdsc/Var C p24/erlotinib_gdsc_", input.response_label, "_1to25_parInd")
  input.output_file_name_plot <- paste0("erl_gdsc_", input.response_label, ".png")
  input.output_file_name_WS <- paste0("erl_gdsc_", input.response_label, ".RData")
  # create WS for plot ------------------------------------------------------
  create_plot_WS(input.seq, input.num_partitions, input.x_axis, 
                 input.WS_paths, input.output_file_name_WS, output_dir, input.type_measure = "auc")
  # create plot -------------------------------------------------------------
  load(paste0(output_dir, input.output_file_name_WS))
  create_plot(input.output_file_name_plot, varying_training_matrix, input.xlab, output_dir)
}

# Erlotinib best var p plot ----------------------------------------------
load(paste0(output_dir, "erl_p2p.RData"))
p2p.results <- varying_training_matrix
p2p.x_axis <- input.x_axis 
rm("varying_training_matrix")
load(paste0(output_dir, "erl_gdsc_cp2p_slope.RData"))
cp2p.results <- varying_training_matrix
rm("varying_training_matrix")
load(paste0(output_dir, "erl_gdsc_c2p_slope_p24.RData"))
c2p.results <- varying_training_matrix
rm("varying_training_matrix")
load(paste0(output_dir, "erl_gdsc_cp2p_slope_lung.RData"))
cp2p_lung.results <- varying_training_matrix
rm("varying_training_matrix")

cols <- brewer.pal(n = 4, name = 'Dark2')
input.color_manual <- c("P2P" = cols[1], "C2P" = cols[2], "CP2P" = cols[3], "CP2P Lung" = cols[4])
input.shape_manual <- c("P2P" = 16, "C2P" = NA, "CP2P" = 17, "CP2P Lung" = 18)

best <- create_plot_best(input.x_axis = p2p.x_axis, 
                         input_x_labels = "# of Patient Samples in Training",
                         input.results = list(p2p.results, cp2p.results, cp2p_lung.results), 
                         input.results_constant = list(c2p.results),
                         input.output_file_name_plot = "erlotinib_gdsc_var_p_best.png", 
                         output_dir = output_dir, 
                         input.labels = c("P2P", "CP2P", "CP2P Lung"),
                         input.labels_constant = c("C2P"),
                         include_std = FALSE,
                         input.color_manual,
                         input.shape_manual,
                         input.constant_best_ind = dim(c2p.results$mean)[1])

temp.best <- c(best$best_model, names(best$best_model_constant))
for (i in 1:4) {
  create_plot_best(input.x_axis = p2p.x_axis, 
                   input_x_labels = "# of Patient Samples in Training",
                   input.results = list(p2p.results, cp2p.results, cp2p_lung.results), 
                   input.results_constant = list(c2p.results),
                   input.output_file_name_plot = paste0("erlotinib_gdsc_var_p_best_", temp.best[i],".png"), 
                   output_dir = output_dir, 
                   input.labels = c("P2P", "CP2P", "CP2P Lung"),
                   input.labels_constant = c("C2P"),
                   include_std = FALSE,
                   input.color_manual,
                   input.shape_manual,
                   input.constant_best_ind = dim(c2p.results$mean)[1],
                   input.best_model_ind = list(which(colnames(p2p.results$mean) == temp.best[i]),
                                               which(colnames(cp2p.results$mean) == temp.best[i]),
                                               which(colnames(cp2p_lung.results$mean) == temp.best[i])),
                   input.constant_best_model_ind = list(which(colnames(c2p.results$mean) == temp.best[i])))
}


# Erlotinib best var c p24 plot ----------------------------------------------
rm(cp2p.results)
load(paste0(output_dir, "erl_p2p.RData"))
p2p.results <- varying_training_matrix
rm("varying_training_matrix")
load(paste0(output_dir, "erl_gdsc_cp2p_slope_p24.RData"))
cp2p.results <- varying_training_matrix
rm("varying_training_matrix")
load(paste0(output_dir, "erl_gdsc_c2p_slope_p24.RData"))
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
                 input.output_file_name_plot = "erlotinib_gdsc_var_c_p24_best.png", 
                 output_dir = output_dir, 
                 input.labels = c("C2P", "CP2P"),
                 input.labels_constant = c("P2P"),
                 include_std = FALSE,
                 input.color_manual,
                 input.shape_manual,
                 input.constant_best_ind = dim(p2p.results$mean)[1])
