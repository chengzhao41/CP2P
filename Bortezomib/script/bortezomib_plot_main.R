# load library ------------------------------------------------------------
library("doParallel")
library("RColorBrewer")
source('Common/plot_data_helper.R')

output_dir <- "Bortezomib/plots/"

# Bortezomib p2p ----------------------------------------------------------
rm(list=ls(pattern="input.*"))
input.seq <- seq(3, 15, 1)
input.num_partitions <- rep(list(1:100), length(input.seq))
input.x_axis <- seq(30, 150, 10)
stopifnot(length(input.seq) == length(input.x_axis))
input.WS_path <- "Bortezomib/computed_WS/No homo/bortezomib_p2p_1to100_parInd"
input.output_file_name_plot <- "bor_p2p.png"
input.output_file_name_WS <- "bor_p2p.RData"
input.xlab <- "# of Patient Samples in Training Data"
# create WS for plot ------------------------------------------------------
create_plot_WS(input.seq, input.num_partitions, input.x_axis, 
               input.WS_path, input.output_file_name_WS, output_dir, input.type_measure = "acc")
# create plot -------------------------------------------------------------
load(paste0(output_dir, input.output_file_name_WS))
create_plot(input.output_file_name_plot, varying_training_matrix, input.xlab, output_dir)

# Bortezomib cp2p var p ------------------------------------------------------
rm(list=ls(pattern="input.*"))
input.seq <- seq(2, 13, 1)
input.num_partitions <- rep(list(1:100), length(input.seq))
input.x_axis <- seq(20, 130, 10)
stopifnot(length(input.seq) == length(input.x_axis))
input.xlab <- "# of Patient Samples in Training Data"
output_dir <- "Bortezomib/plots/"

input.response_labels <- c("ic50", "auc", "slope")
for (i in 1:length(input.response_labels)) {
  input.response_label <- input.response_labels[i]
  input.WS_path <- paste0("Bortezomib/computed_WS/Exp/bortezomib_cp2p_", input.response_label, "_1to100_parInd")
  input.output_file_name_plot <- paste0("bor_cp2p_", input.response_label, ".png")
  input.output_file_name_WS <- paste0("bor_cp2p_", input.response_label, ".RData")
  # create WS for plot ------------------------------------------------------
  create_plot_WS(input.seq, input.num_partitions, input.x_axis, 
                 input.WS_path, input.output_file_name_WS, output_dir, input.type_measure = "auc")
  # create plot -------------------------------------------------------------
  load(paste0(output_dir, input.output_file_name_WS))
  create_plot(input.output_file_name_plot, varying_training_matrix, input.xlab, output_dir)
}

# Bortezomib cp2p var c for p50 ---------------------------------------------------
rm(list=ls(pattern="input.*"))
input.seq <- seq(2, 15, 1)
input.num_partitions <- rep(list(1:100), length(input.seq))
input.x_axis <- seq(60, 300, 20)
stopifnot(length(input.seq) == length(input.x_axis))
input.xlab <- "# of Cell Line Samples in Training Data"
output_dir <- "Bortezomib/plots/"

input.response_labels <- c("cp2p_ic50_p50","cp2p_auc_p50","cp2p_slope_p50",
                           "c2p_ic50_p50","c2p_auc_p50","c2p_slope_p50")
for (i in 1:length(input.response_labels)) {
  input.response_label <- input.response_labels[i]
  input.WS_path <- paste0("Bortezomib/computed_WS/Var C p50/bortezomib_", input.response_label, "_1to100_parInd")
  input.output_file_name_plot <- paste0("bor_", input.response_label, ".png")
  input.output_file_name_WS <- paste0("bor_", input.response_label, ".RData")
  # create WS for plot ------------------------------------------------------
  create_plot_WS(input.seq, input.num_partitions, input.x_axis, 
                 input.WS_path, input.output_file_name_WS, output_dir)
  # create plot -------------------------------------------------------------
  load(paste0(output_dir, input.output_file_name_WS))
  create_plot(input.output_file_name_plot, varying_training_matrix, input.xlab, output_dir)
}

# Bortezomib cp2p var c for p100 ---------------------------------------------------
rm(list=ls(pattern="input.*"))
input.seq <- seq(3, 15, 1)
input.num_partitions <- rep(list(1:100), length(input.seq))
input.x_axis <- seq(60, 300, 20)
stopifnot(length(input.seq) == length(input.x_axis))
input.xlab <- "# of Cell Line Samples in Training Data"
output_dir <- "Bortezomib/plots/"

input.response_labels <- c("cp2p_ic50_p100","cp2p_auc_p100","cp2p_slope_p100",
                           "c2p_ic50_p100","c2p_auc_p100","c2p_slope_p100")
for (i in 1:length(input.response_labels)) {
  input.response_label <- input.response_labels[i]
  # based on response_label
  input.WS_path <- paste0("Bortezomib/computed_WS/Var C p100/bortezomib_", input.response_label, "_1to100_parInd")
  input.output_file_name_plot <- paste0("bor_", input.response_label, ".png")
  input.output_file_name_WS <- paste0("bor_", input.response_label, ".RData")
  # create WS for plot ------------------------------------------------------
  create_plot_WS(input.seq, input.num_partitions, input.x_axis, 
                 input.WS_path, input.output_file_name_WS, output_dir, input.type_measure = "auc")
  # create plot -------------------------------------------------------------
  load(paste0(output_dir, input.output_file_name_WS))
  create_plot(input.output_file_name_plot, varying_training_matrix, input.xlab, output_dir)
}

# Bortezomib best var p plot ----------------------------------------------
load(paste0(output_dir, "bor_p2p.RData"))
p2p.results <- varying_training_matrix
p2p.x_axis <- input.x_axis 
rm("varying_training_matrix")
load(paste0(output_dir, "bor_cp2p_slope.RData"))
cp2p.results <- varying_training_matrix
rm("varying_training_matrix")
load(paste0(output_dir, "bor_c2p_slope_p50.RData"))
c2p.results <- varying_training_matrix
rm("varying_training_matrix")
input.output_file_name_plot = "bortezomib_var_p_best.png"
create_plot_best_var_p(p2p.x_axis, p2p.results, cp2p.results, 
                       c2p.results, input.x_axis_label,
                       input.output_file_name_plot,
                       output_dir)

# Bortezomib best var c p50 plot ----------------------------------------------
rm(cp2p.results)
load(paste0(output_dir, "bor_p2p.RData"))
p2p.results <- varying_training_matrix
rm("varying_training_matrix")
load(paste0(output_dir, "bor_cp2p_slope_p50.RData"))
cp2p.results <- varying_training_matrix
rm("varying_training_matrix")
load(paste0(output_dir, "bor_c2p_slope_p50.RData"))
c2p.results <- varying_training_matrix
c2p.x_axis <- input.x_axis 
rm("varying_training_matrix")
input.output_file_name_plot = "bortezomib_var_c_p50_best.png"

create_plot_best_var_c(c2p.x_axis, p2p.results, cp2p.results, 
                       c2p.results,
                       input.output_file_name_plot,
                       output_dir)

# Bortezomib best var c p100 plot ----------------------------------------------
rm(cp2p.results)
load(paste0(output_dir, "bor_p2p.RData"))
p2p.results <- varying_training_matrix
rm("varying_training_matrix")
load(paste0(output_dir, "bor_cp2p_slope_p100.RData"))
cp2p.results <- varying_training_matrix
rm("varying_training_matrix")
load(paste0(output_dir, "bor_c2p_slope_p100.RData"))
c2p.results <- varying_training_matrix
c2p.x_axis <- input.x_axis 
rm("varying_training_matrix")
input.output_file_name_plot = "bortezomib_var_c_p100_best.png"

create_plot_best_var_c(c2p.x_axis, p2p.results, cp2p.results, 
                       c2p.results,
                       input.output_file_name_plot,
                       output_dir)

