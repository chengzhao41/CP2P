# load library ------------------------------------------------------------
library("doParallel")
library("RColorBrewer")

source('Common/plot_data_helper.R')

# input parameters --------------------------------------------------------

# for Bortezomib p2p
input.seq <- seq(3, 15, 1)
input.num_partitions <- 1:100
input.x_axis <- seq(30, 150, 10)
input.WS_path <- "Bortezomib/computed_WS/bortezomib_p2p_1to100_parInd"
input.output_file_name_plot <- "bor_p2p.png"
input.output_file_name_WS <- "bor_p2p.RData"
input.xlab <- "# of Patient Samples in Training Data"
output_dir <- "Bortezomib/plots/"

# for Bortezomib cp2p

input.x_axis <- seq(10, 150, 10)
input.WS_path <- "bortezomib_cp2p_ic50_1to100_parInd"
input.output_file_name_plot <- "bor_c2p_auc_var.png"
input.output_file_name_WS <- "bor_c2p_auc_var.RData"

# create WS for plot
create_plot_WS(input.seq, input.num_partitions, input.x_axis, 
               input.WS_path, input.output_file_name_WS, output_dir)

# create the plot
load(paste0(output_dir, input.output_file_name_WS))
create_plot(input.output_file_name_plot, varying_training_matrix, input.xlab, output_dir)
