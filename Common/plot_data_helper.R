create_plot <- function(input.output_file_name_plot, 
                        varying_training_matrix, 
                        input.xlab, 
                        output_dir
                        ) {
  stopifnot(!is.null(input.output_file_name_plot))
  stopifnot(!is.null(varying_training_matrix))
  stopifnot(!is.null(input.x_axis))
  stopifnot(!is.null(output_dir))
  
  # create the png file
  png(filename = paste0(output_dir, input.output_file_name_plot), width = 800, height = 800)
  input.ylim = c((round(range(varying_training_matrix$mean)*10)/10)[1] - 0.05, 
                 (round(range(varying_training_matrix$mean)*10)/10)[2] + 0.15)
  input.yaxp = c(input.ylim, round((input.ylim[2] - input.ylim[1]) / 0.05))
  
  par(mar=c(5.1,7,4.5,3.5))
  plot(0, 0, xlim = range(input.x_axis), ylim = input.ylim, type = "n", 
       ylab = "",
       xlab = input.xlab,
       las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = input.yaxp)
  title(ylab = "Test AUROC", cex.lab = 1.5, line =4.5)
  axis(side = 4, las = 1, cex.axis = 1.2, yaxp = input.yaxp)
  
  cl <- rainbow(7)
  temp.points <- rep(1, 7)
  temp.points <- c(temp.points, rep(2, 7))
  temp.points <- c(temp.points, rep(3, 7))
  for (temp.ind in 1:21) {
    lines(input.x_axis, varying_training_matrix$mean[, temp.ind], col = cl[(temp.ind - 1) %% 7 + 1], type = 'b', pch=temp.points[temp.ind])  
  }
  legend("topleft", legend = colnames(varying_training_matrix$mean), col=cl, pch=temp.points, ncol=3, pt.cex = 1, cex = 1.3)
  dev.off()
}

create_plot_WS <- function(input.seq, input.num_partitions, input.x_axis, 
                           input.WS_paths, input.output_file_name_WS, output_dir, input.type_measure)
{
  stopifnot(!is.null(input.seq))
  stopifnot(!is.null(input.num_partitions))
  stopifnot(!is.null(input.x_axis))
  stopifnot(!is.null(input.WS_paths))
  stopifnot(!is.null(input.output_file_name_WS))
  stopifnot(!is.null(output_dir))
  stopifnot(input.type_measure %in% c("acc", "auc"))
  
  # processing each WS ------------------------------------------------------
  varying_training_matrix <- list()
  varying_training_matrix$mean <- matrix(nrow = 0, ncol = 21)
  varying_training_matrix$std <- matrix(nrow = 0, ncol = 21)
  
  temp.loop_ind = 0
  for (temp.num in input.seq) {
    temp.loop_ind = temp.loop_ind + 1
    if (length(input.WS_paths) == 1) {
      load(paste0(input.WS_paths, temp.num, ".RData"))
    } else {
      stopifnot(length(input.WS_paths) == length(input.seq))
      load(paste0(input.WS_paths[temp.loop_ind], temp.num, ".RData"))
    }
    if (length(input.num_partitions) == 1) {
      input.num_partition = input.num_partitions
    } else {
      stopifnot(length(input.num_partitions) == length(input.seq))
      input.num_partition = input.num_partitions[[temp.loop_ind]]
    }
    
    # Boxplot Data ---------------------------------
    if (input.type_measure == "auc") {
      boxplot_matrix <- foreach (temp.ind = input.num_partition, .combine=rbind) %do% {
        c(
          other_model.l1000[[temp.ind]]$elasticNet$test_auc
          , other_model.l1000[[temp.ind]]$lasso$test_auc
          , other_model.l1000[[temp.ind]]$ridge$test_auc
          , other_model.l1000[[temp.ind]]$rf$test_auc
          , other_model.l1000[[temp.ind]]$svmLinear$test_auc
          , other_model.l1000[[temp.ind]]$svmRadial$test_auc
          , snf.single.l1000[[temp.ind]]$test.auc    
          
          , other_model.all[[temp.ind]]$elasticNet$test_auc
          , other_model.all[[temp.ind]]$lasso$test_auc
          , other_model.all[[temp.ind]]$ridge$test_auc
          , other_model.all[[temp.ind]]$rf$test_auc
          , other_model.all[[temp.ind]]$svmLinear$test_auc
          , other_model.all[[temp.ind]]$svmRadial$test_auc
          , snf.single.all[[temp.ind]]$test.auc
          
          , other_model.mRMR1000[[temp.ind]]$elasticNet$test_auc
          , other_model.mRMR1000[[temp.ind]]$lasso$test_auc
          , other_model.mRMR1000[[temp.ind]]$ridge$test_auc
          , other_model.mRMR1000[[temp.ind]]$rf$test_auc
          , other_model.mRMR1000[[temp.ind]]$svmLinear$test_auc
          , other_model.mRMR1000[[temp.ind]]$svmRadial$test_auc
          , snf.single.mRMR1000[[temp.ind]]$test.auc    
        )
      }
    } else if (input.type_measure == "acc") {
      temp.predictions <- matrix(nrow = 0, ncol = 21)
      for (temp.ind in 1:length(input.num_partition)) {
        temp.row <- c(
          other_model.l1000[[temp.ind]]$elasticNet$predict
          , other_model.l1000[[temp.ind]]$lasso$predict
          , other_model.l1000[[temp.ind]]$ridge$predict
          , other_model.l1000[[temp.ind]]$rf$predict
          , other_model.l1000[[temp.ind]]$svmLinear$predict
          , other_model.l1000[[temp.ind]]$svmRadial$predict
          , snf.single.l1000[[temp.ind]]$test.prediction    
          
          , other_model.all[[temp.ind]]$elasticNet$predict
          , other_model.all[[temp.ind]]$lasso$predict
          , other_model.all[[temp.ind]]$ridge$predict
          , other_model.all[[temp.ind]]$rf$predict
          , other_model.all[[temp.ind]]$svmLinear$predict
          , other_model.all[[temp.ind]]$svmRadial$predict
          , snf.single.all[[temp.ind]]$test.prediction
          
          , other_model.mRMR1000[[temp.ind]]$elasticNet$predict
          , other_model.mRMR1000[[temp.ind]]$lasso$predict
          , other_model.mRMR1000[[temp.ind]]$ridge$predict
          , other_model.mRMR1000[[temp.ind]]$rf$predict
          , other_model.mRMR1000[[temp.ind]]$svmLinear$predict
          , other_model.mRMR1000[[temp.ind]]$svmRadial$predict
          , snf.single.mRMR1000[[temp.ind]]$test.prediction    
        )
        
        stopifnot(length(temp.row) %% 21 == 0)       
        temp.predictions = rbind(temp.predictions, matrix(temp.row, length(temp.row) / 21, 21))
      }
    } else {
      stop("invalid input.type_measure")
    }

    if (input.type_measure == "auc") {
      temp.boxplot_data <- data.frame(boxplot_matrix)
      temp.mean_results <- apply(temp.boxplot_data, 2, mean)
      temp.std_results <- apply(temp.boxplot_data, 2, sd)
      varying_training_matrix$mean <- rbind(varying_training_matrix$mean, temp.mean_results)
      varying_training_matrix$std <- rbind(varying_training_matrix$std, temp.std_results)
    } else {
      require("ROCR")
      temp.test_index <- foreach (temp.ind = 1:length(input_partition), .combine=c, .errorhandling = 'remove') %do% {
        input_partition[[temp.ind]]$test_index
      }
      stopifnot(dim(temp.predictions)[1] == length(temp.test_index))
      temp.test_auc <- foreach (temp.model_ind = 1:21, .combine=cbind) %do% {
        temp.predict <- prediction(temp.predictions[, temp.model_ind], input_label[temp.test_index])
        unlist(slot(performance(temp.predict, "auc"), "y.values"))  
      }
      varying_training_matrix$mean <- rbind(varying_training_matrix$mean, temp.test_auc)
    }
  }
  stopifnot(dim(varying_training_matrix$mean)[1] == length(input.seq))
  stopifnot(dim(varying_training_matrix$mean)[2] == 21)
  
  temp.names <- c(
    "Elastic Net L1000" 
    , "Lasso L1000" 
    , "Ridge L1000" 
    , "RF L1000" 
    , "SVM lin L1000" 
    , "SVM rbf L1000"           
    , "SNF L1000"   
    
    , "Elastic Net all" 
    , "Lasso all"   
    , "Ridge all" 
    , "RF all" 
    , "SVM lin all" 
    , "SVM rbf all"     
    , "SNF all"
    
    , "Elastic Net mRMR1000" 
    , "Lasso mRMR1000"   
    , "Ridge mRMR1000" 
    , "RF mRMR1000" 
    , "SVM lin mRMR1000" 
    , "SVM rbf mRMR1000"       
    , "SNF mRMR1000"
  )
  colnames(varying_training_matrix$mean) <- temp.names
  if (!is.null(varying_training_matrix$std)) {
    colnames(varying_training_matrix$std) <- temp.names
  }
  
  # save the WS for the plot
  save(varying_training_matrix, input.x_axis, file = paste0(output_dir, input.output_file_name_WS))
}

create_plot_best_var_p <- function(p2p.x_axis, p2p.results, cp2p.results, 
                                   c2p.results, input.x_axis_label,
                                   input.output_file_name_plot, 
                                   output_dir) 
{
  stopifnot(!is.null(p2p.x_axis))
  stopifnot(!is.null(p2p.results))
  stopifnot(!is.null(cp2p.results))
  stopifnot(!is.null(c2p.results))
  stopifnot(!is.null(input.x_axis_label))
  
  require("ggplot2")
  require("grid")
  require("RColorBrewer")
  best_model <- list()
  
  # p2p
  last_3_ind = (dim(p2p.results$mean)[1] - 2) : dim(p2p.results$mean)[1]
  stopifnot(length(last_3_ind) == 3)
  best_model$p2p_ind <- which(max(apply(p2p.results$mean[last_3_ind, ], 2, mean)) == apply(p2p.results$mean[last_3_ind, ], 2, mean))
  stopifnot(length(best_model$p2p_ind) == 1)
  # cp2p
  last_3_ind = (dim(cp2p.results$mean)[1] - 2) : dim(cp2p.results$mean)[1]
  stopifnot(length(last_3_ind) == 3)
  best_model$cp2p_ind <- which(max(apply(cp2p.results$mean[last_3_ind, ], 2, mean)) == apply(cp2p.results$mean[last_3_ind, ], 2, mean))
  stopifnot(length(best_model$cp2p_ind) == 1)
  # c2p
  last_3_ind = (dim(c2p.results$mean)[1] - 2) : dim(c2p.results$mean)[1]
  stopifnot(length(last_3_ind) == 3)
  best_model$c2p_ind <- which(max(apply(c2p.results$mean[last_3_ind, ], 2, mean)) == apply(c2p.results$mean[last_3_ind, ], 2, mean))
  stopifnot(length(best_model$cp2p_ind) == 1)
  
  temp.cp2p <- data.frame(cbind(x_values = p2p.x_axis, y_values = cp2p.results$mean[, best_model$cp2p_ind], 
                                std_values = cp2p.results$std[, best_model$cp2p_ind]))
  colnames(temp.cp2p) <- c("xval", "yval", "se")
  
  temp.p2p <- cbind(p2p.x_axis, p2p.results$mean[, best_model$p2p_ind], p2p.results$std[, best_model$p2p_ind])
  colnames(temp.p2p) <- c("xval", "yval", "se")
  
  all.m <- rbind(temp.cp2p, temp.p2p)
  all.m$Approach <- c(rep("CP2P", dim(temp.cp2p)[1]), rep("P2P", dim(temp.p2p)[1]))
  all.m$Approach <- factor(all.m$Approach, levels = c("C2P", "P2P", "CP2P"))
  
  all.m[nrow(all.m) + 1, ] <- all.m[nrow(all.m) + 1, ] # Creates a new row filled with NAs
  all.m$Approach[nrow(all.m)] <- "C2P"
  
  # The errorbars overlapped, so use position_dodge to move them horizontally
  pd <- position_dodge(2) # move them .05 to the left and right
  cols <- brewer.pal(n = 3, name = 'Dark2')
  
  plot_best <- ggplot(all.m, aes(x=xval, y=yval, colour = Approach, ymax = 0.95)) + theme_bw() + 
    geom_errorbar(aes(ymin= yval - se, ymax = yval + se), width=5, position=pd) + 
    geom_line(position=pd) + 
    geom_point(aes(shape=Approach, colour = Approach), size = 3, na.rm = TRUE) + 
    geom_hline(aes(yintercept = c2p.results$mean[dim(c2p.results$mean)[1], best_model$c2p_ind], colour = "C2P")) + 
    
    scale_color_manual(values = c("C2P" = cols[1], "P2P" = cols[2], "CP2P" = cols[3])) + 
    scale_shape_manual(values = c("C2P" = NA, "P2P" = 16, "CP2P" = 17)) +
    scale_y_continuous(breaks = seq(0.4, 0.95, 0.05), "Test AUROC") +
    scale_x_continuous(breaks = p2p.x_axis, "# Number of Patient Samples in Training")
  plot_best <- plot_best + theme(legend.direction = 'horizontal', 
                       legend.position = 'top', 
                       plot.margin = unit(c(5.1, 7, 4.5, 3.5)/2, "lines"), 
                       text = element_text(size=15), axis.title.x=element_text(vjust=-1.5), axis.title.y=element_text(vjust=2))   
  png(filename = paste0(output_dir, input.output_file_name_plot), width = 800, height = 800)
  print(plot_best)
  dev.off()
}

create_plot_best_var_c <- function(c2p.x_axis, p2p.results, cp2p.results, 
                                   c2p.results,
                                   input.output_file_name_plot, 
                                   output_dir) 
{
  stopifnot(!is.null(c2p.x_axis))
  stopifnot(!is.null(p2p.results))
  stopifnot(!is.null(cp2p.results))
  stopifnot(!is.null(c2p.results))
  stopifnot(!is.null(input.output_file_name_plot))
  stopifnot(!is.null(output_dir))
  
  require("ggplot2")
  require("grid")
  require("RColorBrewer")
  best_model <- list()
  
  # p2p
  last_3_ind = (dim(p2p.results$mean)[1] - 2) : dim(p2p.results$mean)[1]
  stopifnot(length(last_3_ind) == 3)
  best_model$p2p_ind <- which(max(apply(p2p.results$mean[last_3_ind, ], 2, mean)) == apply(p2p.results$mean[last_3_ind, ], 2, mean))
  stopifnot(length(best_model$p2p_ind) == 1)
  # cp2p
  last_3_ind = (dim(cp2p.results$mean)[1] - 2) : dim(cp2p.results$mean)[1]
  stopifnot(length(last_3_ind) == 3)
  best_model$cp2p_ind <- which(max(apply(cp2p.results$mean[last_3_ind, ], 2, mean)) == apply(cp2p.results$mean[last_3_ind, ], 2, mean))
  stopifnot(length(best_model$cp2p_ind) == 1)
  # c2p
  last_3_ind = (dim(c2p.results$mean)[1] - 2) : dim(c2p.results$mean)[1]
  stopifnot(length(last_3_ind) == 3)
  best_model$c2p_ind <- which(max(apply(c2p.results$mean[last_3_ind, ], 2, mean)) == apply(c2p.results$mean[last_3_ind, ], 2, mean))
  stopifnot(length(best_model$cp2p_ind) == 1)
  
  temp.cp2p <- data.frame(cbind(x_values = c2p.x_axis, y_values = cp2p.results$mean[, best_model$cp2p_ind], 
                                std_values = cp2p.results$std[, best_model$cp2p_ind]))
  colnames(temp.cp2p) <- c("xval", "yval", "se")
  
  temp.c2p <- cbind(c2p.x_axis, c2p.results$mean[, best_model$c2p_ind], c2p.results$std[, best_model$c2p_ind])
  colnames(temp.c2p) <- c("xval", "yval", "se")
  
  all.m <- rbind(temp.cp2p, temp.c2p)
  all.m$Approach <- c(rep("CP2P", dim(temp.cp2p)[1]), rep("C2P", dim(temp.c2p)[1]))
  all.m$Approach <- factor(all.m$Approach, levels = c("P2P", "C2P", "CP2P"))
  
  all.m[nrow(all.m) + 1, ] <- all.m[nrow(all.m) + 1, ] # Creates a new row filled with NAs
  all.m$Approach[nrow(all.m)] <- "P2P"
  
  # The errorbars overlapped, so use position_dodge to move them horizontally
  pd <- position_dodge(2) # move them .05 to the left and right
  cols <- brewer.pal(n = 3, name = 'Dark2')
  temp.min_value <- min(
    min(c2p.results$mean[, best_model$c2p_ind] - c2p.results$std[, best_model$c2p_ind]),
    p2p.results$mean[dim(p2p.results$mean)[1], best_model$p2p_ind],
    min(cp2p.results$mean[, best_model$cp2p_ind] - cp2p.results$std[, best_model$cp2p_ind]))
  temp.min_value <- round(temp.min_value*10)/10-0.05
  temp.max_value <- max(
    max(c2p.results$mean[, best_model$c2p_ind] + c2p.results$std[, best_model$c2p_ind]),
    p2p.results$mean[dim(p2p.results$mean)[1], best_model$p2p_ind],
    max(cp2p.results$mean[, best_model$cp2p_ind] + cp2p.results$std[, best_model$cp2p_ind]))
  temp.max_value <- round(temp.max_value*10)/10
  
  plot_best <- ggplot(all.m, aes(x=xval, y=yval, colour = Approach, ymax = temp.max_value)) + theme_bw() + 
    geom_errorbar(aes(ymin= yval - se, ymax = yval + se), width=5, position=pd) + 
    geom_line(position=pd) + 
    geom_point(aes(shape=Approach, colour = Approach), size = 3, na.rm = TRUE) + 
    geom_hline(aes(yintercept = p2p.results$mean[dim(p2p.results$mean)[1], best_model$p2p_ind], colour = "P2P")) + 
    
    scale_color_manual(values = c("P2P" = cols[1], "C2P" = cols[2], "CP2P" = cols[3])) + 
    scale_shape_manual(values = c("P2P" = NA, "C2P" = 16, "CP2P" = 17)) +
    scale_y_continuous(breaks = seq(temp.min_value, temp.max_value, 0.05), "Test AUROC") +
    scale_x_continuous(breaks = c2p.x_axis, "# Number of Cell Line Samples in Training")
  plot_best <- plot_best + theme(legend.direction = 'horizontal', 
                                 legend.position = 'top', 
                                 plot.margin = unit(c(5.1, 7, 4.5, 3.5)/2, "lines"), 
                                 text = element_text(size=15), axis.title.x=element_text(vjust=-1.5), axis.title.y=element_text(vjust=2))   
  png(filename = paste0(output_dir, input.output_file_name_plot), width = 800, height = 800)
  print(plot_best)
  dev.off()
}
