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
  input.ylim = c((floor(range(varying_training_matrix$mean)[1]*10)/10) - 0.05, 
                 (ceiling(range(varying_training_matrix$mean)[2]*10)/10) + 0.1)
  input.yaxp = c(input.ylim[1], min(1, input.ylim[2]), round((min(1, input.ylim[2]) - input.ylim[1]) / 0.05))
  if (input.yaxp[3] >= 15) {
    input.ylim[2] = input.ylim[2] + 0.2
  }
  
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
      input.num_partition = input.num_partitions[[1]]
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

create_plot_best <- function(input.x_axis, 
                             input_x_labels,
                             input.results, 
                             input.results_constant,
                             input.output_file_name_plot, 
                             output_dir, 
                             input.labels,
                             input.labels_constant,
                             include_std = TRUE,
                             input.color_manual,
                             input.shape_manual,
                             input.constant_best_ind) 
{
  stopifnot(!is.null(input.x_axis))
  stopifnot(!is.null(input_x_labels))
  stopifnot(!is.null(input.results))
  stopifnot(!is.null(input.results_constant))
  stopifnot(!is.null(input.output_file_name_plot))
  stopifnot(!is.null(output_dir))
  stopifnot(!is.null(input.labels))
  stopifnot(!is.null(input.labels_constant))
  stopifnot(!is.null(include_std))
  stopifnot(!is.null(input.color_manual))
  stopifnot(!is.null(input.shape_manual))
  stopifnot(length(input.color_manual) == length(input.shape_manual))
  stopifnot(!is.null(input.constant_best_ind))
  stopifnot(length(input.labels) == length(input.results))
  stopifnot(length(input.labels_constant) == length(input.results_constant))
  
  require("ggplot2")
  require("grid")
  best_model <- vector()
  best_model_constant <- vector() 
  
  min_value <- 1
  max_value <- 0
  all.m <- data.frame()
  approach <- data.frame()
  
  for (i in 1:length(input.results)) {
    current_results <- input.results[[i]]
    last_3_ind = (dim(current_results$mean)[1] - 2) : dim(current_results$mean)[1]
    stopifnot(length(last_3_ind) == 3)
    current_best <- which(max(apply(current_results$mean[last_3_ind, ], 2, mean)) == apply(current_results$mean[last_3_ind, ], 2, mean))
    if (length(current_best) > 1) {
      warnings("Many best")
      print(current_best)
      current_best <- current_best[1]
    }
    
    if (include_std) {
      format_results <- cbind(input.x_axis, current_results$mean[, current_best], current_results$std[, current_best])
      colnames(format_results) <- c("xval", "yval", "se")
      min_value <- min(min_value, current_results$mean[, current_best] - current_results$std[, current_best], min_value)   
      max_value <- max(max_value, current_results$mean[, current_best] + current_results$std[, current_best], min_value)
    } else {
      format_results <- cbind(input.x_axis, current_results$mean[, current_best])
      colnames(format_results) <- c("xval", "yval")
      min_value <- min(min_value, current_results$mean[, current_best], min_value)   
      max_value <- max(max_value, current_results$mean[, current_best], min_value)
    }
    
    best_model <- c(best_model, names(current_best))
    all.m <- rbind(all.m, format_results)
    approach <- c(approach, rep(input.labels[i], dim(format_results)[1]))
  }
  
  all.m$Approach <- approach
  all.m$Approach <- factor(all.m$Approach, levels = c(input.labels, input.labels_constant))

  for (i in 1:length(input.labels_constant)) {
    all.m[nrow(all.m) + 1, ] <- all.m[nrow(all.m) + 1, ] # Creates a new row filled with NAs
    all.m$Approach[nrow(all.m)] <- input.labels_constant[i]
    
    current_results <- input.results_constant[[i]]
    current_best <- which(current_results$mean[input.constant_best_ind, ] == max(current_results$mean[input.constant_best_ind, ]))
    if (length(current_best) > 1) {
      warnings("Many best")
      print(current_best)
      current_best <- current_best[1]
    }
    min_value <- min(min_value, current_results$mean[, current_best])   
    max_value <- max(max_value, current_results$mean[, current_best])
    best_model_constant <- c(best_model_constant, current_best)
  }

  min_value <- floor(min_value*20)/20
  max_value <- ceiling(max_value*20)/20
  calculated_size = length(input.results) + length(input.results_constant)
  
  # The errorbars overlapped, so use position_dodge to move them horizontally
  pd <- position_dodge(2) # move them .05 to the left and right
  
  # get the colours
  plot_best <- ggplot(all.m, aes(x=xval, y=yval, colour = Approach, ymax = max_value, ymin = min_value)) + 
    theme_bw() + 
    geom_line(position=pd) + 
    geom_point(aes(shape=Approach, colour = Approach), size = calculated_size, na.rm = TRUE)
  
  if (include_std) {
    plot_best <- plot_best + geom_errorbar(aes(ymin= yval - se, ymax = yval + se), width=5, position=pd)
  }
  for (i in 1:length(input.results_constant)) {
    temp.yintercept <- input.results_constant[[i]]$mean[input.constant_best_ind, best_model_constant[i]]
    plot_best <- plot_best + geom_hline(aes(yintercept = temp.yintercept, colour = input.labels_constant[i]))
  }
  
  plot_best <- plot_best + 
    scale_color_manual(values = input.color_manual) + 
    scale_shape_manual(values = input.shape_manual) +
    scale_y_continuous(breaks = seq(min_value, max_value, 0.05), "Test AUROC") +
    scale_x_continuous(breaks = input.x_axis, input_x_labels) + 
    theme(legend.direction = 'horizontal', 
    legend.position = 'top', 
    plot.margin = unit(c(5.1, 7, 4.5, 3.5)/2, "lines"), 
    text = element_text(size=15), axis.title.x=element_text(vjust=0), axis.title.y=element_text(vjust=2))
  
  png(filename = paste0(output_dir, input.output_file_name_plot), width = 800, height = 800)
  print(plot_best)
  dev.off()
}
