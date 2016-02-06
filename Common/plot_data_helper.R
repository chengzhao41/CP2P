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
                           input.WS_path, input.output_file_name_WS, output_dir)
{
  stopifnot(!is.null(input.seq))
  stopifnot(!is.null(input.num_partitions))
  stopifnot(!is.null(input.x_axis))
  stopifnot(!is.null(input.WS_path))
  stopifnot(!is.null(input.output_file_name_WS))
  stopifnot(!is.null(output_dir))
  
  # processing each WS ------------------------------------------------------
  varying_training_matrix <- list()
  varying_training_matrix$mean <- matrix(nrow = 0, ncol = 21)
  varying_training_matrix$std <- matrix(nrow = 0, ncol = 21)
  
  for (temp.num in input.seq) {
    load(paste0(input.WS_path, temp.num, ".RData"))
    
    # Boxplot Data ---------------------------------
    boxplot_matrix <- foreach (temp.ind = input.num_partitions, .combine=rbind) %do% {
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
    temp.boxplot_data <- data.frame(boxplot_matrix)
    
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
    temp.order <- order(apply(temp.boxplot_data, 2, mean))
    
    boxplot(temp.boxplot_data, las=2, 
            par(mar = c(10, 5, 4, 2) + 0.1), 
            ylab = "AUROC on test set",
            names = temp.names,
            at = rank(apply(temp.boxplot_data, 2, mean), ties.method = c("first")))
    axis(las=1, side = 4)
    temp.order <- order(apply(temp.boxplot_data, 2, mean))
    temp.mean_results <- apply(temp.boxplot_data, 2, mean)
    temp.std_results <- apply(temp.boxplot_data, 2, sd)
    names(temp.mean_results) <- temp.names
    temp.mean_results[temp.order]
    
    temp.mean_results
    varying_training_matrix$mean <- rbind(varying_training_matrix$mean, temp.mean_results)
    varying_training_matrix$std <- rbind(varying_training_matrix$std, temp.std_results)
  }
  stopifnot(dim(varying_training_matrix$mean)[1] == length(input.seq))
  stopifnot(dim(varying_training_matrix$mean)[2] == 21)
  # save the WS for the plot
  save(varying_training_matrix, input.x_axis, file = paste0(output_dir, input.output_file_name_WS))
}
