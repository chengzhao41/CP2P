library("doParallel")
library("RColorBrewer")
library("ROCR")

cpp.varying_training_matrix <- matrix(nrow = 0, ncol = 21)

label.type <- "auc_p100"

for (temp.num in 1:7) {
  #load(paste0("~/output_temp/Epirubicin/cpp_var/", label.type, "/epirubicin_cpp_1to100_", label.type, "_var", temp.num, ".RData"))   
  
  load(paste0("/Users/chengzhao/Dropbox/output_WS/epirubicin_cpp_1to100_", label.type, "_var", temp.num, ".RData")) 
  
  # Boxplot Data ---------------------------------
  cpp.boxplot_matrix <- foreach (temp.ind = 1:100, .combine=rbind) %do% {
    c(
      cpp.other_model.l1000[[temp.ind]]$elasticNet$test_auc
      , cpp.other_model.l1000[[temp.ind]]$lasso$test_auc
      , cpp.other_model.l1000[[temp.ind]]$ridge$test_auc
      , cpp.other_model.l1000[[temp.ind]]$rf$test_auc
      , cpp.other_model.l1000[[temp.ind]]$svmLinear$test_auc
      , cpp.other_model.l1000[[temp.ind]]$svmRadial$test_auc
      , cpp.snf.single.l1000[[temp.ind]]$test.auc    
      
      , cpp.other_model.all[[temp.ind]]$elasticNet$test_auc
      , cpp.other_model.all[[temp.ind]]$lasso$test_auc
      , cpp.other_model.all[[temp.ind]]$ridge$test_auc
      , cpp.other_model.all[[temp.ind]]$rf$test_auc
      , cpp.other_model.all[[temp.ind]]$svmLinear$test_auc
      , cpp.other_model.all[[temp.ind]]$svmRadial$test_auc
      , cpp.snf.single.all[[temp.ind]]$test.auc
      
      , cpp.other_model.mRMR1000[[temp.ind]]$elasticNet$test_auc
      , cpp.other_model.mRMR1000[[temp.ind]]$lasso$test_auc
      , cpp.other_model.mRMR1000[[temp.ind]]$ridge$test_auc
      , cpp.other_model.mRMR1000[[temp.ind]]$rf$test_auc
      , cpp.other_model.mRMR1000[[temp.ind]]$svmLinear$test_auc
      , cpp.other_model.mRMR1000[[temp.ind]]$svmRadial$test_auc
      , cpp.snf.single.mRMR1000[[temp.ind]]$test.auc    
    )  
  }
  temp.boxplot_data <- data.frame(
    cpp.boxplot_matrix)
  
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
  temp.order <- order(apply(temp.boxplot_data, 2, median))
  #with(temp.boxplot_data, reorder(TYPE, temp.boxplot_data, median))
  
  #dev.off()
  boxplot(temp.boxplot_data, las=2, 
          par(mar = c(10, 5, 4, 2) + 0.1), 
          ylab = "AUC on test set",
          names = temp.names,
          at = rank(apply(temp.boxplot_data, 2, median), ties.method = c("first")))
  axis(las=1, side = 4)
  temp.order <- order(apply(temp.boxplot_data, 2, median))
  temp.median_results <- apply(temp.boxplot_data, 2, median)
  names(temp.median_results) <- temp.names
  temp.median_results[temp.order]
  
  temp.median_results
  cpp.varying_training_matrix <- rbind(cpp.varying_training_matrix, temp.median_results)
}

num_of_cell_lines = seq(5, 35, 5)
stopifnot(dim(cpp.varying_training_matrix)[1] == length(num_of_cell_lines))
save(cpp.varying_training_matrix, num_of_cell_lines, label.type, file = paste0("Epirubicin/output_WS/cpp_var_", label.type, ".RData"))

###
stopifnot(dim(cpp.varying_training_matrix)[2] == 21)

range(cpp.varying_training_matrix)

png(filename = paste0("Epirubicin/plots/epi_cpp_var_", label.type, ".png"), width = 800, height = 800)

if (label.type == "auc_p50") {
  ylim = c(0.4, 0.75)
  yaxp = c(0.4, 0.65, 5)
} else if (label.type == "auc_p100") {
  ylim = c(0.35, 0.85)
  yaxp = c(0.35, 0.75, 8)
} else if (label.type == "slope_p50") {
  ylim = c(0.4, 0.75)
  yaxp = c(0.4, 0.65, 5)
} else if (label.type == "slope_p100") {
  ylim = c(0.35, 0.85)
  yaxp = c(0.35, 0.75, 8)
} else {
  stop("label not defined!")
}

par(mar=c(5.1,7,4.5,3.5))
plot(0, 0, xlim = range(num_of_cell_lines), ylim = ylim, type = "n", 
     ylab = "",
     xlab = "# of Cell Line Samples in Training Data",
     #las = 1, cex.axis = 1.2, cex.lab = 1.5)
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = yaxp)

title(ylab = "Test AUROC", cex.lab = 1.5, line = 4.5)

#axis(side = 4, las = 1, cex.axis = 1.2)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = yaxp)

cl <- rainbow(7)
temp.points <- rep(1, 7)
temp.points <- c(temp.points, rep(2, 7))
temp.points <- c(temp.points, rep(3, 7))
for (temp.ind in 1:21) {
  lines(num_of_cell_lines, cpp.varying_training_matrix[, temp.ind], col = cl[(temp.ind - 1) %% 7 + 1], type = 'b', pch=temp.points[temp.ind])
}
legend("topleft", legend = colnames(cpp.varying_training_matrix), col=cl, pch=temp.points, ncol=3, pt.cex = 1, cex = 1.3) # optional legend

dev.off()

