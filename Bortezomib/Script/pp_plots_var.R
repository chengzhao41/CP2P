stop("Set working directory to current source file")
#setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Code/Bortezomib/Script")

library("doParallel")
library("RColorBrewer")

pp.varying_training_matrix <- matrix(nrow = 0, ncol = 21)

for (temp.num in 2:15) {
  load(paste0("~/output_temp/Bortezomib/pp_var/bortezomib_pp_1to100_", temp.num, "_var", temp.num, ".RData"))   
  load(paste0("~/output_temp/Bortezomib/pp_var/bortezomib_pp_1to100_", temp.num, "_var", temp.num, "_l1000.RData"))   
  # Boxplot Data ---------------------------------
  pp.boxplot_matrix <- foreach (temp.ind = 1:100, .combine=rbind) %do% {
    c(
      pp.other_model.l1000[[temp.ind]]$elasticNet$test_auc
      , pp.other_model.l1000[[temp.ind]]$lasso$test_auc
      , pp.other_model.l1000[[temp.ind]]$ridge$test_auc
      , pp.other_model.l1000[[temp.ind]]$rf$test_auc
      , pp.other_model.l1000[[temp.ind]]$svmLinear$test_auc
      , pp.other_model.l1000[[temp.ind]]$svmRadial$test_auc
      , pp.snf.single.l1000[[temp.ind]]$test.auc    
      
      , pp.other_model.all[[temp.ind]]$elasticNet$test_auc
      , pp.other_model.all[[temp.ind]]$lasso$test_auc
      , pp.other_model.all[[temp.ind]]$ridge$test_auc
      , pp.other_model.all[[temp.ind]]$rf$test_auc
      , pp.other_model.all[[temp.ind]]$svmLinear$test_auc
      , pp.other_model.all[[temp.ind]]$svmRadial$test_auc
      , pp.snf.single.all[[temp.ind]]$test.auc
      
      , pp.other_model.mRMR1000[[temp.ind]]$elasticNet$test_auc
      , pp.other_model.mRMR1000[[temp.ind]]$lasso$test_auc
      , pp.other_model.mRMR1000[[temp.ind]]$ridge$test_auc
      , pp.other_model.mRMR1000[[temp.ind]]$rf$test_auc
      , pp.other_model.mRMR1000[[temp.ind]]$svmLinear$test_auc
      , pp.other_model.mRMR1000[[temp.ind]]$svmRadial$test_auc
      , pp.snf.single.mRMR1000[[temp.ind]]$test.auc    
    )  
  }
  temp.boxplot_data <- data.frame(
    pp.boxplot_matrix)
  
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
  pp.varying_training_matrix <- rbind(pp.varying_training_matrix, temp.median_results)
}

setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Bortezomib/WS")
num_of_patients = seq(20, 150, 10)
save(pp.varying_training_matrix, num_of_patients, file = "pp_var.RData")


#load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Bortezomib/WS/pp_var.RData")
stopifnot(dim(pp.varying_training_matrix)[2] == 21)
dev.off()

png(filename = paste0("bor_pp_var.png"), width = 800, height = 800)

setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Bortezomib/Results/")

par(mar=c(5.1,7,4.5,3.5))
plot(0, 0, xlim = range(num_of_patients), ylim = c(0.5, 0.9), type = "n", 
     ylab = "",
     xlab = "# of Patient Samples in Training Data",
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = c(0.5, 0.8, 6))
     
title(ylab = "Test AUROC", cex.lab = 1.5, line =4.5)

axis(side = 4, las = 1, cex.axis = 1.2, yaxp = c(0.55, 0.8, 5))

cl <- rainbow(7)
temp.points <- rep(1, 7)
temp.points <- c(temp.points, rep(2, 7))
temp.points <- c(temp.points, rep(3, 7))
for (temp.ind in 1:21) {
  if (temp.ind == 17) {
    lines(num_of_patients[c(1, 4:14)], pp.varying_training_matrix[c(1, 4:14), temp.ind], col = cl[(temp.ind - 1) %% 7 + 1], type = 'b', pch=temp.points[temp.ind])  
  } else {
    lines(num_of_patients, pp.varying_training_matrix[, temp.ind], col = cl[(temp.ind - 1) %% 7 + 1], type = 'b', pch=temp.points[temp.ind])  
  }  
}
legend("topleft", legend = colnames(pp.varying_training_matrix), col=cl, pch=temp.points, ncol=3, pt.cex = 1, cex = 1.3)

dev.off()
