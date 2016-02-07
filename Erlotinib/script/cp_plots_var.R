library("doParallel")
library("ROCR")
cp.varying_training_matrix <- matrix(nrow = 0, ncol = 21)
label.type <- "slope"

for (temp.num in 1:10) {
  load(paste0("~/output_temp/Erlotinib/cp_var/erlotinib_cp_1to100_", label.type, "_var", temp.num, ".RData"))  
  
  cp.boxplot_matrix <- foreach (temp.ind = 1:100, .combine=rbind) %do% {
    c(
      cp.other_model.l1000[[temp.ind]]$elasticNet$test_auc
      , cp.other_model.l1000[[temp.ind]]$lasso$test_auc
      , cp.other_model.l1000[[temp.ind]]$ridge$test_auc
      , cp.other_model.l1000[[temp.ind]]$rf$test_auc
      , cp.other_model.l1000[[temp.ind]]$svmLinear$test_auc
      , cp.other_model.l1000[[temp.ind]]$svmRadial$test_auc
      , cp.snf.single.l1000[[temp.ind]]$test.auc    
      
      , cp.other_model.all[[temp.ind]]$elasticNet$test_auc
      , cp.other_model.all[[temp.ind]]$lasso$test_auc
      , cp.other_model.all[[temp.ind]]$ridge$test_auc
      , cp.other_model.all[[temp.ind]]$rf$test_auc
      , cp.other_model.all[[temp.ind]]$svmLinear$test_auc
      , cp.other_model.all[[temp.ind]]$svmRadial$test_auc
      , cp.snf.single.all[[temp.ind]]$test.auc
      
      , cp.other_model.mRMR1000[[temp.ind]]$elasticNet$test_auc
      , cp.other_model.mRMR1000[[temp.ind]]$lasso$test_auc
      , cp.other_model.mRMR1000[[temp.ind]]$ridge$test_auc
      , cp.other_model.mRMR1000[[temp.ind]]$rf$test_auc
      , cp.other_model.mRMR1000[[temp.ind]]$svmLinear$test_auc
      , cp.other_model.mRMR1000[[temp.ind]]$svmRadial$test_auc
      , cp.snf.single.mRMR1000[[temp.ind]]$test.auc    
    )  
  }
  temp.boxplot_data <- data.frame(
    cp.boxplot_matrix)
  
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
  cp.varying_training_matrix <- rbind(cp.varying_training_matrix, temp.median_results)
}

setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS")
num_of_cell_lines = seq(30, 300, 30)
stopifnot(dim(cp.varying_training_matrix)[1] == length(num_of_cell_lines))
save(cp.varying_training_matrix, num_of_cell_lines, label.type, file = paste0("cp_var_", label.type, ".RData"))

###
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/cp_var_ic50.RData")
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/cp_var_auc.RData")
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/cp_var_slope.RData")

stopifnot(dim(cp.varying_training_matrix)[2] == 21)
setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/Results/")

min(cp.varying_training_matrix)
max(cp.varying_training_matrix)

png(filename = paste0("erl_cp_var_", label.type, ".png"), width = 800, height = 800)

if (label.type == "ic50") {
  ylim = c(0.4, 1.1)
} else if (label.type == "auc") {
  ylim = c(0.3, 0.9)
} else if (label.type == "slope") {
  ylim = c(0.25, 0.95)
  #par(yaxp = c(0.4, 0.95, 10))
} else {
  stop("label not defined!")
}

par(mar=c(5.1,7,4.5,3.5))
plot(0, 0, xlim = range(num_of_cell_lines), ylim = ylim, type = "n", 
     ylab = "",
     xlab = "# of Cell line Samples in Training Data",
     #las = 1, cex.axis = 1.2, cex.lab = 1.5)
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = c(0.25, 0.8, 11))
title(ylab = "Test AUROC", cex.lab = 1.5, line = 4.5)

#axis(side = 4, las = 1, cex.axis = 1.2)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = c(0.25, 0.8, 11))

cl <- rainbow(7)
temp.points <- rep(1, 7)
temp.points <- c(temp.points, rep(2, 7))
temp.points <- c(temp.points, rep(3, 7))
for (temp.ind in 1:21) {
  lines(num_of_cell_lines, cp.varying_training_matrix[, temp.ind], col = cl[(temp.ind - 1) %% 7 + 1], type = 'b', pch=temp.points[temp.ind])
}
legend("topleft", legend = colnames(cp.varying_training_matrix), col=cl, pch=temp.points, ncol=3, pt.cex = 1, cex = 1.3) # optional legend

dev.off()

