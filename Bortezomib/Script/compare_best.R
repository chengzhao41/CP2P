library("doParallel")

#### get the best patient data ########
pp.best <- list()
pp.best$matrix <- matrix(nrow = 100, ncol = 0)
for (temp.num in 2:15) {
  load(paste0("~/output_temp/Bortezomib/pp_var/bortezomib_pp_1to100_", temp.num, "_var", temp.num, ".RData"))   
  pp.current_temp <- foreach (temp.ind = 1:100, .combine=rbind) %do% {
    c(pp.other_model.all[[temp.ind]]$svmLinear$test_auc)  
  }     
  pp.best$matrix <- cbind(pp.best$matrix, pp.current_temp)
}

pp.best$medium <- apply(pp.best$matrix, 2, quantile, 0.5)
pp.best$low <- apply(pp.best$matrix, 2, quantile, 0.25)
pp.best$high <- apply(pp.best$matrix, 2, quantile, 0.75)
pp.best$mean <- apply(pp.best$matrix, 2, mean)
pp.best$sdev <- apply(pp.best$matrix, 2, sd)

pp.best$x_values <- seq(20, 150, 10)

#############################################
label.type <- "slope"
cpp.best <- list()
cpp.best$slope <- matrix(nrow = 100, ncol = 0)

for (temp.num in 1:15) {
  load(paste0("~/output_temp/Bortezomib/cpp_var/slope/bortezomib_cpp_1to100_slope_var", temp.num, ".RData"))   
  cpp.temp_slope <- foreach (temp.ind = 1:100, .combine=rbind) %do% {
    cpp.other_model.all[[temp.ind]]$svmLinear$test_auc
  } 
  cpp.best$slope <- cbind(cpp.best$slope, cpp.temp_slope)
}

cpp.best$slope.medium <- apply(cpp.best$slope, 2, quantile, 0.5)
cpp.best$slope.low <- apply(cpp.best$slope, 2, quantile, 0.25)
cpp.best$slope.high <- apply(cpp.best$slope, 2, quantile, 0.75)
cpp.best$slope.mean <- apply(cpp.best$slope, 2, mean)
cpp.best$slope.sdev <- apply(cpp.best$slope, 2, sd)
cpp.best$x_values <- seq(10, 150, 10)

#############################################
label.type <- "auc"

cpp.best$auc <- matrix(nrow = 100, ncol = 0)

for (temp.num in 1:15) {
  load(paste0("~/output_temp/Bortezomib/cpp_var/auc/bortezomib_cpp_1to100_auc_var", temp.num, ".RData"))   
  cpp.temp_auc <- foreach (temp.ind = 1:100, .combine=rbind) %do% {
    cpp.other_model.all[[temp.ind]]$svmLinear$test_auc
  } 
  cpp.best$auc <- cbind(cpp.best$auc, cpp.temp_auc)
}

cpp.best$auc.medium <- apply(cpp.best$auc, 2, quantile, 0.5)
cpp.best$auc.low <- apply(cpp.best$auc, 2, quantile, 0.25)
cpp.best$auc.high <- apply(cpp.best$auc, 2, quantile, 0.75)
cpp.best$auc.mean <- apply(cpp.best$auc, 2, mean)
cpp.best$auc.sdev <- apply(cpp.best$auc, 2, sd)
cpp.best$x_values <- seq(10, 150, 10)

##################################
label.type <- "ic50"

cpp.best$ic50 <- matrix(nrow = 100, ncol = 0)

for (temp.num in 1:15) {
  load(paste0("~/output_temp/Bortezomib/cpp_var/ic50/bortezomib_cpp_1to100_ic50_var", temp.num, ".RData"))   
  cpp.temp_ic50 <- foreach (temp.ind = 1:100, .combine=rbind) %do% {
    cpp.other_model.all[[temp.ind]]$svmLinear$test_auc
  } 
  cpp.best$ic50 <- cbind(cpp.best$ic50, cpp.temp_ic50)
}

cpp.best$ic50.medium <- apply(cpp.best$ic50, 2, quantile, 0.5)
cpp.best$ic50.low <- apply(cpp.best$ic50, 2, quantile, 0.25)
cpp.best$ic50.high <- apply(cpp.best$ic50, 2, quantile, 0.75)
cpp.best$ic50.mean <- apply(cpp.best$ic50, 2, mean)
cpp.best$ic50.sdev <- apply(cpp.best$ic50, 2, sd)
cpp.best$x_values <- seq(10, 150, 10)

##############################################
cp.best <- list()

load("~/output_temp/Bortezomib/cp_var/bortezomib_cp_1to1_slope_once.RData")
load("~/output_temp/Bortezomib/cp_var/bortezomib_cp_1to1_auc_once.RData")
load("~/output_temp/Bortezomib/cp_var/bortezomib_cp_1to1_ic50_once.RData")
cp.boxplot_matrix <- foreach (temp.ind = 1:1, .combine=rbind) %do% {
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
names(cp.boxplot_matrix) <- c(
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

max(cp.boxplot_matrix)
cp.boxplot_matrix[which(max(cp.boxplot_matrix) == cp.boxplot_matrix)]

cp.best$slope <- cp.boxplot_matrix[which(max(cp.boxplot_matrix) == cp.boxplot_matrix)]
cp.best$auc <- cp.boxplot_matrix[which(max(cp.boxplot_matrix) == cp.boxplot_matrix)]
cp.best$ic50 <- cp.boxplot_matrix[which(max(cp.boxplot_matrix) == cp.boxplot_matrix)]


##################################################
setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Bortezomib/WS")
save(pp.best, cpp.best, cp.best, file = "bor_best_plot.RData")
##################################################
plot_lines <- function(x_values, y_values, sdev = NULL, high = NULL, low = NULL, input_col, input_type = "b", pch = 1) {
  lines(x = x_values, y = y_values, xlim = range(x_values), type = input_type, col = input_col, pch = pch)  
  if (!is.null(sdev)) {
    arrows(x_values, y_values - sdev, x_values, y_values + sdev, length = 0.15, angle = 90, code = 3, col = input_col)
  } else {
    arrows(x_values, low, x_values, high, length = 0.15, angle = 90, code = 3, col = input_col)
  }
}
####################### SLOPE
label.type <- "slope"
setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Bortezomib/Results/")
png(filename = paste0("bor_best_mean_slope.png"), width = 800, height = 800)

par(mar=c(5.1,7,4.5,3.5))
plot(0, 0, xlim = c(10, 150), ylim = c(0.3, 1), type = "n", 
     ylab = "",
     xlab = "# of Patient Samples in Training Data",
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = c(0.3, 1, 14))
title(ylab = "Test AUROC", cex.lab = 1.5, line = 4.5)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = c(0.3, 1, 14))
plot_lines(cpp.best$x_values, cpp.best$slope.mean, cpp.best$slope.sdev, input_col = "blue", pch = 1)
plot_lines(x_values = pp.best$x_values, y_values = pp.best$mean, sdev = pp.best$sdev, input_col = "red", pch = 2)
abline(a = cp.best$slope, b = 0, col = "dark green", lty = 2)
legend("topleft", legend = c("CP2P SVM lin all", "P2P SVM lin all", paste("C2P", names(cp.best$slope))), col = c("blue", "red", "dark green"), pch = c(1,2,45), ncol=3, pt.cex = 1, cex = 1.3)
dev.off()

### median
setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Bortezomib/Results/")
png(filename = paste0("bor_best_median_slope.png"), width = 800, height = 800)

par(mar=c(5.1,7,4.5,3.5))
plot(0, 0, xlim = c(10, 150), ylim = c(0.3, 1), type = "n", 
     ylab = "",
     xlab = "# of Patient Samples in Training Data",
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = c(0.3, 1, 14))
title(ylab = "Test AUC", cex.lab = 1.5, line = 4.5)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = c(0.3, 1, 14))
plot_lines(cpp.best$x_values, cpp.best$slope.medium, low = cpp.best$slope.low, high = cpp.best$slope.high, input_col = "blue", pch = 1)
plot_lines(x_values = pp.best$x_values, y_values = pp.best$medium, low = pp.best$low, high = pp.best$high, input_col = "red", pch = 2)
abline(a = cp.best$slope, b = 0, col = "dark green", lty = 2)
legend("topleft", legend = c("CP2P SVM lin all", "P2P SVM lin all", paste("C2P", names(cp.best$slope))), col = c("blue", "red", "dark green"), pch = c(1,2,45), ncol=3, pt.cex = 1, cex = 1.3)
dev.off()

############# AUC PLOT
setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Bortezomib/Results/")
png(filename = paste0("bor_best_mean_auc.png"), width = 800, height = 800)

par(mar=c(5.1,7,4.5,3.5))
plot(0, 0, xlim = c(10, 150), ylim = c(0.3, 1), type = "n", 
     ylab = "",
     xlab = "# of Patient Samples in Training Data",
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = c(0.3, 1, 14))
title(ylab = "Test AUROC", cex.lab = 1.5, line = 4.5)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = c(0.3, 1, 14))
plot_lines(cpp.best$x_values, cpp.best$auc.mean, cpp.best$auc.sdev, input_col = "blue", pch = 1)
plot_lines(x_values = pp.best$x_values, y_values = pp.best$mean, sdev = pp.best$sdev, input_col = "red", pch = 2)
abline(a = cp.best$auc, b = 0, col = "dark green", lty = 2)
legend("topleft", legend = c("CP2P SVM lin all", "P2P SVM lin all", paste("C2P", names(cp.best$auc))), col = c("blue", "red", "dark green"), pch = c(1,2,45), ncol=3, pt.cex = 1, cex = 1.3)
dev.off()

# medium
setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Bortezomib/Results/")
png(filename = paste0("bor_best_median_auc.png"), width = 800, height = 800)

par(mar=c(5.1,7,4.5,3.5))
plot(0, 0, xlim = c(10, 150), ylim = c(0.3, 1), type = "n", 
     ylab = "",
     xlab = "# of Patient Samples in Training Data",
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = c(0.3, 1, 14))
title(ylab = "Test AUC", cex.lab = 1.5, line = 4.5)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = c(0.3, 1, 14))
plot_lines(cpp.best$x_values, cpp.best$auc.medium, low = cpp.best$auc.low, high = cpp.best$auc.high, input_col = "blue", pch = 1)
plot_lines(x_values = pp.best$x_values, y_values = pp.best$medium, low = pp.best$low, high = pp.best$high, input_col = "red", pch = 2)
abline(a = cp.best$auc, b = 0, col = "dark green", lty = 2)
legend("topleft", legend = c("CP2P SVM lin all", "P2P SVM lin all", paste("C2P", names(cp.best$auc))), col = c("blue", "red", "dark green"), pch = c(1,2,3), ncol=3, pt.cex = 1, cex = 1.3)
dev.off()


############# ic50 plot
label.type <- "ic50"

setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Bortezomib/Results/")
png(filename = paste0("bor_best_mean_ic50.png"), width = 800, height = 800)

par(mar=c(5.1,7,4.5,3.5))
plot(0, 0, xlim = c(10, 150), ylim = c(0.3, 1), type = "n", 
     ylab = "",
     xlab = "# of Patient Samples in Training Data",
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = c(0.3, 1, 14))
title(ylab = "Test AUROC", cex.lab = 1.5, line = 4.5)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = c(0.3, 1, 14))
plot_lines(cpp.best$x_values, cpp.best$ic50.mean, cpp.best$ic50.sdev, input_col = "blue", pch = 1)
plot_lines(x_values = pp.best$x_values, y_values = pp.best$mean, sdev = pp.best$sdev, input_col = "red", pch = 2)
abline(a = cp.best$ic50, b = 0, col = "dark green", lty = 2)
legend("topleft", legend = c("CP2P SVM lin all", "P2P SVM lin all", paste("C2P", names(cp.best$ic50))), col = c("blue", "red", "dark green"), pch = c(1,2,45), ncol=3, pt.cex = 1, cex = 1.3)
dev.off()

# medium
setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Bortezomib/Results/")
png(filename = paste0("bor_best_median_ic50.png"), width = 800, height = 800)

par(mar=c(5.1,7,4.5,3.5))
plot(0, 0, xlim = c(10, 150), ylim = c(0.3, 1), type = "n", 
     ylab = "",
     xlab = "# of Patient Samples in Training Data",
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = c(0.3, 1, 14))
title(ylab = "Test AUC", cex.lab = 1.5, line = 4.5)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = c(0.3, 1, 14))
plot_lines(cpp.best$x_values, cpp.best$ic50.medium, low = cpp.best$ic50.low, high = cpp.best$ic50.high, input_col = "blue", pch = 1)
plot_lines(x_values = pp.best$x_values, y_values = pp.best$medium, low = pp.best$low, high = pp.best$high, input_col = "red", pch = 2)
abline(a = cp.best$ic50, b = 0, col = "dark green", lty = 2)
legend("topleft", legend = c("CP2P SVM lin all", "P2P SVM lin all", paste("C2P", names(cp.best$ic50))), col = c("blue", "red", "dark green"), pch = c(1,2,3), ncol=3, pt.cex = 1, cex = 1.3)
dev.off()

########################### compare IC50, AUC, and slope ################################
plot_lines <- function(x_values, y_values, sdev = NULL, high = NULL, low = NULL, input_col, input_type = "b", pch = 1) {
  lines(x = x_values, y = y_values, xlim = range(x_values), type = input_type, col = input_col, pch = pch)    
}

setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Bortezomib/Results/")
png(filename = paste0("bor_best_mean_labels.png"), width = 800, height = 800)

par(mar=c(5.1,7,4.5,3.5))
plot(0, 0, xlim = c(10, 150), ylim = c(0.55, 0.95), type = "n", 
     ylab = "",
     xlab = "# of Patient Samples in Training Data",
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = c(0.55, 0.85, 6))
title(ylab = "Test AUROC", cex.lab = 1.5, line = 4.5)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = c(0.55, 0.85, 6))
abline(a = cp.best$ic50, b = 0, col = "dark blue", lty = 2)
abline(a = cp.best$auc, b = 0, col = "dark red", lty = 2)
abline(a = cp.best$slope, b = 0, col = "dark green", lty = 2)
plot_lines(cpp.best$x_values, cpp.best$ic50.mean, input_col = "blue", pch = 1)
plot_lines(cpp.best$x_values, cpp.best$auc.mean, input_col = "red", pch = 1)
plot_lines(cpp.best$x_values, cpp.best$slope.mean, input_col = "green", pch = 1)
plot_lines(x_values = pp.best$x_values, y_values = pp.best$mean, input_col = "black", pch = 2)

legend("topleft", legend = c(
  "CP2P-IC50 SVM lin all", 
  "CP2P-AUC SVM lin all", 
  "CP2P-slope SVM lin all", 
  paste("C2P-IC50", names(cp.best$ic50)), 
  paste("C2P-AUC", names(cp.best$auc)), 
  paste("C2P-slope", names(cp.best$slope)),
  "P2P SVM lin all"), 
col = c("blue", "red", "green", "dark blue", "dark red", "dark green", "black"), pch = c(1,1,1,45,45,45,2), ncol=1, pt.cex = 1, cex = 1.3)
dev.off()

