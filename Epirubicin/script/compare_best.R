library("doParallel")

#### get the best patient data ########
pp.best <- list()
pp.best$matrix <- matrix(nrow = 100, ncol = 0)
#setwd("E:/AeroFS/temp_output/Epirubicin")
setwd("/home/cheng/output_temp/Epirubicin/")
for (temp.num in 1:9) {
  load(paste0("pp_var/epirubicin_pp_1to100_", temp.num, "_var", temp.num, ".RData"))   
  pp.current_temp <- foreach (temp.ind = 1:100, .combine=rbind) %do% {
    c(pp.other_model.mRMR1000[[temp.ind]]$svmLinear$test_auc)  
  }     
  pp.best$matrix <- cbind(pp.best$matrix, pp.current_temp)
}

pp.best$medium <- apply(pp.best$matrix, 2, quantile, 0.5)
pp.best$low <- apply(pp.best$matrix, 2, quantile, 0.25)
pp.best$high <- apply(pp.best$matrix, 2, quantile, 0.75)
pp.best$mean <- apply(pp.best$matrix, 2, mean)
pp.best$sdev <- apply(pp.best$matrix, 2, sd)

pp.best$x_values <- seq(20, 100, 10)
pp.best$name <- "P2P SVM lin mRMR1000"

#############################################
label.type <- "slope"
cpp.best <- list()
cpp.best$slope <- matrix(nrow = 100, ncol = 0)

for (temp.num in 1:10) {
  load(paste0("~/output_temp/Epirubicin/cpp_var/slope/epirubicin_cpp_1to100_slope_var", temp.num, ".RData"))   
  cpp.temp_slope <- foreach (temp.ind = 1:100, .combine=rbind) %do% {
    cpp.snf.single.mRMR1000[[temp.ind]]$test.auc
  } 
  cpp.best$slope <- cbind(cpp.best$slope, cpp.temp_slope)
}

cpp.best$slope.medium <- apply(cpp.best$slope, 2, quantile, 0.5)
cpp.best$slope.low <- apply(cpp.best$slope, 2, quantile, 0.25)
cpp.best$slope.high <- apply(cpp.best$slope, 2, quantile, 0.75)
cpp.best$slope.mean <- apply(cpp.best$slope, 2, mean)
cpp.best$slope.sdev <- apply(cpp.best$slope, 2, sd)
cpp.best$x_values <- seq(10, 100, 10)

cpp.best$slope.name <- "CP2P-slope SNF mRMR1000"

#############################################
label.type <- "auc"

cpp.best$auc <- matrix(nrow = 100, ncol = 0)

for (temp.num in 1:10) {
  load(paste0("~/output_temp/Epirubicin/cpp_var/auc/epirubicin_cpp_1to100_auc_var", temp.num, ".RData"))   
  cpp.temp_auc <- foreach (temp.ind = 1:100, .combine=rbind) %do% {
    cpp.snf.single.mRMR1000[[temp.ind]]$test.auc
  } 
  cpp.best$auc <- cbind(cpp.best$auc, cpp.temp_auc)
}

cpp.best$auc.medium <- apply(cpp.best$auc, 2, quantile, 0.5)
cpp.best$auc.low <- apply(cpp.best$auc, 2, quantile, 0.25)
cpp.best$auc.high <- apply(cpp.best$auc, 2, quantile, 0.75)
cpp.best$auc.mean <- apply(cpp.best$auc, 2, mean)
cpp.best$auc.sdev <- apply(cpp.best$auc, 2, sd)
cpp.best$auc.name <- "CP2P-AUC SNF mRMR1000"

##############################################
cp.best <- list()

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Epirubicin/WS/cp_var_auc.RData")
cp.best$auc <- max(cp.varying_training_matrix[4, ])
cp.best$auc.name <- paste("C2P-AUC", names(which(max(cp.varying_training_matrix[4, ]) == cp.varying_training_matrix[4, ])))

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Epirubicin/WS/cp_var_slope.RData")
cp.best$slope <- max(cp.varying_training_matrix[4, ])
cp.best$slope.name <- paste("C2P-slope", names(which(max(cp.varying_training_matrix[4, ]) == cp.varying_training_matrix[4, ])))

##################################################
setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Epirubicin/WS")
save(pp.best, cpp.best, cp.best, file = "epi_best_plot.RData")
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
setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Epirubicin/Results/")
png(filename = paste0("epi_best_mean_slope.png"), width = 800, height = 800)

temp.min <- min(c(cpp.best$slope.mean - cpp.best$slope.sdev, cp.best$slope, pp.best$mean - pp.best$sdev))
temp.max <- max(c(cpp.best$slope.mean + cpp.best$slope.sdev, cp.best$slope, pp.best$mean + pp.best$sdev))
temp.yaxp = c(floor(temp.min*10)/10, ceiling(temp.max* 10) / 10, (ceiling(temp.max* 10) - floor(temp.min*10)) / 0.5)

par(mar=c(5.1,7,4.5,3.5))

plot(0, 0, xlim = range(cpp.best$x_values), ylim = c(temp.min, temp.max + 0.1), type = "n", 
     ylab = "",
     xlab = "# of Patient Samples in Training Data",
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = temp.yaxp)
title(ylab = "Test AUROC", cex.lab = 1.5, line = 4.5)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = temp.yaxp)

plot_lines(cpp.best$x_values, cpp.best$slope.mean, cpp.best$slope.sdev, input_col = "blue", pch = 1)
plot_lines(pp.best$x_values, pp.best$mean, pp.best$sdev, input_col = "red", pch = 2)
abline(a = cp.best$slope, b = 0, col = "dark green", lty = 2)
legend("topleft", legend = c(cpp.best$slope.name, pp.best$name, cp.best$slope.name), col = c("blue", "red", "dark green"), pch = c(1,2,45), ncol=1, pt.cex = 1, cex = 1.3)
dev.off()

############# AUC PLOT
setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Epirubicin/Results/")
png(filename = paste0("epi_best_mean_auc.png"), width = 800, height = 800)

temp.min <- min(c(cpp.best$auc.mean - cpp.best$auc.sdev, cp.best$auc, pp.best$mean - pp.best$sdev))
temp.max <- max(c(cpp.best$auc.mean + cpp.best$auc.sdev, cp.best$auc, pp.best$mean + pp.best$sdev))
temp.yaxp = c(floor(temp.min*10)/10, ceiling(temp.max* 10) / 10, (ceiling(temp.max* 10) - floor(temp.min*10)) / 0.5)

par(mar=c(5.1,7,4.5,3.5))

plot(0, 0, xlim = range(cpp.best$x_values), ylim = c(temp.min, temp.max + 0.1), type = "n", 
     ylab = "",
     xlab = "# of Patient Samples in Training Data",
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = temp.yaxp)
title(ylab = "Test AUROC", cex.lab = 1.5, line = 4.5)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = temp.yaxp)

plot_lines(cpp.best$x_values, cpp.best$auc.mean, cpp.best$auc.sdev, input_col = "blue", pch = 1)
plot_lines(pp.best$x_values, pp.best$mean, pp.best$sdev, input_col = "red", pch = 2)
abline(a = cp.best$auc, b = 0, col = "dark green", lty = 2)
legend("topleft", legend = c(cpp.best$auc.name, pp.best$name, cp.best$auc.name), col = c("blue", "red", "dark green"), pch = c(1,2,45), ncol=1, pt.cex = 1, cex = 1.3)
dev.off()


########################### compare IC50, AUC, and slope ################################
plot_lines <- function(x_values, y_values, sdev = NULL, high = NULL, low = NULL, input_col, input_type = "b", pch = 1) {
  lines(x = x_values, y = y_values, xlim = range(x_values), type = input_type, col = input_col, pch = pch)    
}

setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Epirubicin/Results/")
png(filename = paste0("epi_best_mean_labels.png"), width = 800, height = 800)

par(mar=c(5.1,7,4.5,3.5))

temp.min <- min(c(cpp.best$auc.mean, cpp.best$slope.mean, cp.best$auc, pp.best$mean))
temp.max <- max(c(cpp.best$auc.mean, cpp.best$slope.mean, cp.best$auc, pp.best$mean))

temp.yaxp = c(floor(temp.min*10)/10, ceiling(temp.max* 10) / 10, (ceiling(temp.max* 10) - floor(temp.min*10)) / 0.5)

plot(0, 0, xlim = c(10, 100), ylim = c(temp.min, temp.max + 0.1), type = "n", 
     ylab = "",
     xlab = "# of Patient Samples in Training Data",
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = temp.yaxp)
title(ylab = "Test AUROC", cex.lab = 1.5, line = 4.5)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = temp.yaxp)

abline(a = cp.best$auc, b = 0, col = "dark red", lty = 2)
abline(a = cp.best$slope, b = 0, col = "dark green", lty = 2)
#abline(a = 0.78, b = 0, col = "dark green", lty = 2)
plot_lines(cpp.best$x_values, cpp.best$auc.mean, input_col = "red", pch = 1)
plot_lines(cpp.best$x_values, cpp.best$slope.mean, input_col = "green", pch = 1)
plot_lines(x_values = pp.best$x_values, y_values = pp.best$mean, input_col = "black", pch = 2)

legend("topleft", legend = c(
  
  cpp.best$auc.name,
  cpp.best$slope.name,
  
  cp.best$auc.name,
  cp.best$slope.name,  
  pp.best$name), 
  col = c("red", "green", "dark red", "dark green", "black"), pch = c(1,1,45,45,2), ncol=1, pt.cex = 1, cex = 1.3)
dev.off()
