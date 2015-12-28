library("doParallel")

setwd("/Users/chengzhao/Google Drive (chengzhao41@gmail.com)/DrugResponse Paper Results (R Workspace)/")

#### get the best patient data ########
load("Bortezomib/pp_var/bortezomib_pp_1to100_10_var10.RData")

pp.current_temp <- foreach (temp.ind = 1:100, .combine=rbind) %do% {
  c(pp.other_model.all[[temp.ind]]$svmLinear$test_auc)  
}

pp.best.mean <- mean(pp.current_temp)

#############################################
label.type <- "slope"
cp.best <- list()
cp.best$slope <- matrix(nrow = 100, ncol = 0)

for (temp.num in 3:30) {
  load(paste0("Bortezomib/cp_var/slope/bortezomib_cp_1to100_slope_var", temp.num, ".RData"))   
  cp.temp_slope <- foreach (temp.ind = 1:100, .combine=rbind) %do% {
    cp.other_model.all[[temp.ind]]$svmLinear$test_auc
  } 
  cp.best$slope <- cbind(cp.best$slope, cp.temp_slope)
}

cp.best$slope.medium <- apply(cp.best$slope, 2, quantile, 0.5)
cp.best$slope.low <- apply(cp.best$slope, 2, quantile, 0.25)
cp.best$slope.high <- apply(cp.best$slope, 2, quantile, 0.75)
cp.best$slope.mean <- apply(cp.best$slope, 2, mean)
cp.best$slope.sdev <- apply(cp.best$slope, 2, sd)
cp.best$x_values <- seq(30, 300, 10)

#############################################

#setwd("/Users/chengzhao/Dropbox/output_WS/Bortezomib/")

label.type <- "slope_p100"
cpp.best <- list()
cpp.best$slope <- matrix(nrow = 100, ncol = 0)

for (temp.num in 1:30) {
  load(paste0("bortezomib_cpp_1to100_", label.type, "_var", temp.num, ".RData"))
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
cpp.best$x_values <- seq(10, 300, 10)


##############################################
setwd("/Users/chengzhao/Git/CP2P/Bortezomib/output_WS/")
save(pp.best.mean, cpp.best, cp.best, file = "bor_varC_best_plot.RData")

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
setwd("/Users/chengzhao/Git/CP2P/Bortezomib/plots/")
png(filename = paste0("bor_best_mean_slope.png"), width = 800, height = 800)

par(mar=c(5.1,7,4.5,3.5))
plot(0, 0, xlim = c(10, 300), ylim = c(0.4, 0.9), type = "n", 
     ylab = "",
     xlab = "# of Cell Line Samples in Training Data",
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = c(0.4, 0.9, 10))
title(ylab = "Test AUROC", cex.lab = 1.5, line = 4.5)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = c(0.4, 0.9, 10))
plot_lines(x_values = cpp.best$x_values, y_values = cpp.best$slope.mean, sdev = cpp.best$slope.sdev, input_col = "blue", pch = 1)
plot_lines(x_values = cp.best$x_values, y_values = cp.best$slope.mean, sdev = cp.best$slope.sdev, input_col = "dark green", pch = 2)
abline(a = pp.best.mean, b = 0, col = "red", lty = 2)
legend("topleft", legend = c("CP2P SVM lin all", "C2P SVM lin all", paste("P2P", "P2P SVM lin all")), col = c("blue", "dark green", "red"), pch = c(1,2,45), ncol=3, pt.cex = 1, cex = 1.3)
dev.off()


