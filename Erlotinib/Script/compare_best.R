library("doParallel")

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/pp_var.RData")
temp.mean <- apply(pp.varying_training_matrix, 2, mean)
which(max(temp.mean) == temp.mean)

temp.median <- apply(pp.varying_training_matrix, 2, median)
which(max(temp.median) == temp.median)

#### get the best patient data ########
pp.best <- list()

pp.best$test_auc <- pp.varying_training_matrix[, which(max(temp.mean) == temp.mean)]
pp.best$x_values <- num_of_patients
pp.best$name <- paste("P2P", names(which(max(temp.mean) == temp.mean)))

#############################################
label.type <- "slope"

setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/")
load(paste0("cpp_var_", label.type, ".RData"))

cpp.best <- list()

temp.mean <- apply(cpp.varying_training_matrix2[20:24, ], 2, mean)
temp.picked <- which(max(temp.mean) == temp.mean)
temp.picked
#temp.picked <- which(max(cpp.varying_training_matrix2[24, ]) == cpp.varying_training_matrix2[24, ])
cpp.best$slope.name <- paste0("CP2P-slope ", names(temp.picked))
cpp.best$slope.test_auc = cpp.varying_training_matrix2[, temp.picked]
cpp.best$x_values <- num_of_patients

#############################################
label.type <- "auc"
rm(temp.mean, cpp.varying_training_matrix2)
load(paste0("cpp_var_", label.type, ".RData"))

temp.mean <- apply(cpp.varying_training_matrix2[20:24, ], 2, mean)
which(max(temp.mean) == temp.mean)
cpp.best$AUC.name <- paste0("CP2P-AUC ", names(which(max(temp.mean) == temp.mean)))
cpp.best$AUC.test_auc = cpp.varying_training_matrix2[, which(max(temp.mean) == temp.mean)]

##############################################
cp.best <- list()

label.type <- "slope"
load(paste0("cp_var_", label.type, ".RData"))
temp.cp_best <- which(cp.varying_training_matrix[10, ] == max(cp.varying_training_matrix[10, ]))
cp.best$slope <- cp.varying_training_matrix[10, temp.cp_best]
cp.best$slope.name <- paste0("C2P-slope ", names(temp.cp_best))

rm(cp.varying_training_matrix)
label.type <- "auc"
load(paste0("cp_var_", label.type, ".RData"))
temp.cp_best <- which(cp.varying_training_matrix[10, ] == max(cp.varying_training_matrix[10, ]))
temp.cp_best
cp.best$AUC.name <- paste0("C2P-AUC ", names(temp.cp_best))
cp.best$AUC <- cp.varying_training_matrix[10, temp.cp_best]

##################################################
setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS")
setwd("C:/Users/Cheng/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/Results/")

#save(pp.best, cpp.best, cp.best, file = "erl_best_plot.RData")
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/erl_best_plot.RData")

##################################################
plot_lines <- function(x_values, y_values, input_col, input_type = "b", pch = 1) {
  lines(x = x_values, y = y_values, xlim = range(x_values), type = input_type, col = input_col, pch = pch)    
}

####################### SLOPE
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/erl_cp_lung_only.RData")
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/erl_cpp_lung_only.RData")

label.type <- "slope"
#setwd("C:/Users/Cheng/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/Results/")
setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/Results/")
png(filename = paste0("erl_best_slope.png"), width = 800, height = 800)

par(mar=c(5.1,7,4.5,3.5))

temp.min <- min(c(cpp.best$slope.test_auc, cp.best$slope, pp.best$test_auc, cp.best$slope.breast_only))
temp.max <- max(c(cpp.best$slope.test_auc, cp.best$slope, pp.best$test_auc))
temp.yaxp = c(floor(temp.min*10)/10, ceiling(temp.max* 10) / 10, (ceiling(temp.max* 10) - floor(temp.min*10)) / 0.5)

plot(0, 0, xlim = c(1, 24), ylim = c(temp.min, temp.max + 0.1), type = "n", 
     ylab = "",
     xlab = "# of Patient Samples in Training Data",
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = temp.yaxp)
title(ylab = "Test AUROC", cex.lab = 1.5, line = 4.5)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = temp.yaxp)

plot_lines(cpp.best$x_values, cpp.best$slope.test_auc, input_col = "blue", pch = 1)
abline(a = cpp.best_lung_only.slope, b = 0, col = "purple", lty = 2)
plot_lines(pp.best$x_values, pp.best$test_auc, input_col = "black", pch = 2)
abline(a = cp.best$slope, b = 0, col = "dark red", lty = 2)
abline(a = cp.best_lung_only.slope, b = 0, col = "dark green", lty = 2)

temp.name_cpp <- paste("CP2P-slope lung only", names(cpp.best_lung_only.slope))
temp.name_cp <- paste("C2P-slope lung only", names(cp.best_lung_only.slope))
legend("topleft", legend = c(cpp.best$slope.name, temp.name_cpp, pp.best$name, cp.best$slope.name, temp.name_cp), 
       col = c("blue", "purple", "black", "red", "dark green"), pch = c(1,45,2,45,45), ncol=1, pt.cex = 1, cex = 1.3)
dev.off()


############# AUC PLOT
png(filename = paste0("erl_best_AUC.png"), width = 800, height = 800)

par(mar=c(5.1,7,4.5,3.5))

temp.min <- min(c(cpp.best$AUC.test_auc, cp.best$AUC, pp.best$test_auc))
temp.max <- max(c(cpp.best$AUC.test_auc, cp.best$AUC, pp.best$test_auc))

temp.yaxp = c(floor(temp.min*10)/10, ceiling(temp.max* 10) / 10, (ceiling(temp.max* 10) - floor(temp.min*10)) / 0.5)

plot(0, 0, xlim = c(1, 24), ylim = c(temp.min, temp.max + 0.1), type = "n", 
     ylab = "",
     xlab = "# of Patient Samples in Training Data",
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = temp.yaxp)
title(ylab = "Test AUROC", cex.lab = 1.5, line = 4.5)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = temp.yaxp)

plot_lines(cpp.best$x_values, cpp.best$AUC.test_auc, input_col = "blue", pch = 1)
abline(a = cpp.best_lung_only.ic50, b = 0, col = "purple", lty = 2)
plot_lines(pp.best$x_values, pp.best$test_auc, input_col = "black", pch = 2)
abline(a = cp.best$AUC, b = 0, col = "dark red", lty = 2)
abline(a = cp.best_lung_only.ic50, b = 0, col = "dark green", lty = 2)

temp.name_cpp <- paste("CP2P-AUC lung only", names(cpp.best_lung_only.ic50))
temp.name_cp <- paste("C2P-AUC lung only", names(cp.best_lung_only.ic50))
legend("topleft", legend = c(cpp.best$AUC.name, temp.name_cpp, pp.best$name, cp.best$AUC.name, temp.name_cp), 
       col = c("blue", "purple", "black", "red", "dark green"), pch = c(1,45,2,45,45), ncol=1, pt.cex = 1, cex = 1.3)
dev.off()

########################### compare IC50, AUC, and slope ################################
plot_lines <- function(x_values, y_values, sdev = NULL, high = NULL, low = NULL, input_col, input_type = "b", pch = 1) {
  lines(x = x_values, y = y_values, xlim = range(x_values), type = input_type, col = input_col, pch = pch)    
}

setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/Results/")
png(filename = paste0("erl_best_mean_labels.png"), width = 800, height = 800)

par(mar=c(5.1,7,4.5,3.5))

temp.min <- min(c(cpp.best$AUC.test_auc, cpp.best$slope.test_auc, cp.best$AUC, pp.best$test_auc))
temp.max <- max(c(cpp.best$AUC.test_auc, cp.best$AUC, pp.best$test_auc))

temp.yaxp = c(floor(temp.min*10)/10, ceiling(temp.max* 10) / 10, (ceiling(temp.max* 10) - floor(temp.min*10)) / 0.5)

plot(0, 0, xlim = c(1, 24), ylim = c(temp.min, temp.max + 0.1), type = "n", 
     ylab = "",
     xlab = "# of Patient Samples in Training Data",
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = temp.yaxp)
title(ylab = "Test AUROC", cex.lab = 1.5, line = 4.5)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = temp.yaxp)

abline(a = cp.best$AUC, b = 0, col = "dark red", lty = 2)
abline(a = cp.best$slope, b = 0, col = "dark green", lty = 2)
#abline(a = 0.78, b = 0, col = "dark green", lty = 2)
plot_lines(cpp.best$x_values, cpp.best$AUC.test_auc, input_col = "red", pch = 1)
plot_lines(cpp.best$x_values, cpp.best$slope.test_auc, input_col = "green", pch = 1)
plot_lines(x_values = pp.best$x_values, y_values = pp.best$test_auc, input_col = "black", pch = 2)

legend("topleft", legend = c(
  
  cpp.best$AUC.name,
  cpp.best$slope.name,
  
  cp.best$AUC.name,
  cp.best$slope.name,  
  "P2P SVM rbf mRMR1000"), 
  col = c("red", "green", "dark red", "dark green", "black"), pch = c(1,1,45,45,2), ncol=1, pt.cex = 1, cex = 1.3)
dev.off()

