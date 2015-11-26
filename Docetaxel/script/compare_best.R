library("doParallel")

load("C:/Users/Cheng/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/pp_var.RData")
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/pp_var.RData")
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
#setwd("C:/Users/Cheng/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/")
setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS")
load(paste0("cpp_var_", label.type, ".RData"))

cpp.best <- list()

temp.mean <- apply(cpp.varying_training_matrix[19:23, ], 2, mean)
temp.best <- which(max(temp.mean) == temp.mean)

cpp.best$slope.name <- paste0("CP2P-slope ", names(temp.best))
cpp.best$slope.test_auc = cpp.varying_training_matrix[,temp.best]
cpp.best$x_values <- num_of_patients

#############################################
label.type <- "auc"
rm(temp.mean, cpp.varying_training_matrix)
load(paste0("cpp_var_", label.type, ".RData"))

temp.mean <- apply(cpp.varying_training_matrix[19:23, ], 2, mean)
which(max(temp.mean) == temp.mean)
cpp.best$AUC.name <- paste0("CP2P-AUC ", names(which(max(temp.mean) == temp.mean)))
cpp.best$AUC.test_auc = cpp.varying_training_matrix[, which(max(temp.mean) == temp.mean)]

##################################
label.type <- "ic50"
rm(temp.mean, cpp.varying_training_matrix)
load(paste0("cpp_var_", label.type, ".RData"))

temp.mean <- apply(cpp.varying_training_matrix[19:23, ], 2, mean)
which(max(temp.mean) == temp.mean)
cpp.best$IC50.name <- paste0("CP2P-IC50 ", names(which(max(temp.mean) == temp.mean)))
cpp.best$IC50.test_auc = cpp.varying_training_matrix[, which(max(temp.mean) == temp.mean)]

##############################################
cp.best <- list()

label.type <- "slope"
load(paste0("cp_var_", label.type, ".RData"))
temp.cp_best <- which(cp.varying_training_matrix[19, ] == max(cp.varying_training_matrix[19, ]))
cp.best$slope <- cp.varying_training_matrix[19, temp.cp_best]
cp.best$slope.name <- paste0("C2P-slope ", names(temp.cp_best))

label.type <- "auc"
load(paste0("cp_var_", label.type, ".RData"))
temp.cp_best <- which(cp.varying_training_matrix[19, ] == max(cp.varying_training_matrix[19, ]))
temp.cp_best
cp.best$AUC.name <- paste0("C2P-AUC ", names(temp.cp_best))
cp.best$AUC <- cp.varying_training_matrix[19, temp.cp_best]

label.type <- "ic50"
load(paste0("cp_var_", label.type, ".RData"))
temp.cp_best <- which(cp.varying_training_matrix[19, ] == max(cp.varying_training_matrix[19, ]))
temp.cp_best
cp.best$IC50.name <- paste0("C2P-IC50 ", names(temp.cp_best[1]))
cp.best$IC50 <- cp.varying_training_matrix[19, temp.cp_best[1]]

##################################################
setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS")
#setwd("C:/Users/Cheng/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/Results/")

#save(pp.best, cpp.best, cp.best, file = "doc_best_plot.RData")
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/doc_best_plot.RData")

doc_cp_ic50_breast <- read.csv("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Paper/Table/doc_cp_ic50_breast.csv")
max(doc_cp_ic50_breast[, 2])
cp.best$IC50.breast_only <- max(doc_cp_ic50_breast[, 2])
cp.best$IC50.breast_only.name <- paste("C2P breast only", doc_cp_ic50_breast[which(max(doc_cp_ic50_breast[, 2]) == doc_cp_ic50_breast[, 2]), 1])[1]

doc_cp_auc_breast <- read.csv("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Paper/Table/doc_cp_auc_breast.csv")
max(doc_cp_auc_breast[, 2])
cp.best$auc.breast_only <- max(doc_cp_auc_breast[, 2])
cp.best$auc.breast_only.name <- paste("C2P breast only", doc_cp_auc_breast[which(max(doc_cp_auc_breast[, 2]) == doc_cp_auc_breast[, 2]), 1])[1]

doc_cp_slope_breast <- read.csv("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Paper/Table/doc_cp_slope_breast.csv")
max(doc_cp_slope_breast[, 2])
cp.best$slope.breast_only <- max(doc_cp_slope_breast[, 2])
cp.best$slope.breast_only.name <- paste("C2P breast only", doc_cp_slope_breast[which(max(doc_cp_slope_breast[, 2]) == doc_cp_slope_breast[, 2]), 1])[1]

doc_cpp_ic50_breast <- read.csv("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Paper/Table/doc_cpp_ic50_breast.csv")
max(doc_cpp_ic50_breast[, 2])
cpp.best$IC50.breast_only <- max(doc_cpp_ic50_breast[, 2])
cpp.best$IC50.breast_only.name <- paste("CP2P breast only", doc_cpp_ic50_breast[which(max(doc_cpp_ic50_breast[, 2]) == doc_cpp_ic50_breast[, 2]), 1])[1]

doc_cpp_auc_breast <- read.csv("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Paper/Table/doc_cpp_auc_breast.csv")
max(doc_cpp_auc_breast[, 2])
cpp.best$auc.breast_only <- max(doc_cpp_auc_breast[, 2])
cpp.best$auc.breast_only.name <- paste("CP2P breast only", doc_cpp_auc_breast[which(max(doc_cpp_auc_breast[, 2]) == doc_cpp_auc_breast[, 2]), 1])[1]

doc_cpp_slope_breast <- read.csv("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Paper/Table/doc_cpp_slope_breast.csv")
max(doc_cpp_slope_breast[, 2])
cpp.best$slope.breast_only <- max(doc_cpp_slope_breast[, 2])
cpp.best$slope.breast_only.name <- paste("CP2P breast only", doc_cpp_slope_breast[which(max(doc_cpp_slope_breast[, 2]) == doc_cpp_slope_breast[, 2]), 1])[1]

save(pp.best, cpp.best, cp.best, file = "doc_best_plot.RData")

##################################################
plot_lines <- function(x_values, y_values, input_col, input_type = "b", pch = 1) {
  lines(x = x_values, y = y_values, xlim = range(x_values), type = input_type, col = input_col, pch = pch)    
}
####################### SLOPE
label.type <- "slope"
#setwd("C:/Users/Cheng/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/Results/")
setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/Results/")
png(filename = paste0("doc_best_slope.png"), width = 800, height = 800)

par(mar=c(5.1,7,4.5,3.5))

temp.min <- min(c(cpp.best$slope.test_auc, cp.best$slope, pp.best$test_auc, cp.best$slope.breast_only))
temp.max <- max(c(cpp.best$slope.test_auc, cp.best$slope, pp.best$test_auc))
temp.yaxp = c(floor(temp.min*10)/10, ceiling(temp.max* 10) / 10, (ceiling(temp.max* 10) - floor(temp.min*10)) / 0.5)

plot(0, 0, xlim = c(1, 23), ylim = c(temp.min, temp.max + 0.1), type = "n", 
     ylab = "",
     xlab = "# of Patient Samples in Training Data",
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = temp.yaxp)
title(ylab = "Test AUROC", cex.lab = 1.5, line = 4.5)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = temp.yaxp)

plot_lines(cpp.best$x_values, cpp.best$slope.test_auc, input_col = "blue", pch = 1)
abline(a = cpp.best$slope.breast_only, b = 0, col = "purple", lty = 2)
plot_lines(pp.best$x_values, pp.best$test_auc, input_col = "black", pch = 2)
abline(a = cp.best$slope, b = 0, col = "dark red", lty = 2)
abline(a = cp.best$slope.breast_only, b = 0, col = "dark green", lty = 2)
legend("topleft", legend = c(cpp.best$slope.name, cpp.best$slope.breast_only.name, pp.best$name, cp.best$slope.name, cp.best$slope.breast_only.name), 
       col = c("blue", "purple", "black", "red", "dark green"), pch = c(1,45,2,45,45), ncol=1, pt.cex = 1, cex = 1.3)
dev.off()

############# AUC PLOT
png(filename = paste0("doc_best_AUC.png"), width = 800, height = 800)

par(mar=c(5.1,7,4.5,3.5))

temp.min <- min(c(cpp.best$AUC.test_auc, cp.best$AUC, pp.best$test_auc))
temp.max <- max(c(cpp.best$AUC.test_auc, cp.best$AUC, pp.best$test_auc))

temp.yaxp = c(floor(temp.min*10)/10, ceiling(temp.max* 10) / 10, (ceiling(temp.max* 10) - floor(temp.min*10)) / 0.5)

plot(0, 0, xlim = c(1, 23), ylim = c(temp.min, temp.max + 0.1), type = "n", 
     ylab = "",
     xlab = "# of Patient Samples in Training Data",
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = temp.yaxp)
title(ylab = "Test AUROC", cex.lab = 1.5, line = 4.5)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = temp.yaxp)

plot_lines(cpp.best$x_values, cpp.best$AUC.test_auc, input_col = "blue", pch = 1)
abline(a = cpp.best$auc.breast_only, b = 0, col = "purple", lty = 2)
plot_lines(pp.best$x_values, pp.best$test_auc, input_col = "black", pch = 2)
abline(a = cp.best$AUC, b = 0, col = "dark red", lty = 2)
abline(a = cp.best$auc.breast_only, b = 0, col = "dark green", lty = 2)
legend("topleft", legend = c(cpp.best$AUC.name, cpp.best$auc.breast_only.name, pp.best$name, cp.best$AUC.name, cp.best$auc.breast_only.name), 
       col = c("blue", "purple", "black", "red", "dark green"), pch = c(1,45,2,45,45), ncol=1, pt.cex = 1, cex = 1.3)
dev.off()

############# ic50 plot
png(filename = paste0("doc_best_IC50.png"), width = 800, height = 800)

par(mar=c(5.1,7,4.5,3.5))

temp.min <- min(c(cpp.best$IC50.test_auc, cp.best$IC50, pp.best$test_auc))
temp.max <- max(c(cpp.best$IC50.test_auc, cp.best$IC50, pp.best$test_auc))

temp.yaxp = c(floor(temp.min*10)/10, ceiling(temp.max* 10) / 10, (ceiling(temp.max* 10) - floor(temp.min*10)) / 0.5)

plot(0, 0, xlim = c(1, 23), ylim = c(temp.min, temp.max + 0.1), type = "n", 
     ylab = "",
     xlab = "# of Patient Samples in Training Data",
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = temp.yaxp)
title(ylab = "Test AUROC", cex.lab = 1.5, line = 4.5)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = temp.yaxp)

plot_lines(cpp.best$x_values, cpp.best$IC50.test_auc, input_col = "blue", pch = 1)
abline(a = cpp.best$IC50.breast_only, b = 0, col = "purple", lty = 2)
plot_lines(pp.best$x_values, pp.best$test_auc, input_col = "black", pch = 2)
abline(a = cp.best$IC50, b = 0, col = "dark red", lty = 2)
abline(a = cp.best$IC50.breast_only, b = 0, col = "dark green", lty = 2)
legend("topleft", legend = c(cpp.best$IC50.name, cpp.best$IC50.breast_only.name, pp.best$name, cp.best$IC50.name, cp.best$IC50.breast_only.name), 
       col = c("blue", "purple", "black", "red", "dark green"), pch = c(1,45,2,45,45), ncol=1, pt.cex = 1, cex = 1.3)
dev.off()

########################### compare IC50, AUC, and slope ################################
plot_lines <- function(x_values, y_values, sdev = NULL, high = NULL, low = NULL, input_col, input_type = "b", pch = 1) {
  lines(x = x_values, y = y_values, xlim = range(x_values), type = input_type, col = input_col, pch = pch)    
}

setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/Results/")
png(filename = paste0("doc_best_mean_labels.png"), width = 800, height = 800)

par(mar=c(5.1,7,4.5,3.5))
plot(0, 0, xlim = c(1, 23), ylim = c(0.65, 1.05), type = "n", 
     ylab = "",
     xlab = "# of Patient Samples in Training Data",
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = c(0.65, 0.95, 6))
title(ylab = "Test AUROC", cex.lab = 1.5, line = 4.5)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = c(0.65, 0.95, 6))
abline(a = cp.best$IC50, b = 0, col = "dark blue", lty = 2)
abline(a = cp.best$AUC, b = 0, col = "dark red", lty = 2)
abline(a = cp.best$slope, b = 0, col = "dark green", lty = 2)
plot_lines(cpp.best$x_values, cpp.best$IC50.test_auc, input_col = "blue", pch = 1)
plot_lines(cpp.best$x_values, cpp.best$AUC.test_auc, input_col = "red", pch = 1)
plot_lines(cpp.best$x_values, cpp.best$slope.test_auc, input_col = "green", pch = 1)
plot_lines(x_values = pp.best$x_values, y_values = pp.best$test_auc, input_col = "black", pch = 2)

legend("topleft", legend = c(
  cpp.best$IC50.name,
  cpp.best$AUC.name,
  cpp.best$slope.name,
  cp.best$IC50.name,
  cp.best$AUC.name,
  cp.best$slope.name,  
  "P2P SVM rbf mRMR1000"), 
  col = c("blue", "red", "green", "dark blue", "dark red", "dark green", "black"), pch = c(1,1,1,45,45,45,2), ncol=1, pt.cex = 1, cex = 1.3)
dev.off()

# ########################### compare IC50, AUC, and slope ################################
# plot_lines <- function(x_values, y_values, sdev = NULL, high = NULL, low = NULL, input_col, input_type = "b", pch = 1) {
#   lines(x = x_values, y = y_values, xlim = range(x_values), type = input_type, col = input_col, pch = pch)    
# }
# 
# setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/Results/")
# png(filename = paste0("erl_best_mean_labels.png"), width = 800, height = 800)
# 
# par(mar=c(5.1,7,4.5,3.5))
# 
# temp.min <- min(c(cpp.best$AUC.test_auc, cpp.best$slope.test_auc, cp.best$AUC, pp.best$test_auc))
# temp.max <- max(c(cpp.best$AUC.test_auc, cp.best$AUC, pp.best$test_auc))
# 
# temp.yaxp = c(floor(temp.min*10)/10, ceiling(temp.max* 10) / 10, (ceiling(temp.max* 10) - floor(temp.min*10)) / 0.5)
# 
# plot(0, 0, xlim = c(1, 24), ylim = c(temp.min, temp.max + 0.1), type = "n", 
#      ylab = "",
#      xlab = "# of Patient Samples in Training Data",
#      las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = temp.yaxp)
# title(ylab = "Test AUC", cex.lab = 1.5, line = 4.5)
# axis(side = 4, las = 1, cex.axis = 1.2, yaxp = temp.yaxp)
# 
# abline(a = cp.best$AUC, b = 0, col = "dark red", lty = 2)
# abline(a = cp.best$slope, b = 0, col = "dark green", lty = 2)
# #abline(a = 0.78, b = 0, col = "dark green", lty = 2)
# plot_lines(cpp.best$x_values, cpp.best$AUC.test_auc, input_col = "red", pch = 1)
# plot_lines(cpp.best$x_values, cpp.best$slope.test_auc, input_col = "green", pch = 1)
# plot_lines(x_values = pp.best$x_values, y_values = pp.best$test_auc, input_col = "black", pch = 2)
# 
# legend("topleft", legend = c(
#   
#   cpp.best$AUC.name,
#   cpp.best$slope.name,
#   
#   cp.best$AUC.name,
#   cp.best$slope.name,  
#   "P2P SVM rbf mRMR1000"), 
#   col = c("red", "green", "dark red", "dark green", "black"), pch = c(1,1,3,3,2), ncol=1, pt.cex = 1, cex = 1.3)
# dev.off()
