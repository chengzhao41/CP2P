load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Epirubicin/WS/epi_best_plot.RData")

# AUC
for (temp.alt_ind in 1:9) {    
  temp.wtest <- wilcox.test(cpp.best$auc[, temp.alt_ind + 1], pp.best$matrix[, temp.alt_ind], paired = FALSE, alternative = "great")
  print(temp.wtest$p.value)  
}

# slope
for (temp.alt_ind in 1:9) {    
  temp.wtest <- wilcox.test(cpp.best$slope[, temp.alt_ind + 1], pp.best$matrix[, temp.alt_ind], paired = FALSE, alternative = "great")
  print(temp.wtest$p.value)  
}

# AUC vs. slope
for (temp.alt_ind in 1:10) {    
  temp.wtest <- wilcox.test(cpp.best$auc[, temp.alt_ind], cpp.best$slope[, temp.alt_ind], paired = TRUE, alternative = "t")
  print(temp.wtest$p.value)  
}

temp.wtest <- wilcox.test(cpp.best$auc[, 6], cpp.best$slope[, 6], paired = TRUE, alternative = "g")
print(temp.wtest$p.value)  

## C2P AUC vs. slope
# AUC vs. slope


load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Epirubicin/WS/cp_var_auc.RData")
cp.best_auc <- cp.varying_training_matrix
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Epirubicin/WS/cp_var_slope.RData")
cp.best_slope <- cp.varying_training_matrix

for (temp.alt_ind in 1:4) {    
  temp.wtest <- wilcox.test(cp.best_auc, cp.best_slope, paired = TRUE, alternative = "t")
  print(temp.wtest$p.value)  
}

