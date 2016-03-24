load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/pp_var.RData")
which(pp.varying_training_matrix[12, ] == max(pp.varying_training_matrix[12, ]))
max(pp.varying_training_matrix[12, ])

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/cp_var_auc.RData")
which(cp.varying_training_matrix[10, ] == max(cp.varying_training_matrix[10, ]))
max(cp.varying_training_matrix[10, ])

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/cp_var_slope.RData")
which(cp.varying_training_matrix[10, ] == max(cp.varying_training_matrix[10, ]))
max(cp.varying_training_matrix[10, ])

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/cpp_var_auc.RData")
which(cpp.varying_training_matrix2[24, ] == max(cpp.varying_training_matrix2[24, ]))
max(cpp.varying_training_matrix2[24, ])

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/cpp_var_slope.RData")
which(cpp.varying_training_matrix2[24, ] == max(cpp.varying_training_matrix2[24, ]))
max(cpp.varying_training_matrix2[24, ])

### compare CP2P to P2P
wilcox.test(cpp.best$AUC.test_auc[13:24], pp.best$test_auc, paired = TRUE)

wilcox.test(cpp.best$slope.test_auc[13:24], pp.best$test_auc, paired = TRUE)

wilcox.test(cpp.best$slope.test_auc[13:24], cpp.best$AUC.test_auc[13:24], paired = TRUE)


load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/cp_var_auc.RData")
temp.best_ind <- which(cp.varying_training_matrix[10, ] == max(cp.varying_training_matrix[10, ]))
temp.cp.best_auc <- cp.varying_training_matrix[, temp.best_ind]

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/cp_var_slope.RData")
temp.best_ind <- which(cp.varying_training_matrix[10, ] == max(cp.varying_training_matrix[10, ]))
temp.cp.best_slope <- cp.varying_training_matrix[, temp.best_ind]

wilcox.test(temp.cp.best_auc, temp.cp.best_slope, paired = TRUE)


### C2P
