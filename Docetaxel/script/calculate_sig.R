load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/cp_var_ic50.RData")
which(cp.varying_training_matrix[19, ] == max(cp.varying_training_matrix[19, ]))
max(which(cp.varying_training_matrix[19, ] == max(cp.varying_training_matrix[19, ])))

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/cp_var_auc.RData")
which(cp.varying_training_matrix[19, ] == max(cp.varying_training_matrix[19, ]))
max(cp.varying_training_matrix[19, ])

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/cp_var_slope.RData")
which(cp.varying_training_matrix[19, ] == max(cp.varying_training_matrix[19, ]))
max(cp.varying_training_matrix[19, ])

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/cpp_var_ic50.RData")
which(cpp.varying_training_matrix[23, ] == max(cpp.varying_training_matrix[23, ]))
max(cpp.varying_training_matrix[23, ])


load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/cpp_var_auc.RData")
which(cpp.varying_training_matrix[23, ] == max(cpp.varying_training_matrix[23, ]))
max(cpp.varying_training_matrix[23, ])

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/cpp_var_slope.RData")
which(cpp.varying_training_matrix[23, ] == max(cpp.varying_training_matrix[23, ]))
max(cpp.varying_training_matrix[23, ])

### compare CP2P to P2P
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/pp_var.RData")
temp.pp_best <- which(pp.varying_training_matrix[10, ] == max(pp.varying_training_matrix[10, ]))
temp.pp_best.auc <- pp.varying_training_matrix[, temp.pp_best]

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/cpp_var_ic50.RData")
temp.cpp_best <- which(cpp.varying_training_matrix[23, ] == max(cpp.varying_training_matrix[23, ]))
temp.cpp_best.auc <- cpp.varying_training_matrix[14:23, temp.cpp_best]
wilcox.test(temp.cpp_best.auc, temp.pp_best.auc, alternative = "g", paired=TRUE)

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/cpp_var_auc.RData")
temp.cpp_best <- which(cpp.varying_training_matrix[23, ] == max(cpp.varying_training_matrix[23, ]))
temp.cpp_best.auc <- cpp.varying_training_matrix[14:23, temp.cpp_best]
wilcox.test(temp.cpp_best.auc, temp.pp_best.auc, alternative = "g", paired = TRUE)

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/cpp_var_slope.RData")
temp.cpp_best <- which(cpp.varying_training_matrix[23, ] == max(cpp.varying_training_matrix[23, ]))
temp.cpp_best.auc <- cpp.varying_training_matrix[14:23, temp.cpp_best]
wilcox.test(temp.cpp_best.auc, temp.pp_best.auc, alternative = "g", paired = TRUE)

### compare CP2P
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/cpp_var_ic50.RData")
temp.cpp_best <- which(cpp.varying_training_matrix[23, ] == max(cpp.varying_training_matrix[23, ]))
temp.cpp_best.ic50 <- cpp.varying_training_matrix[, temp.cpp_best]

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/cpp_var_auc.RData")
temp.cpp_best <- which(cpp.varying_training_matrix[23, ] == max(cpp.varying_training_matrix[23, ]))
temp.cpp_best.auc <- cpp.varying_training_matrix[, temp.cpp_best]

range(temp.cpp_best.auc)
range(temp.cpp_best.ic50)
wilcox.test(temp.cpp_best.auc, temp.cpp_best.ic50, alternative = "t", paired=TRUE)

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/cpp_var_slope.RData")
temp.cpp_best.slope <- cpp.varying_training_matrix[, temp.cpp_best]
wilcox.test(temp.cpp_best.auc, temp.cpp_best.slope, alternative = "g", paired=TRUE)


### C2P
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/cp_var_ic50.RData")
temp.cp_best <- which(cp.varying_training_matrix[19, ] == max(cp.varying_training_matrix[19, ]))
temp.cp_best.ic50 <- cp.varying_training_matrix[, temp.cp_best[1]]

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/cp_var_auc.RData")
temp.cp_best <- which(cp.varying_training_matrix[19, ] == max(cp.varying_training_matrix[19, ]))
temp.cp_best.auc <- cp.varying_training_matrix[, temp.cp_best]

range(temp.cp_best.auc)
range(temp.cp_best.ic50)
wilcox.test(temp.cp_best.ic50, temp.cp_best.auc, alternative = "g", paired = TRUE)


