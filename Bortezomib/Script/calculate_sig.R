library("doParallel")

pp.test_auc <- list()

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
  
  colnames(pp.boxplot_matrix) <- c(
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
  pp.test_auc[[temp.num]] <- pp.boxplot_matrix
}

temp.order <- order(apply(pp.test_auc[[15]], 2, median), decreasing = TRUE)
best1 <- pp.test_auc[[15]][, temp.order[1]]
best1 <- c(best1, pp.test_auc[[14]][, temp.order[1]])

for (temp.alt_ind in 2:21) {
  best2 <- pp.test_auc[[15]][, temp.order[temp.alt_ind]]
  best2 <- c(best2, pp.test_auc[[14]][, temp.order[temp.alt_ind]])
  temp.wtest <- wilcox.test(best1, best2, paired = TRUE, alternative = "great")
  print(temp.wtest$p.value)
  stopifnot(temp.wtest$p.value < 0.05)
}

####################################### C2P ##############################
cp.best_300 <- list()
cp.varying_training_matrix <- matrix(nrow = 0, ncol = 21)
label.types <- c("ic50", "auc", "slope")

for (temp.label_ind in 1:3) {
  temp.num = 30
  load(paste0("~/output_temp/Bortezomib/cp_var/", label.types[temp.label_ind], "/bortezomib_cp_1to100_", label.types[temp.label_ind], "_var", temp.num, ".RData"))   
  # Boxplot Data ---------------------------------
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
  
  colnames(cp.boxplot_matrix) <- c(
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
  
  if (temp.label_ind == 1) {
    cp.best_300$ic50 <- cp.boxplot_matrix   
  } else if (temp.label_ind == 2) {
    cp.best_300$auc <- cp.boxplot_matrix     
  } else if (temp.label_ind == 3) {  
    cp.best_300$slope <- cp.boxplot_matrix   
  } else {
    stop("wrong temp.label_ind")
  }
}

# IC50
temp.order <- order(apply(cp.best_300$ic50, 2, median), decreasing = TRUE)
best1 <- cp.best_300$ic50[, temp.order[1]]
for (temp.alt_ind in 2:21) {  
  best2 <- cp.best_300$ic50[, temp.order[temp.alt_ind]]    
  temp.wtest <- wilcox.test(best1, best2, paired = TRUE, alternative = "great")
  print(temp.wtest$p.value)
  stopifnot(temp.wtest$p.value < 0.05)  
}

# AUC
temp.order <- order(apply(cp.best_300$auc, 2, median), decreasing = TRUE)
best1 <- cp.best_300$auc[, temp.order[1]]

for (temp.alt_ind in 2:21) {  
  best2 <- cp.best_300$auc[, temp.order[temp.alt_ind]]    
  print(mean(best2))
  temp.wtest <- wilcox.test(best1, best2, paired = TRUE, alternative = "great")
  print(temp.wtest$p.value)
  stopifnot(temp.wtest$p.value < 0.05)  
}

# AUC
temp.order <- order(apply(cp.best_300$slope, 2, median), decreasing = TRUE)
best1 <- cp.best_300$slope[, temp.order[1]]

for (temp.alt_ind in 2:21) {  
  best2 <- cp.best_300$slope[, temp.order[temp.alt_ind]]    
  print(mean(best2))
  temp.wtest <- wilcox.test(best1, best2, paired = TRUE, alternative = "great")
  print(temp.wtest$p.value)
  stopifnot(temp.wtest$p.value < 0.05)  
}

####################################### CP2P ##############################
cpp.best <- list()
cpp.varying_training_matrix <- matrix(nrow = 0, ncol = 21)
label.types <- c("ic50", "auc", "slope")

for (temp.label_ind in 1:3) {  
  for (temp.num in 1:15) {    
    load(paste0("~/output_temp/Bortezomib/cpp_var/", label.types[temp.label_ind], "/bortezomib_cpp_1to100_", label.types[temp.label_ind], "_var", temp.num, ".RData"))   
    
    # Boxplot Data ---------------------------------
    cpp.boxplot_matrix <- foreach (temp.ind = 1:100, .combine=rbind) %do% {
      c(
        cpp.other_model.l1000[[temp.ind]]$elasticNet$test_auc
        , cpp.other_model.l1000[[temp.ind]]$lasso$test_auc
        , cpp.other_model.l1000[[temp.ind]]$ridge$test_auc
        , cpp.other_model.l1000[[temp.ind]]$rf$test_auc
        , cpp.other_model.l1000[[temp.ind]]$svmLinear$test_auc
        , cpp.other_model.l1000[[temp.ind]]$svmRadial$test_auc
        , cpp.snf.single.l1000[[temp.ind]]$test.auc    
        
        , cpp.other_model.all[[temp.ind]]$elasticNet$test_auc
        , cpp.other_model.all[[temp.ind]]$lasso$test_auc
        , cpp.other_model.all[[temp.ind]]$ridge$test_auc
        , cpp.other_model.all[[temp.ind]]$rf$test_auc
        , cpp.other_model.all[[temp.ind]]$svmLinear$test_auc
        , cpp.other_model.all[[temp.ind]]$svmRadial$test_auc
        , cpp.snf.single.all[[temp.ind]]$test.auc
        
        , cpp.other_model.mRMR1000[[temp.ind]]$elasticNet$test_auc
        , cpp.other_model.mRMR1000[[temp.ind]]$lasso$test_auc
        , cpp.other_model.mRMR1000[[temp.ind]]$ridge$test_auc
        , cpp.other_model.mRMR1000[[temp.ind]]$rf$test_auc
        , cpp.other_model.mRMR1000[[temp.ind]]$svmLinear$test_auc
        , cpp.other_model.mRMR1000[[temp.ind]]$svmRadial$test_auc
        , cpp.snf.single.mRMR1000[[temp.ind]]$test.auc    
      )  
    }
    
    colnames(cpp.boxplot_matrix) <- c(
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
    
    if (temp.label_ind == 1) {
      cpp.best$ic50[[temp.num]] <- cpp.boxplot_matrix   
    } else if (temp.label_ind == 2) {
      cpp.best$auc[[temp.num]] <- cpp.boxplot_matrix     
    } else if (temp.label_ind == 3) {  
      cpp.best$slope[[temp.num]] <- cpp.boxplot_matrix   
    } else {
      stop("wrong temp.label_ind")
    }
  }
}

# IC50
temp.order <- order(apply(cpp.best$ic50[[15]], 2, median), decreasing = TRUE)
best1 <- cpp.best$ic50[[15]][, temp.order[1]]
for (temp.alt_ind in 3:21) {  
  best2 <- cpp.best$ic50[[15]][, temp.order[temp.alt_ind]]    
  temp.wtest <- wilcox.test(best1, best2, paired = TRUE, alternative = "great")
  print(temp.wtest$p.value)
  stopifnot(temp.wtest$p.value < 0.05)  
}

# AUC
temp.order <- order(apply(cpp.best$auc[[15]], 2, median), decreasing = TRUE)
best1 <- cpp.best$auc[[15]][, temp.order[1]]

for (temp.alt_ind in 3:21) {  
  best2 <- cpp.best$auc[[15]][, temp.order[temp.alt_ind]]    
  print(mean(best2))
  temp.wtest <- wilcox.test(best1, best2, paired = TRUE, alternative = "great")
  print(temp.wtest$p.value)
  stopifnot(temp.wtest$p.value < 0.05)  
}

# slope
temp.order <- order(apply(cpp.best$slope[[15]], 2, median), decreasing = TRUE)
best1 <- cpp.best$slope[[15]][, temp.order[1]]

for (temp.alt_ind in 2:21) {  
  best2 <- cpp.best$slope[[15]][, temp.order[temp.alt_ind]]    
  print(mean(best2))
  temp.wtest <- wilcox.test(best1, best2, paired = TRUE, alternative = "great")
  print(temp.wtest$p.value)
  stopifnot(temp.wtest$p.value < 0.05)  
}

### COMPARISON #####
# IC50
temp.order <- order(apply(cpp.best$ic50[[15]], 2, median), decreasing = TRUE)
for (temp.t_ind in 2:6) {
  best.cpp <- cpp.best$ic50[[temp.t_ind]][, temp.order[1]]
  best.pp <- pp.test_auc[[temp.t_ind]][, temp.order[1]]
  temp.wtest <- wilcox.test(best.cpp, best.pp, paired = FALSE, alternative = "great")
  print(temp.wtest$p.value)
}

best.cpp <- vector()
best.pp <- vector()
for (temp.t_ind in 2:6) {
  best.cpp <- c(best.cpp, cpp.best$ic50[[temp.t_ind]][, temp.order[1]])
  best.pp <- c(best.pp, pp.test_auc[[temp.t_ind]][, temp.order[1]])
}

temp.wtest <- wilcox.test(best.cpp, best.pp, paired = FALSE, alternative = "great")
print(temp.wtest$p.value, digits = 3)

# C2P and P2P 
temp.order <- order(apply(cp.best_300$ic50, 2, median), decreasing = TRUE)
best.cp <- cp.best_300$ic50[, temp.order[1]]

temp.order <- order(apply(cpp.best$ic50[[15]], 2, median), decreasing = TRUE)
for (temp.t_ind in 2:3) {  
  best.pp <- pp.test_auc[[temp.t_ind]][, temp.order[1]]  
  temp.wtest <- wilcox.test(best.cp, best.pp, paired = FALSE, alternative = "great")
  print(temp.wtest$p.value)
}

# compare with itself
best.pp <- pp.test_auc[[15]][, temp.order[1]]
for (temp.t_ind in 2:15) {
  best.pp2 <- pp.test_auc[[temp.t_ind]][, temp.order[1]]
  temp.wtest <- wilcox.test(best.pp, best.pp2, paired = FALSE, alternative = "t")
  if (temp.wtest$p.value > 0.05) {
    print(temp.wtest$p.value, digits = 3)  
    print(temp.t_ind)
  }
}

best.cpp <- cpp.best$ic50[[15]][, temp.order[1]]
for (temp.t_ind in 2:15) {
  best.cpp2 <- cpp.best$ic50[[temp.t_ind]][, temp.order[1]]
  temp.wtest <- wilcox.test(best.cpp, best.cpp2, paired = FALSE, alternative = "t")
  if (temp.wtest$p.value > 0.05) {
    print(temp.wtest$p.value, digits = 3)  
    print(temp.t_ind)
  }
}



#### AUC
temp.order <- order(apply(cpp.best$auc[[15]], 2, median), decreasing = TRUE)
for (temp.t_ind in 2:15) {
  best.cpp <- cpp.best$auc[[temp.t_ind]][, temp.order[1]]
  best.pp <- pp.test_auc[[temp.t_ind]][, temp.order[1]]
  temp.wtest <- wilcox.test(best.cpp, best.pp, paired = FALSE, alternative = "great")
  print(temp.wtest$p.value)  
}

# C2P and P2P 
temp.order <- order(apply(cp.best_300$auc, 2, median), decreasing = TRUE)
best.cp <- cp.best_300$auc[, temp.order[1]]

temp.order <- order(apply(pp.test_auc[[15]], 2, median), decreasing = TRUE)
for (temp.t_ind in 2:5) {  
  best.pp <- pp.test_auc[[temp.t_ind]][, temp.order[1]]  
  temp.wtest <- wilcox.test(best.cp, best.pp, paired = FALSE, alternative = "great")
  print(temp.wtest$p.value)
}

# compare with itself
best.cpp <- cpp.best$auc[[15]][, temp.order[1]]
for (temp.t_ind in 2:15) {
  best.cpp2 <- cpp.best$auc[[temp.t_ind]][, temp.order[1]]
  temp.wtest <- wilcox.test(best.cpp, best.cpp2, paired = FALSE, alternative = "t")
  if (temp.wtest$p.value > 0.05) {
    print(temp.wtest$p.value, digits = 3)  
    print(temp.t_ind)
  }
}

# 90 patients similar performance
best.pp <- pp.test_auc[[9]][, temp.order[1]]
for (temp.t_ind in 2:15) {
  best.cpp <- cpp.best$auc[[temp.t_ind]][, temp.order[1]]
  temp.wtest <- wilcox.test(best.pp, best.cpp, paired = FALSE, alternative = "t")
  if (temp.wtest$p.value > 0.05) {
    print(temp.wtest$p.value, digits = 3)  
    print(temp.t_ind)
  }
}

#### Slope
temp.order <- order(apply(cpp.best$slope[[15]], 2, median), decreasing = TRUE)
for (temp.t_ind in 2:15) {
  best.cpp <- cpp.best$slope[[temp.t_ind]][, temp.order[1]]
  best.pp <- pp.test_auc[[temp.t_ind]][, temp.order[1]]
  temp.wtest <- wilcox.test(best.cpp, best.pp, paired = FALSE, alternative = "great")
  print(temp.wtest$p.value)
  
  if (temp.wtest$p.value > 0.05) {
    print(temp.t_ind)
  }
}

# C2P and P2P 
temp.order <- order(apply(cp.best_300$slope, 2, median), decreasing = TRUE)
best.cp <- cp.best_300$slope[, temp.order[1]]

temp.order <- order(apply(pp.test_auc[[15]], 2, median), decreasing = TRUE)
for (temp.t_ind in 2:5) {  
  best.pp <- pp.test_auc[[temp.t_ind]][, temp.order[1]]  
  temp.wtest <- wilcox.test(best.cp, best.pp, paired = FALSE, alternative = "great")
  print(temp.wtest$p.value)
}

# compare with itself
best.cpp <- cpp.best$slope[[15]][, temp.order[1]]
for (temp.t_ind in 2:15) {
  best.cpp2 <- cpp.best$slope[[temp.t_ind]][, temp.order[1]]
  temp.wtest <- wilcox.test(best.cpp, best.cpp2, paired = FALSE, alternative = "t")
  if (temp.wtest$p.value > 0.05) {
    print(temp.wtest$p.value, digits = 3)  
    print(temp.t_ind)
  }
}

# 90 patients similar performance
best.pp <- pp.test_auc[[9]][, temp.order[1]]
for (temp.t_ind in 2:15) {
  best.cpp <- cpp.best$slope[[temp.t_ind]][, temp.order[1]]
  temp.wtest <- wilcox.test(best.pp, best.cpp, paired = FALSE, alternative = "t")
  if (temp.wtest$p.value > 0.05) {
    print(temp.wtest$p.value, digits = 3)  
    print(temp.t_ind)
  }
}

# better than 90 patient performance
best.pp <- pp.test_auc[[9]][, temp.order[1]]
for (temp.t_ind in 1:15) {
  best.cpp <- cpp.best$slope[[temp.t_ind]][, temp.order[1]]
  temp.wtest <- wilcox.test(best.cpp, best.pp, paired = FALSE, alternative = "g")
  #if (temp.wtest$p.value < 0.05) {
    print(temp.wtest$p.value, digits = 3)  
    print(temp.t_ind)
  #}
}