library("doParallel")
library("ROCR")
# Boxplot Data ---------------------------------
cp.predictions <- foreach (temp.ind = 1:length(cp.other_model.l1000), .combine=rbind) %do% {
  c(
    cp.other_model.l1000[[temp.ind]]$elasticNet$predict
    , cp.other_model.l1000[[temp.ind]]$lasso$predict
    , cp.other_model.l1000[[temp.ind]]$ridge$predict
    , cp.other_model.l1000[[temp.ind]]$rf$predict
    , cp.other_model.l1000[[temp.ind]]$svmLinear$predict
    , cp.other_model.l1000[[temp.ind]]$svmRadial$predict
    , cp.snf.single.l1000[[temp.ind]]$test.prediction
    
    , cp.other_model.all[[temp.ind]]$elasticNet$predict
    , cp.other_model.all[[temp.ind]]$lasso$predict
    , cp.other_model.all[[temp.ind]]$ridge$predict
    , cp.other_model.all[[temp.ind]]$rf$predict
    , cp.other_model.all[[temp.ind]]$svmLinear$predict
    , cp.other_model.all[[temp.ind]]$svmRadial$predict
    , cp.snf.single.all[[temp.ind]]$test.prediction
    
    , cp.other_model.mRMR1000[[temp.ind]]$elasticNet$predict
    , cp.other_model.mRMR1000[[temp.ind]]$lasso$predict
    , cp.other_model.mRMR1000[[temp.ind]]$ridge$predict
    , cp.other_model.mRMR1000[[temp.ind]]$rf$predict
    , cp.other_model.mRMR1000[[temp.ind]]$svmLinear$predict
    , cp.other_model.mRMR1000[[temp.ind]]$svmRadial$predict
    , cp.snf.single.mRMR1000[[temp.ind]]$test.prediction    
  )  
}

cp.test_auc <- foreach (temp.ind = 1:dim(cp.predictions)[2], .combine=rbind) %do% {
  temp.predict <- prediction(cp.predictions[, temp.ind], erlotinib.labels$patient)
  unlist(slot(performance(temp.predict, "auc"), "y.values"))  
}

cp.test_auc <- as.vector(cp.test_auc)

names(cp.test_auc) <- c(
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

sort(cp.test_auc)


##################



library("doParallel")

# Boxplot Data ---------------------------------
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

cp.boxplot_matrix[which(cp.boxplot_matrix == max(cp.boxplot_matrix))]

cp.best_lung_only.ic50 <- cp.boxplot_matrix[which(cp.boxplot_matrix == max(cp.boxplot_matrix))]
cp.best_lung_only.slope <- cp.boxplot_matrix[which(cp.boxplot_matrix == max(cp.boxplot_matrix))]
setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS")
save(cp.best_lung_only.ic50, cp.best_lung_only.slope, file="erl_cp_lung_only.RData")
