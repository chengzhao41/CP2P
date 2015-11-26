library("doParallel")
library("ROCR")

#load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/erlotinib_cpp_1to25.RData")

# Boxplot Data ---------------------------------
pp.predictions <- foreach (temp.ind = 1:length(pp.other_model.l1000), .combine=rbind) %do% {
  c(
    pp.other_model.l1000[[temp.ind]]$elasticNet$predict
    , pp.other_model.l1000[[temp.ind]]$lasso$predict
    , pp.other_model.l1000[[temp.ind]]$ridge$predict
    , pp.other_model.l1000[[temp.ind]]$rf$predict
    , pp.other_model.l1000[[temp.ind]]$svmLinear$predict
    , pp.other_model.l1000[[temp.ind]]$svmRadial$predict
    , pp.snf.single.l1000[[temp.ind]]$test.prediction
    
    , pp.other_model.all[[temp.ind]]$elasticNet$predict
    , pp.other_model.all[[temp.ind]]$lasso$predict
    , pp.other_model.all[[temp.ind]]$ridge$predict
    , pp.other_model.all[[temp.ind]]$rf$predict
    , pp.other_model.all[[temp.ind]]$svmLinear$predict
    , pp.other_model.all[[temp.ind]]$svmRadial$predict
    , pp.snf.single.all[[temp.ind]]$test.prediction
    
    , pp.other_model.mRMR1000[[temp.ind]]$elasticNet$predict
    , pp.other_model.mRMR1000[[temp.ind]]$lasso$predict
    , pp.other_model.mRMR1000[[temp.ind]]$ridge$predict
    , pp.other_model.mRMR1000[[temp.ind]]$rf$predict
    , pp.other_model.mRMR1000[[temp.ind]]$svmLinear$predict
    , pp.other_model.mRMR1000[[temp.ind]]$svmRadial$predict
    , pp.snf.single.mRMR1000[[temp.ind]]$test.prediction    
  )  
}

pp.test_auc <- foreach (temp.ind = 1:dim(pp.predictions)[2], .combine=rbind) %do% {
  temp.predict <- prediction(pp.predictions[, temp.ind], erlotinib.labels$patient)
  unlist(slot(performance(temp.predict, "auc"), "y.values"))  
}

pp.test_auc <- as.vector(pp.test_auc)

names(pp.test_auc) <- c(
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

sort(pp.test_auc)
