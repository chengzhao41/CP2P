library("doParallel")
library("ROCR")

# Boxplot Data ---------------------------------
cpp.predictions <- foreach (temp.ind = 1:length(cpp.other_model.l1000), .combine=rbind) %do% {
  c(
    cpp.other_model.l1000[[temp.ind]]$elasticNet$predict
    , cpp.other_model.l1000[[temp.ind]]$lasso$predict
    , cpp.other_model.l1000[[temp.ind]]$ridge$predict
    , cpp.other_model.l1000[[temp.ind]]$rf$predict
    , cpp.other_model.l1000[[temp.ind]]$svmLinear$predict
    , cpp.other_model.l1000[[temp.ind]]$svmRadial$predict
    , cpp.snf.single.l1000[[temp.ind]]$test.prediction
    
    , cpp.other_model.all[[temp.ind]]$elasticNet$predict
    , cpp.other_model.all[[temp.ind]]$lasso$predict
    , cpp.other_model.all[[temp.ind]]$ridge$predict
    , cpp.other_model.all[[temp.ind]]$rf$predict
    , cpp.other_model.all[[temp.ind]]$svmLinear$predict
    , cpp.other_model.all[[temp.ind]]$svmRadial$predict
    , cpp.snf.single.all[[temp.ind]]$test.prediction
    
    , cpp.other_model.mRMR1000[[temp.ind]]$elasticNet$predict
    , cpp.other_model.mRMR1000[[temp.ind]]$lasso$predict
    , cpp.other_model.mRMR1000[[temp.ind]]$ridge$predict
    , cpp.other_model.mRMR1000[[temp.ind]]$rf$predict
    , cpp.other_model.mRMR1000[[temp.ind]]$svmLinear$predict
    , cpp.other_model.mRMR1000[[temp.ind]]$svmRadial$predict
    , cpp.snf.single.mRMR1000[[temp.ind]]$test.prediction    
  )  
}

cpp.test_auc <- foreach (temp.ind = 1:dim(cpp.predictions)[2], .combine=rbind) %do% {
  temp.predict <- prediction(cpp.predictions[, temp.ind], docetaxel.labels$patient)
  unlist(slot(performance(temp.predict, "auc"), "y.values"))  
}

cpp.test_auc <- as.vector(cpp.test_auc)

names(cpp.test_auc) <- c(
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

sort(cpp.test_auc)

setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Paper/Table")

temp.label_type <- "ic50"
temp.label_type <- "auc"
temp.label_type <- "slope"

temp.label_type <- "ic50_breast"
temp.label_type <- "auc_breast"
temp.label_type <- "slope_breast"


write.csv(sort(cpp.test_auc), file=paste0("doc_cpp_", temp.label_type, ".csv"))


