library("doParallel")

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
temp.boxplot_data <- data.frame(
  cp.boxplot_matrix)

temp.names <- c(
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
temp.order <- order(apply(temp.boxplot_data, 2, median))
#with(temp.boxplot_data, reorder(TYPE, temp.boxplot_data, median))

dev.off()
boxplot(temp.boxplot_data, las=2, 
        par(mar = c(10, 5, 4, 2) + 0.1), 
        ylab = "AUC on test set",
        names = temp.names,
        at = rank(apply(temp.boxplot_data, 2, median), ties.method = c("first")))
axis(side = 4)

temp.median_results <- apply(temp.boxplot_data, 2, median)[temp.order]
names(temp.median_results) <- temp.names[temp.order]
temp.median_results

setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Paper/Table")

temp.label_type <- "ic50"
temp.label_type <- "auc"
temp.label_type <- "slope"

temp.label_type <- "ic50_breast"
temp.label_type <- "auc_breast"
temp.label_type <- "slope_breast"

write.csv(temp.median_results, file=paste0("doc_cp_", temp.label_type, ".csv"))
