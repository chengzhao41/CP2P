library("doParallel")
library("RColorBrewer")

# Boxplot Data ---------------------------------
temp.boxplot_matrix <- foreach (temp.ind = 1:100, .combine=rbind) %do% {
  c(
    cpp.snf.single.l1000[[temp.ind]]$test.auc
    , cpp.other_model.l1000[[temp.ind]]$elasticNet$test_auc
    , cpp.other_model.l1000[[temp.ind]]$lasso$test_auc
    , cpp.other_model.l1000[[temp.ind]]$ridge$test_auc
    , cpp.other_model.l1000[[temp.ind]]$rf$test_auc
    , cpp.other_model.l1000[[temp.ind]]$svmLinear$test_auc
    , cpp.other_model.l1000[[temp.ind]]$svmRadial$test_auc    
    
    , cpp.snf.single.all[[temp.ind]]$test.auc
    , cpp.other_model.all[[temp.ind]]$elasticNet$test_auc
    , cpp.other_model.all[[temp.ind]]$lasso$test_auc
    , cpp.other_model.all[[temp.ind]]$ridge$test_auc
    , cpp.other_model.all[[temp.ind]]$rf$test_auc
    , cpp.other_model.all[[temp.ind]]$svmLinear$test_auc
    , cpp.other_model.all[[temp.ind]]$svmRadial$test_auc

    , cpp.snf.single.mRMR1000[[temp.ind]]$test.auc
    , cpp.other_model.mRMR1000[[temp.ind]]$elasticNet$test_auc
    , cpp.other_model.mRMR1000[[temp.ind]]$lasso$test_auc
    , cpp.other_model.mRMR1000[[temp.ind]]$ridge$test_auc
    , cpp.other_model.mRMR1000[[temp.ind]]$rf$test_auc
    , cpp.other_model.mRMR1000[[temp.ind]]$svmLinear$test_auc
    , cpp.other_model.mRMR1000[[temp.ind]]$svmRadial$test_auc    
  )  
}

## error checking
for (temp.ind in 2:100) {
  stopifnot(length(cpp.other_model.mRMR1000[[temp.ind - 1]]$svmRadial$predict) 
  == length(cpp.other_model.mRMR1000[[temp.ind]]$elasticNet$predict))
  stopifnot(length(cpp.snf.single.all[[temp.ind]]$test.prediction) == 
              length(cpp.snf.single.all[[temp.ind - 1]]$test.prediction)
            )
}

temp.boxplot_data <- data.frame(
  temp.boxplot_matrix)

temp.names <- c(
  "SNF L1000"
  , "Elastic Net L1000" 
  , "Lasso L1000" 
  , "Ridge L1000" 
  , "RF L1000" 
  , "SVM lin L1000" 
  , "SVM rbf L1000"             
  
  , "SNF all"
  , "Elastic Net all" 
  , "Lasso all"           
  , "Ridge all"
  , "RF all" 
  , "SVM lin all"
  , "SVM rbf all"     
  
  , "SNF mRMR1000"
  , "Elastic Net mRMR1000"
  , "Lasso mRMR1000"
  , "Ridge mRMR1000"
  , "RF mRMR1000"
  , "SVM lin mRMR1000"
  , "SVM rbf mRMR1000"
)

dev.off()

cl <- brewer.pal(n = 3, name = 'Dark2')
boxplot(temp.boxplot_data, las=2, 
        par(mar = c(10, 4, 2, 2.5) + 0.1), 
        ylab = "Test AUROC",        
        names = temp.names,
        at = rank(apply(temp.boxplot_data, 2, median)), col = c(rep(cl[1], 7), rep(cl[2], 7), rep(cl[3], 7))) 
axis(side = 4, las = 1, cex.lab = 1.5)
par(xpd=TRUE)
legend(x = 0, y = 1.1 , legend=c("all","L1000","mRMR1000"),fill=c(cl[2], cl[1], cl[3]), cex = 0.75, ncol = 3)

temp.order <- order(apply(temp.boxplot_data, 2, median))
temp.median_results <- apply(temp.boxplot_data, 2, median)
names(temp.median_results) <- temp.names
temp.median_results[temp.order]

setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Paper/Table")

temp.label_type <- "ic50"
temp.label_type <- "auc"
temp.label_type <- "slope"
write.csv(temp.median_results[temp.order], file=paste0("bor_cpp_", temp.label_type, ".csv"))

