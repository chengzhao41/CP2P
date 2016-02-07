library("doParallel")

pp.varying_training_matrix <- matrix(nrow = 0, ncol = 21)

for (temp.num in 12:23) {
  
  if (temp.num != 23) {
    load(paste0("~/output_temp/Erlotinib/pp_var/erlotinib_pp_1to100_var", temp.num, ".RData"))   
  } else {
    load("~/output_temp/Erlotinib/Used/erlotinib_pp_1to25.RData")
  }
  
  # Boxplot Data ---------------------------------
  pp.predictions <- matrix(nrow = 0, ncol = 21)
  temp.failed_runs <- vector()
  for (temp.ind in 1:length(pp.other_model.l1000)) {
    temp.row <- c(
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
    if (length(temp.row) %% 21 != 0) {
      print(temp.ind)
      temp.failed_runs <- c(temp.failed_runs, temp.ind)
      next()
    }
    stopifnot(length(temp.row) %% 21 == 0)       
    pp.predictions = rbind(pp.predictions, matrix(temp.row, length(temp.row) / 21, 21))
  }   
  
  if (temp.num != 23) {
    pp.test_index <- foreach (temp.ind = 1:length(input_partition$pp), .combine=c, .errorhandling = 'remove') %do% {
      if (temp.ind %in% temp.failed_runs) {
        print(temp.ind)
        stop()
      }
      input_partition$pp[[temp.ind]]$test_index
    }  
  } else {
    pp.test_index = 1:25
  }
  stopifnot(dim(pp.predictions)[1] == length(pp.test_index))
  pp.test_auc <- foreach (temp.ind = 1:21, .combine=rbind) %do% {
    temp.predict <- prediction(pp.predictions[, temp.ind], erlotinib.labels$patient[pp.test_index])
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
  
  pp.varying_training_matrix <- rbind(pp.varying_training_matrix, pp.test_auc)
}

setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS")
num_of_patients = seq(13, 24, 1)
save(pp.varying_training_matrix, num_of_patients, file = "pp_var.RData")

#load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/pp_var.RData")
stopifnot(dim(pp.varying_training_matrix)[2] == 21)
dev.off()

png(filename = paste0("erl_pp_var.png"), width = 800, height = 800)

setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/Results")

par(mar=c(5.1,7,4.5,3.5))
plot(0, 0, xlim = range(num_of_patients), ylim = c(0.5, 1.1), type = "n", 
     ylab = "",
     xlab = "# of Patient Samples in Training Data",
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = c(0.5, 1, 10))

title(ylab = "Test AUROC", cex.lab = 1.5, line =4.5)

axis(side = 4, las = 1, cex.axis = 1.2, yaxp = c(0.5, 1, 10))

cl <- rainbow(7)
temp.points <- rep(1, 7)
temp.points <- c(temp.points, rep(2, 7))
temp.points <- c(temp.points, rep(3, 7))
for (temp.ind in 1:21) {  
  lines(num_of_patients, pp.varying_training_matrix[, temp.ind], col = cl[(temp.ind - 1) %% 7 + 1], type = 'b', pch=temp.points[temp.ind])    
}
legend("topleft", legend = colnames(pp.varying_training_matrix), col=cl, pch=temp.points, ncol=3, pt.cex = 1, cex = 1.3)

dev.off()
