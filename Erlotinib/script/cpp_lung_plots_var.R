library("doParallel")
library("ROCR")
cpp.varying_training_matrix2 <- matrix(nrow = 0, ncol = 21)
label.type <- "lung_slope"

for (temp.num in 13:24) {
  # Boxplot Data ---------------------------------
  if (temp.num != 24) {
    load(paste0("~/output_temp/Erlotinib/cpp_var/lung_only/erlotinib_cpp_1to100_", label.type, "_var", temp.num, ".RData"))       
  } else {
    load(paste0("~/output_temp/Erlotinib/Used/erlotinib_cpp_1to25_", label.type, ".RData"))
  }
  
  # Boxplot Data ---------------------------------
  cpp.predictions <- matrix(nrow = 0, ncol = 21)
  temp.failed_runs <- vector()
  for (temp.ind in 1:length(cpp.other_model.l1000)) {
    temp.row <- c(
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
    if (length(temp.row) %% 21 != 0) {
      print("FAILED RUN!")
      print(temp.ind)
      temp.failed_runs <- c(temp.failed_runs, temp.ind)
      next()
    }
    stopifnot(length(temp.row) %% 21 == 0)       
    cpp.predictions = rbind(cpp.predictions, matrix(temp.row, length(temp.row) / 21, 21))
  }   
  
  if (temp.num != 24) {
    cpp.test_index <- foreach (temp.ind = 1:length(cpp.other_model.all), .combine=c, .errorhandling = 'remove') %do% {
      if (temp.ind %in% temp.failed_runs) {
        print(temp.ind)
        stop()
      }
      input_partition$cpp[[temp.ind]]$test_index
    }  
  } else {
    cpp.test_index = length(input_partition$cell_lines) + 1:length(input_partition$patients)
  }
  stopifnot(dim(cpp.predictions)[1] == length(cpp.test_index))
  cpp.test_auc <- foreach (temp.ind = 1:21, .combine=rbind) %do% {
    if (label.type == "lung_slope") {
      temp.predict <- prediction(cpp.predictions[, temp.ind], erlotinib.labels$lung_all.slope[cpp.test_index])
    } else if (label.type == "lung_auc") {
      temp.predict <- prediction(cpp.predictions[, temp.ind], erlotinib.labels$lung_all.IC50[cpp.test_index])
    } else if (label.type == "lung_ic50") {
      temp.predict <- prediction(cpp.predictions[, temp.ind], erlotinib.labels$lung_all.IC50[cpp.test_index])
    } else {
      stop("Wrong labels")
    }
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
  
  cpp.varying_training_matrix2 <- rbind(cpp.varying_training_matrix2, cpp.test_auc)
  print(dim(cpp.varying_training_matrix2))
}

setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS")
num_of_patients = seq(13, 24, 1)
stopifnot(dim(cpp.varying_training_matrix2)[1] == length(num_of_patients))
save(cpp.varying_training_matrix2, num_of_patients, label.type, file = paste0("cpp_lung_var_", label.type, ".RData"))

###
#load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/cpp_var_ic50.RData")
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/cpp_var_auc.RData")
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/cpp_var_slope.RData")

min(cpp.varying_training_matrix2)
max(cpp.varying_training_matrix2)

stopifnot(dim(cpp.varying_training_matrix2)[2] == 21)
setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/Results/")

png(filename = paste0("erl_cpp_var_", label.type, ".png"), width = 800, height = 800)

if (label.type == "ic50") {
  ylim = c(0.45, 1.1)
} else if (label.type == "auc") {
  ylim = c(0.4, 1.05)
} else if (label.type == "slope") {
  ylim = c(0.35, 1.05)
  #par(yaxp = c(0.4, 0.95, 10))
} else {
  stop("label not defined!")
}

par(mar=c(5.1,7,4.5,3.5))
plot(0, 0, xlim = range(num_of_patients), ylim = ylim, type = "n", 
     ylab = "",
     xlab = "# of Patient Samples in Training Data",
     #las = 1, cex.axis = 1.2, cex.lab = 1.5)
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = c(0.35, 0.9, 11))

title(ylab = "Test AUC", cex.lab = 1.5, line = 4.5)

#axis(side = 4, las = 1, cex.axis = 1.2)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = c(0.35, 0.9, 11))

cl <- rainbow(7)
temp.points <- rep(1, 7)
temp.points <- c(temp.points, rep(2, 7))
temp.points <- c(temp.points, rep(3, 7))
for (temp.ind in 1:21) {
  lines(num_of_patients, cpp.varying_training_matrix2[, temp.ind], col = cl[(temp.ind - 1) %% 7 + 1], type = 'b', pch=temp.points[temp.ind])
}
legend("topleft", legend = colnames(cpp.varying_training_matrix2), col=cl, pch=temp.points, ncol=3, pt.cex = 1, cex = 1.3) # optional legend

dev.off()

