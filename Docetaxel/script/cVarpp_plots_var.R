library("doParallel")
library("RColorBrewer")
library("ROCR")

### USER SUPPLIED PARAMETERS
label.type <- "slope_p20"
### END

cpp.varying_training_matrix <- matrix(nrow = 0, ncol = 21)

for (temp.num in 1:20) {
  # Boxplot Data ---------------------------------
  load(paste0("/Users/chengzhao/Dropbox/output_WS/docetaxel_cpp_100to1_", label.type, "_var", temp.num, ".RData"))  
  
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
      print(temp.ind)
      temp.failed_runs <- c(temp.failed_runs, temp.ind)
      next()
    }
    stopifnot(length(temp.row) %% 21 == 0)       
    cpp.predictions = rbind(cpp.predictions, matrix(temp.row, length(temp.row) / 21, 21))
  }   
  
  cpp.test_index <- foreach (temp.ind = 1:length(input_partition$cpp), .combine=c, .errorhandling = 'remove') %do% {
    if (temp.ind %in% temp.failed_runs) {
      print(temp.ind)
      stop()
    }
    input_partition$cpp[[temp.ind]]$test_index
  }  
  
  stopifnot(dim(cpp.predictions)[1] == length(cpp.test_index))
  cpp.test_auc <- foreach (temp.ind = 1:21, .combine=rbind) %do% {
    if (label.type == "slope" || label.type == "slope_p20") {
      temp.predict <- prediction(cpp.predictions[, temp.ind], docetaxel.labels$slope_combined[cpp.test_index])
    } else if (label.type == "auc" || label.type == "auc_p20") {
      temp.predict <- prediction(cpp.predictions[, temp.ind], docetaxel.labels$AUC_combined[cpp.test_index])
    } else if (label.type == "ic50" || label.type == "ic50_p20") {
      temp.predict <- prediction(cpp.predictions[, temp.ind], docetaxel.labels$IC50_combined[cpp.test_index])
    } else {
      stop("Provide the label.type logic")
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
  
  cpp.varying_training_matrix <- rbind(cpp.varying_training_matrix, cpp.test_auc)
}

num_of_cell_lines = seq(30, 600, 30)
stopifnot(dim(cpp.varying_training_matrix)[1] == length(num_of_cell_lines))
save(cpp.varying_training_matrix, num_of_cell_lines, label.type, file = paste0("Docetaxel/output_WS/cpp_var_", label.type, ".RData"))

###
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/cpp_var_ic50.RData")
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/cpp_var_auc.RData")
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/cpp_var_slope.RData")


stopifnot(dim(cpp.varying_training_matrix)[2] == 21)

png(filename = paste0("Docetaxel/plots/doc_cpp_var_", label.type, ".png"), width = 800, height = 800)

if (label.type == "ic50") {
  ylim = c(0.45, 1.1)
} else if (label.type == "auc") {
  ylim = c(0.5, 1.1)
} else if (label.type == "slope") {
  ylim = c(0.4, 1.05)
  #par(yaxp = c(0.4, 0.95, 10))
} else if (label.type == "ic50_p20") {
  ylim = c(0.65, 1)  
} else if (label.type == "auc_p20") {
  ylim = c(0.7, 1)  
} else if (label.type == "slope_p20") {
  ylim = c(0.7, 0.95)  
} else {
  stop("label not defined!")
}

par(mar=c(5.1,7,4.5,3.5))
plot(0, 0, xlim = range(num_of_cell_lines), ylim = ylim, type = "n", 
     ylab = "",
     xlab = "# of Cell Line Samples in Training Data",
     #las = 1, cex.axis = 1.2, cex.lab = 1.5)
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxp = c(0.4, 1, 12))

title(ylab = "Test AUROC", cex.lab = 1.5, line = 4.5)

#axis(side = 4, las = 1, cex.axis = 1.2)
axis(side = 4, las = 1, cex.axis = 1.2, yaxp = c(0.4, 1, 12))

cl <- rainbow(7)
temp.points <- rep(1, 7)
temp.points <- c(temp.points, rep(2, 7))
temp.points <- c(temp.points, rep(3, 7))
for (temp.ind in 1:21) {
  lines(num_of_cell_lines, cpp.varying_training_matrix[, temp.ind], col = cl[(temp.ind - 1) %% 7 + 1], type = 'b', pch=temp.points[temp.ind])
}
legend("topleft", legend = colnames(cpp.varying_training_matrix), col=cl, pch=temp.points, ncol=3, pt.cex = 1, cex = 1.3) # optional legend

dev.off()

