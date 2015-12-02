sva_combine <- function(batch, label, input_data, method = "fast", num_sv_method = "leek", n.sv = 0) {
  require("sva")
  
  stopifnot(length(label) == length(batch))
  stopifnot(dim(input_data)[1] == length(label))
  
  temp.batch <- as.factor(batch)
  levels(temp.batch) <- 1:length(levels(temp.batch))
  
  temp.pheno = data.frame(batch = temp.batch, label=label)
  temp.mod = model.matrix(~as.factor(label), data=temp.pheno)
  temp.mod0 = model.matrix(~1, data=temp.pheno) 
  
  if (n.sv == 0) {
    if (num_sv_method == "leek") {
      temp.n.sv = num.sv(dat = t(input_data), mod = temp.mod, method="leek")
      print(paste("leek n.sv:", temp.n.sv))
    }    
    if (num_sv_method == "be") {
      stopifnot(FALSE)
      temp.n.sv = num.sv(t(input_data), temp.mod, method="be")
      print(paste("be n.sv:", temp.n.sv))
    }
  } else {
    temp.n.sv = n.sv
  }
  
  temp.svobj = sva(t(input_data), mod=temp.mod, mod0=temp.mod0, n.sv=temp.n.sv)
  temp.sva <- fsva(dbdat=t(input_data), mod=temp.mod, sv=temp.svobj, newdat=t(input_data), method = method)
  return(t(temp.sva$db))  
}

fsva_combine <- function(batch, label, input_data, training_ind, test_ind, method = "fast", num_sv_method = "leek", n.sv = 0) {
  require("sva")  
  
  stopifnot(length(label) == length(batch))
  stopifnot(dim(input_data)[1] == length(label))  
  stopifnot(length(training_ind) + length(test_ind) == dim(input_data)[1])
  
  temp.batch <- as.factor(batch)
  levels(temp.batch) <- 1:length(levels(temp.batch))  
  
  input_data.training <- input_data[training_ind, ]
  input_data.test <- input_data[test_ind, ]
  batch.training <- batch[training_ind]
  
  label.training <- label[training_ind]
  
  temp.pheno = data.frame(batch = batch.training, label=label.training)
  temp.mod = model.matrix(~as.factor(label.training), data=temp.pheno)
  temp.mod0 = model.matrix(~1, data=temp.pheno)  
  
  if (n.sv == 0) {
    if (num_sv_method == "leek") {
      temp.n.sv = num.sv(dat = t(input_data.training), mod = temp.mod, method="leek")
      print(paste("leek n.sv:", temp.n.sv))
    }  
    if (num_sv_method == "be") {
      temp.n.sv = num.sv(t(input_data.training), temp.mod, method="be")
      print(paste("be n.sv:", temp.n.sv))
    }
  } else {
    temp.n.sv = n.sv
  }
  
  temp.svobj = sva(t(input_data.training), mod=temp.mod, mod0=temp.mod0, n.sv=temp.n.sv)
  temp.sva <- fsva(dbdat=t(input_data.training), mod=temp.mod, sv=temp.svobj, newdat=t(input_data.test), method = method)
  
  temp.output <- input_data
  temp.output[training_ind, ] <- t(temp.sva$db)
  temp.output[test_ind, ] <- t(temp.sva$new)
  
  return(temp.output)
}


ComBat_combine <- function(input_data, label, batch, par.prior = TRUE, prior.plots = TRUE) {    
  require("sva")
  levels(batch) <- 1:length(levels(batch))
  cell_line.pheno = data.frame(batch = batch, label=label)
  cell_line.mod = model.matrix(~as.factor(label), data=cell_line.pheno)
  
  ## Use ComBat  
  temp.ComBat = ComBat(dat=t(input_data), batch = batch, mod = cell_line.mod, par.prior = par.prior, prior.plots=prior.plots)
  return(t(temp.ComBat))  
}


show_pca <- function(input_data, label, choices = c(1,2), pca_line_plot = TRUE, print = TRUE) {
  require(ggbiplot)  
  stopifnot(dim(input_data)[1] == length(label))
  pca <- prcomp(input_data, center = TRUE, scale. = TRUE)
  if (pca_line_plot) {
    plot(pca, type = "l", main="Variances of PCA")
  }
  
  # by response  
  g <- ggbiplot(pca, obs.scale = 1, choices = choices, var.scale = 1, 
                groups = as.factor(label), ellipse = TRUE, 
                circle = TRUE, varname.size=0, var.axes = FALSE, 
                labels.size = 10)
  if (print == TRUE) {
    print(g + theme_bw())  
  }
  return(g + theme_bw())
}

remove_features <- function(input_data, input_ind, input_label) {
  
  ## exclude any features that have all NAs
  temp.exclude_feature_ind <- c()
  for (temp.f in 1:dim(input_data)[2]) {
    if (sum(!is.na(input_data[, temp.f])) == 0) {
      temp.exclude_feature_ind <- c(temp.exclude_feature_ind, temp.f)
    }  
  }
  print(paste("# of features with all NAs: ", length(temp.exclude_feature_ind)))
  stopifnot(length(temp.exclude_feature_ind) == 0)
  
  ## exclude samples with any NAs
  temp.exclude <- c()
  for (temp.samp in 1:dim(input_data)[1]) {
    if (sum(is.na(input_data[temp.samp, ])) > 0) {
      temp.exclude <- c(temp.exclude, temp.samp)
    }  
  }
  print(paste("# of samples with all NAs: ", length(temp.exclude)))    
  
  input_data <- input_data[-temp.exclude, ]
  stopifnot(sum(is.na(input_data)) == 0)
  input_label <- input_label[-temp.exclude]
  input_ind <- input_ind[-temp.exclude]
  stopifnot(dim(input_data)[1] == length(input_label))
  stopifnot(length(input_ind) == length(input_label))
  
  return(
    list(data = input_data,
         ind = input_ind,
         label = input_label))
}


