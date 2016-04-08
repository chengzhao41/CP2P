# generates the random partitions with the test set common to all approaches
generate_random_partition <- function(
  labels.cell_lines, 
  labels.patient,
  num.training.p,
  cell_line_order,
  num.training.c, 
  metric = 'auc',
  input_partition = NULL,
  num.min_labels.training.p = 8,
  num.min_labels.training.c = 8,
  num.min_labels.test = 5,
  num.partitions = 100,
  num.test_size = NULL) {

  stopifnot(length(labels.cell_lines) > 0)
  stopifnot(length(labels.cell_lines) == length(cell_line_order))
  stopifnot(num.training.p >= 10)
  stopifnot(metric %in% c('auc', 'acc'))
  stopifnot(length(labels.patient) > num.training.p)
  stopifnot(num.training.p >= (num.min_labels.training.p * 2))
  
  require("doParallel")

  ind.patients <- 1:length(labels.patient)
  ind.cell_lines <- list()
  
  if (!is.null(labels.cell_lines$slope)) {
    stopifnot(length(labels.cell_lines$slope) >= num.training.c)
    stopifnot(min(table(labels.cell_lines$slope)) >= num.min_labels.training.c)
    
    ind.cell_lines$slope <- length(ind.patients) + 1:length(labels.cell_lines$slope)
    
    ind.cell_lines$slope <- getCellLineBalancedOrder(
      ind = ind.cell_lines$slope,
      order = cell_line_order$slope,
      labels = labels.cell_lines$slope,
      num.min_labels = num.min_labels.training.c)
  } else {
    warning("slope response labels are not provided")
  }
  if (!is.null(labels.cell_lines$AUC)) {
    stopifnot(length(labels.cell_lines$AUC) >= num.training.c)
    stopifnot(min(table(labels.cell_lines$AUC)) >= num.min_labels.training.c)
    ind.cell_lines$AUC <- length(ind.patients) + 1:length(labels.cell_lines$AUC)
    
    ind.cell_lines$AUC <- getCellLineBalancedOrder(
      ind = ind.cell_lines$AUC,
      order = cell_line_order$AUC,
      labels = labels.cell_lines$AUC,
      num.min_labels = num.min_labels.training.c)
  } else {
    warning("AUC response labels are not provided")
  }
  if (!is.null(labels.cell_lines$IC50)) {
    stopifnot(length(labels.cell_lines$IC50) >= num.training.c)
    stopifnot(min(table(labels.cell_lines$IC50)) >= num.min_labels.training.c)
    ind.cell_lines$IC50 <- length(ind.patients) + 1:length(labels.cell_lines$IC50)
    
    ind.cell_lines$IC50 <- getCellLineBalancedOrder(
      ind = ind.cell_lines$IC50,
      order = cell_line_order$IC50,
      labels = labels.cell_lines$IC50,
      num.min_labels = num.min_labels.training.c)
  } else {
    warning("IC50 response labels are not provided")
  }
  
  # get the true and false class labels indices
  ind.patients.true <- ind.patients[which(labels.patient == TRUE)]
  ind.patient.false <- ind.patients[which(labels.patient == FALSE)]

  partition <- list()
  
  if (length(labels.patient) != (num.training.p + 1)) { # it is not leave-one-out
  
    for (ind.partition in 1:num.partitions) {
      
      training_index.cp2p <- list()
      training_index.c2p <- list()
      training_index.p2p <- vector()

      if (is.null(input_partition)) {
        
        training_index.p2p <- getMinBalancedInd(
          num.min_labels = num.min_labels.training.p,
          ind.true = ind.patients.true,
          ind.false = ind.patient.false)
        
        test_index <- getMinBalancedInd(
          num.min_labels = num.min_labels.test,
          ind.true = ind.patients.true[-which(ind.patients.true %in% training_index.p2p)],
          ind.false = ind.patient.false[-which(ind.patient.false %in% training_index.p2p)])
        
        training_index.p2p <- c(training_index.p2p, 
                                sample(ind.patients[-which(ind.patients %in% c(training_index.p2p, test_index))], 
                                       num.training.p - (num.min_labels.training.p * 2)))
        
        stopifnot(min(table(labels.patient[training_index.p2p])) >= num.min_labels.training.p)
        stopifnot(length(training_index.p2p) == num.training.p)
        
        test_index <- setdiff(ind.patients, training_index.p2p)
        
        if (!is.null(num.test_size)) {
          stopifnot(length(test_index) >= num.test_size)
          test_index.true <- test_index[which(labels.patient[test_index] == TRUE)]
          test_index.false <- test_index[which(labels.patient[test_index] == FALSE)]
          
          test_index.min_balanced <- getMinBalancedInd(
            num.min_labels = num.min_labels.test,
            ind.true = test_index.true,
            ind.false = test_index.false)
          
          if ((num.test_size - 2 * num.min_labels.test) > 0) {
            test_index <- sample(test_index[-which(test_index %in% test_index.min_balanced)], num.test_size - 2 * num.min_labels.test) 
            test_index <- c(test_index, test_index.min_balanced)
          } else {
            test_index <- test_index.min_balanced
          }
        }
        
        stopifnot(min(table(labels.patient[test_index])) >= num.min_labels.test)
        
      } else {
        training_index.p2p <- input_partition$p2p[[ind.partition]]$training_index
        training_index.p2p <- training_index.p2p[1:num.training.p]
        stopifnot(length(training_index.p2p) == num.training.p)
        test_index <- input_partition$p2p[[ind.partition]]$test_index
      }
      
      if (length(ind.cell_lines$slope) > 0) {
        training_index.c2p$slope = ind.cell_lines$slope[1:num.training.c]
        training_index.cp2p$slope <- c(training_index.c2p$slope, training_index.p2p)
        stopifnot(length(training_index.c2p$slope) == num.training.c)
        stopifnot(min(table(labels.cell_lines$slope[training_index.c2p$slope - length(ind.patients)])) >= num.min_labels.training.c)
      }
      if (length(ind.cell_lines$AUC) > 0) {
        training_index.c2p$AUC = ind.cell_lines$AUC[1:num.training.c]
        training_index.cp2p$AUC <- c(training_index.c2p$AUC, training_index.p2p)
        stopifnot(length(training_index.c2p$AUC) == num.training.c)
        stopifnot(min(table(labels.cell_lines$AUC[training_index.c2p$AUC - length(ind.patients)])) >= num.min_labels.training.c)
      }
      if (length(ind.cell_lines$IC50) > 0) {
        training_index.c2p$IC50 = ind.cell_lines$IC50[1:num.training.c]
        training_index.cp2p$IC50 <- c(training_index.c2p$IC50, training_index.p2p)
        stopifnot(length(training_index.c2p$IC50) == num.training.c)
        stopifnot(min(table(labels.cell_lines$IC50[training_index.c2p$IC50 - length(ind.patients)])) >= num.min_labels.training.c)
      }
      
      stopifnot(length(table(labels.patient[test_index])) == 2)
      stopifnot(min(table(labels.patient[test_index])) >= num.min_labels.test)
      stopifnot(length(table(labels.patient[training_index.p2p])) == 2)
      stopifnot(min(table(labels.patient[training_index.p2p])) >= num.min_labels.training.p)
      
      stopifnot(length(training_index.p2p) == num.training.p)
      stopifnot(length(intersect(test_index, training_index.p2p)) == 0)
      stopifnot(test_index %in% ind.patients)
      
      if (length(ind.cell_lines$slope) > 0) {
        stopifnot(length(training_index.cp2p$slope) == num.training.p + num.training.c)
        stopifnot(length(intersect(test_index, training_index.cp2p$slope)) == 0)
        stopifnot(length(intersect(test_index, training_index.c2p$slope)) == 0)
      }
      if (length(ind.cell_lines$AUC) > 0) {
        stopifnot(length(training_index.cp2p$AUC) == num.training.p + num.training.c)      
        stopifnot(length(intersect(test_index, training_index.cp2p$AUC)) == 0)      
        stopifnot(length(intersect(test_index, training_index.c2p$AUC)) == 0)
      }
      if (length(ind.cell_lines$IC50) > 0) {
        stopifnot(length(training_index.cp2p$IC50) == num.training.p + num.training.c)
        stopifnot(length(intersect(test_index, training_index.cp2p$IC50)) == 0)
        stopifnot(length(intersect(test_index, training_index.c2p$IC50)) == 0)
      }
      
      partition$p2p[[ind.partition]] <- list(test_index = test_index
                                           , training_index = training_index.p2p)
      
      if (length(ind.cell_lines$slope) > 0) {
        partition$c2p.slope[[ind.partition]] <- list(test_index = test_index
                                             , training_index = training_index.c2p$slope)
        partition$cp2p.slope[[ind.partition]] <- list(test_index = test_index
                                             , training_index = training_index.cp2p$slope)
      }
      if (length(ind.cell_lines$AUC) > 0) {
        partition$c2p.AUC[[ind.partition]] <- list(test_index = test_index
                                                  , training_index = training_index.c2p$AUC)
        partition$cp2p.AUC[[ind.partition]] <- list(test_index = test_index
                                                     , training_index = training_index.cp2p$AUC)
      }
      if (length(ind.cell_lines$IC50) > 0) {
        partition$c2p.IC50[[ind.partition]] <- list(test_index = test_index
                                                   , training_index = training_index.c2p$IC50)
        partition$cp2p.IC50[[ind.partition]] <- list(test_index = test_index
                                                     , training_index = training_index.cp2p$IC50)
      }
    }
    
    if (num.partitions > 1) {
      for (ind.partition in 1:(num.partitions - 1)) {
        stopifnot((partition$p2p[[ind.partition]]$test_index == partition$p2p[[ind.partition + 1]]$test_index) != length(partition$p2p[[ind.partition + 1]]$test_index))
        stopifnot((partition$p2p[[ind.partition]]$training_index == partition$p2p[[ind.partition + 1]]$training_index) != length(partition$p2p[[ind.partition + 1]]$training_index))
        
        if (length(ind.cell_lines$slope) > 0) {
          stopifnot((partition$c2p.slope[[ind.partition]]$training_index == partition$c2p.slope[[ind.partition + 1]]$training_index) != length(partition$c2p.slope[[ind.partition + 1]]$training_index))
          stopifnot((partition$cp2p.slope[[ind.partition]]$training_index == partition$cp2p.slope[[ind.partition + 1]]$training_index) != length(partition$cp2p.slope[[ind.partition + 1]]$training_index))
        }
        if (length(ind.cell_lines$AUC) > 0) {
          stopifnot((partition$c2p.AUC[[ind.partition]]$training_index == partition$c2p.AUC[[ind.partition + 1]]$training_index) != length(partition$c2p.AUC[[ind.partition + 1]]$training_index))
          stopifnot((partition$cp2p.AUC[[ind.partition]]$training_index == partition$cp2p.AUC[[ind.partition + 1]]$training_index) != length(partition$cp2p.AUC[[ind.partition + 1]]$training_index))
        }
        if (length(ind.cell_lines$IC50) > 0) {
          stopifnot((partition$c2p.IC50[[ind.partition]]$training_index == partition$c2p.IC50[[ind.partition + 1]]$training_index) != length(partition$c2p.IC50[[ind.partition + 1]]$training_index))
          stopifnot((partition$cp2p.IC50[[ind.partition]]$training_index == partition$cp2p.IC50[[ind.partition + 1]]$training_index) != length(partition$cp2p.IC50[[ind.partition + 1]]$training_index))
        }
      }
    }
  
  } else {
    warning("Warning doing leave-one-out partitioning")
    for (ind.partition in 1:length(labels.patient)) {
      
      # if the training and test set does not contain at least 5 samples of each label, then resample
      temp.loop_count = 0
      training_index.cp2p <- list()
      training_index.c2p <- list()
      training_index.p2p <- vector()
      
      test_index <- ind.patients[ind.partition]
      training_index.p2p <- setdiff(ind.patients, test_index)
      stopifnot(length(training_index.p2p) == num.training.p)
      
      if (length(ind.cell_lines$slope) > 0) {
        training_index.c2p$slope = ind.cell_lines$slope[1:num.training.c]
        training_index.cp2p$slope <- c(training_index.c2p$slope, training_index.p2p)
        stopifnot(length(training_index.c2p$slope) == num.training.c)
      }
      if (length(ind.cell_lines$AUC) > 0) {
        training_index.c2p$AUC = ind.cell_lines$AUC[1:num.training.c]
        training_index.cp2p$AUC <- c(training_index.c2p$AUC, training_index.p2p)
        stopifnot(length(training_index.c2p$AUC) == num.training.c)
      }
      if (length(ind.cell_lines$IC50) > 0) {
        training_index.c2p$IC50 = ind.cell_lines$IC50[1:num.training.c]
        training_index.cp2p$IC50 <- c(training_index.c2p$IC50, training_index.p2p)
        stopifnot(length(training_index.c2p$IC50) == num.training.c)
      }

      stopifnot(length(training_index.p2p) == num.training.p)
      stopifnot(length(intersect(test_index, training_index.p2p)) == 0)
      stopifnot(test_index %in% ind.patients)
      
      if (length(ind.cell_lines$slope) > 0) {
        stopifnot(length(training_index.cp2p$slope) == num.training.p + num.training.c)
        stopifnot(length(intersect(test_index, training_index.cp2p$slope)) == 0)
        stopifnot(length(intersect(test_index, training_index.c2p$slope)) == 0)
      }
      if (length(ind.cell_lines$AUC) > 0) {
        stopifnot(length(training_index.cp2p$AUC) == num.training.p + num.training.c)      
        stopifnot(length(intersect(test_index, training_index.cp2p$AUC)) == 0)      
        stopifnot(length(intersect(test_index, training_index.c2p$AUC)) == 0)
      }
      if (length(ind.cell_lines$IC50) > 0) {
        stopifnot(length(training_index.cp2p$IC50) == num.training.p + num.training.c)
        stopifnot(length(intersect(test_index, training_index.cp2p$IC50)) == 0)
        stopifnot(length(intersect(test_index, training_index.c2p$IC50)) == 0)
      }
      
      partition$p2p[[ind.partition]] <- list(test_index = test_index
                                            , training_index = training_index.p2p)
      
      if (length(ind.cell_lines$slope) > 0) {
        partition$c2p.slope[[ind.partition]] <- list(test_index = test_index
                                                    , training_index = training_index.c2p$slope)
        partition$cp2p.slope[[ind.partition]] <- list(test_index = test_index
                                                     , training_index = training_index.cp2p$slope)
      }
      if (length(ind.cell_lines$AUC) > 0) {
        partition$c2p.AUC[[ind.partition]] <- list(test_index = test_index
                                                  , training_index = training_index.c2p$AUC)
        partition$cp2p.AUC[[ind.partition]] <- list(test_index = test_index
                                                   , training_index = training_index.cp2p$AUC)
      }
      if (length(ind.cell_lines$IC50) > 0) {
        partition$c2p.IC50[[ind.partition]] <- list(test_index = test_index
                                                   , training_index = training_index.c2p$IC50)
        partition$cp2p.IC50[[ind.partition]] <- list(test_index = test_index
                                                    , training_index = training_index.cp2p$IC50)
      }
    }
    
    for (ind.partition in 1:(length(labels.patient) - 1)) {
      stopifnot((partition$p2p[[ind.partition]]$test_index == partition$p2p[[ind.partition + 1]]$test_index) != length(partition$p2p[[ind.partition + 1]]$test_index))
      stopifnot((partition$p2p[[ind.partition]]$training_index == partition$p2p[[ind.partition + 1]]$training_index) != length(partition$p2p[[ind.partition + 1]]$training_index))
      
      if (length(ind.cell_lines$slope) > 0) {
        stopifnot((partition$c2p.slope[[ind.partition]]$training_index == partition$c2p.slope[[ind.partition + 1]]$training_index) != length(partition$c2p.slope[[ind.partition + 1]]$training_index))
        stopifnot((partition$cp2p.slope[[ind.partition]]$training_index == partition$cp2p.slope[[ind.partition + 1]]$training_index) != length(partition$cp2p.slope[[ind.partition + 1]]$training_index))
      }
      if (length(ind.cell_lines$AUC) > 0) {
        stopifnot((partition$c2p.AUC[[ind.partition]]$training_index == partition$c2p.AUC[[ind.partition + 1]]$training_index) != length(partition$c2p.AUC[[ind.partition + 1]]$training_index))
        stopifnot((partition$cp2p.AUC[[ind.partition]]$training_index == partition$cp2p.AUC[[ind.partition + 1]]$training_index) != length(partition$cp2p.AUC[[ind.partition + 1]]$training_index))
      }
      if (length(ind.cell_lines$IC50) > 0) {
        stopifnot((partition$c2p.IC50[[ind.partition]]$training_index == partition$c2p.IC50[[ind.partition + 1]]$training_index) != length(partition$c2p.IC50[[ind.partition + 1]]$training_index))
        stopifnot((partition$cp2p.IC50[[ind.partition]]$training_index == partition$cp2p.IC50[[ind.partition + 1]]$training_index) != length(partition$cp2p.IC50[[ind.partition + 1]]$training_index))
      }
    }
  }
      
  return(partition)
}

# generates the random partitions with the test set common to all approaches
getMinBalancedInd <- function(
  num.min_labels, 
  ind.true,
  ind.false) {
  stopifnot(num.min_labels >= 1)
  stopifnot(length(intersect(ind.true, ind.false)) == 0)
  stopifnot(length(ind.true) >= num.min_labels)
  stopifnot(length(ind.false) >= num.min_labels)
  
  ind <- sample(ind.true, num.min_labels)
  ind <- c(rbind(ind, sample(ind.false, num.min_labels)))
  
  return (ind)
}

getCellLineBalancedOrder <- function(
  ind,
  order,
  labels,
  num.min_labels) {
  
  stopifnot(length(order) == length(ind))
  stopifnot(length(order) == length(labels))
  stopifnot(min(table(labels)) >= num.min_labels)
  
  ind_in_order <- ind[order]
  labels_in_order <- labels[order]
  
  ind_true_in_order <- ind_in_order[which(labels_in_order == TRUE)]
  ind_false_in_order <- ind_in_order[which(labels_in_order == FALSE)]
  
  ind_balanced <- getMinBalancedInd(
    num.min_labels = num.min_labels,
    ind.true = ind_true_in_order,
    ind.false = ind_false_in_order)
  
  check_size = length(ind_in_order)
  ind_in_order <- ind_in_order[-which(ind_in_order %in% ind_balanced)]
  ind_balanced_and_in_order <- c(ind_balanced, ind_in_order)
  stopifnot(check_size == length(unique(ind_balanced_and_in_order)))
  
  return(ind_balanced_and_in_order)
}
