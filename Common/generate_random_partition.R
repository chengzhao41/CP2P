# generates the random partitions with the test set common to all approaches
generate_random_partition <- function(
  labels_cell_lines, 
  labels_patient,
  training_amount.p,
  leave_one_out = NULL,
  cell_line_order,
  training_amount.c, 
  acc_training = NULL,
  input_partition = NULL,
  force_balance = TRUE # if the test set is given
  ) {
  
  stopifnot(!is.null(acc_training))
  stopifnot(length(labels_cell_lines) > 0)
  stopifnot(length(labels_cell_lines) == length(cell_line_order))
  stopifnot(training_amount.p >= 10)
  stopifnot(length(labels_patient) > training_amount.p)
  stopifnot(min(table(labels_patient)) >= 8)
  
  require("doParallel")

  temp.patients <- 1:length(labels_patient)
  temp.cell_lines <- list()
  
  if (!is.null(labels_cell_lines$slope)) {
    stopifnot(length(labels_cell_lines$slope) >= training_amount.c)
    temp.cell_lines$slope <- length(temp.patients) + 1:length(labels_cell_lines$slope)
  } else {
    warning("slope response labels are not provided")
  }
  if (!is.null(labels_cell_lines$AUC)) {
    stopifnot(length(labels_cell_lines$AUC) >= training_amount.c)
    temp.cell_lines$AUC <- length(temp.patients) + 1:length(labels_cell_lines$AUC)
  } else {
    warning("AUC response labels are not provided")
  }
  if (!is.null(labels_cell_lines$IC50)) {
    stopifnot(length(labels_cell_lines$IC50) >= training_amount.c)
    temp.cell_lines$IC50 <- length(temp.patients) + 1:length(labels_cell_lines$IC50)
  } else {
    warning("IC50 response labels are not provided")
  }
  
  # get the true and false class labels indices
  temp.patient.true <- which(labels_patient == TRUE)
  temp.patient.false <- which(labels_patient == FALSE)
  stopifnot(length(temp.patient.true) >= 8)
  stopifnot(length(temp.patient.false) >= 8)
  

  partition <- list()
  
  if (length(labels_patient) != (training_amount.p + 1)) {
  
    for (temp.run_ind in 1:100) {
      
      # if the training and test set does not contain at least 5 samples of each label, then resample
      temp.training_index.cp2p <- list()
      temp.training_index.c2p <- list()
      temp.training_index.p2p <- vector()
      
      repeat {
        
        if (is.null(input_partition)) {
          if (force_balance) {
            temp.training_index.p2p <- sample(temp.patient.false, 8)
            temp.training_index.p2p <- c(temp.training_index.p2p, sample(temp.patient.true, 8))
            
            temp.test_index <- sample(temp.patient.false[-which(temp.patient.false %in% temp.training_index.p2p)], 8)
            temp.test_index <- c(temp.test_index, sample(temp.patient.true[-which(temp.patient.true %in% temp.training_index.p2p)], 8))
            
            temp.training_index.p2p <- c(temp.training_index.p2p, sample(temp.patients[-c(temp.training_index.p2p, temp.test_index)], training_amount.p - 16))
            stopifnot(min(table(labels_patient[temp.training_index.p2p])) >= 8)
          } else {
            temp.training_index.p2p <- sample(temp.patient.false, 5)
            temp.training_index.p2p <- c(temp.training_index.p2p, sample(temp.patient.true, 5))
            temp.training_index.p2p <- c(temp.training_index.p2p, sample(temp.patients[-c(temp.training_index.p2p)], training_amount.p - 10))
            stopifnot(min(table(labels_patient[temp.training_index.p2p])) >= 5)
          }
          stopifnot(length(temp.training_index.p2p) == training_amount.p)
          temp.test_index <- setdiff(temp.patients, temp.training_index.p2p)
          if (force_balance) {
            stopifnot(min(table(labels_patient[temp.test_index])) >= 5)
          }
        } else {
          temp.training_index.p2p <- input_partition$p2p[[temp.run_ind]]$training_index
          stopifnot(length(temp.training_index.p2p) == training_amount.p)
          temp.test_index <- input_partition$p2p[[temp.run_ind]]$test_index
        }
        
        if (length(temp.cell_lines$slope) > 0) {
          temp.training_index.c2p$slope = temp.cell_lines$slope[cell_line_order$slope][1:training_amount.c]
          temp.training_index.cp2p$slope <- c(temp.training_index.c2p$slope, temp.training_index.p2p)
          stopifnot(length(temp.training_index.c2p$slope) == training_amount.c)
        }
        if (length(temp.cell_lines$AUC) > 0) {
          temp.training_index.c2p$AUC = temp.cell_lines$AUC[cell_line_order$AUC][1:training_amount.c]
          temp.training_index.cp2p$AUC <- c(temp.training_index.c2p$AUC, temp.training_index.p2p)
          stopifnot(length(temp.training_index.c2p$AUC) == training_amount.c)
        }
        if (length(temp.cell_lines$IC50) > 0) {
          temp.training_index.c2p$IC50 = temp.cell_lines$IC50[cell_line_order$IC50][1:training_amount.c]
          temp.training_index.cp2p$IC50 <- c(temp.training_index.c2p$IC50, temp.training_index.p2p)
          stopifnot(length(temp.training_index.c2p$IC50) == training_amount.c)
        }
        
        if (acc_training) {
          if ((length(table(labels_patient[temp.training_index.p2p])) == 2 && min(table(labels_patient[temp.training_index.p2p])) >= 5)) {
            break
          }          
        } else {
          if ((length(table(labels_patient[temp.test_index])) == 2 && min(table(labels_patient[temp.test_index])) >= 5
               && length(table(labels_patient[temp.training_index.p2p])) == 2 && min(table(labels_patient[temp.training_index.p2p])) >= 5)) {
            break
          }
        }
        
        stop("Something went wrong during partitioning!")
      }
      
      stopifnot(length(temp.training_index.p2p) == training_amount.p)
      stopifnot(length(intersect(temp.test_index, temp.training_index.p2p)) == 0)
      stopifnot(temp.test_index %in% temp.patients)
      if (length(temp.cell_lines$slope) > 0) {
        stopifnot(length(temp.training_index.cp2p$slope) == training_amount.p + training_amount.c)
        stopifnot(length(intersect(temp.test_index, temp.training_index.cp2p$slope)) == 0)
        stopifnot(length(intersect(temp.test_index, temp.training_index.c2p$slope)) == 0)
      }
      if (length(temp.cell_lines$AUC) > 0) {
        stopifnot(length(temp.training_index.cp2p$AUC) == training_amount.p + training_amount.c)      
        stopifnot(length(intersect(temp.test_index, temp.training_index.cp2p$AUC)) == 0)      
        stopifnot(length(intersect(temp.test_index, temp.training_index.c2p$AUC)) == 0)
      }
      if (length(temp.cell_lines$IC50) > 0) {
        stopifnot(length(temp.training_index.cp2p$IC50) == training_amount.p + training_amount.c)
        stopifnot(length(intersect(temp.test_index, temp.training_index.cp2p$IC50)) == 0)
        stopifnot(length(intersect(temp.test_index, temp.training_index.c2p$IC50)) == 0)
      }
      
      partition$p2p[[temp.run_ind]] <- list(test_index = temp.test_index
                                           , training_index = temp.training_index.p2p)
      
      if (length(temp.cell_lines$slope) > 0) {
        partition$c2p.slope[[temp.run_ind]] <- list(test_index = temp.test_index
                                             , training_index = temp.training_index.c2p$slope)
        partition$cp2p.slope[[temp.run_ind]] <- list(test_index = temp.test_index
                                             , training_index = temp.training_index.cp2p$slope)
      }
      if (length(temp.cell_lines$AUC) > 0) {
        partition$c2p.AUC[[temp.run_ind]] <- list(test_index = temp.test_index
                                                  , training_index = temp.training_index.c2p$AUC)
        partition$cp2p.AUC[[temp.run_ind]] <- list(test_index = temp.test_index
                                                     , training_index = temp.training_index.cp2p$AUC)
      }
      if (length(temp.cell_lines$IC50) > 0) {
        partition$c2p.IC50[[temp.run_ind]] <- list(test_index = temp.test_index
                                                   , training_index = temp.training_index.c2p$IC50)
        partition$cp2p.IC50[[temp.run_ind]] <- list(test_index = temp.test_index
                                                     , training_index = temp.training_index.cp2p$IC50)
      }
    }
    
    for (temp.run_ind in 1:99) {
      stopifnot((partition$p2p[[temp.run_ind]]$test_index == partition$p2p[[temp.run_ind + 1]]$test_index) != length(partition$p2p[[temp.run_ind + 1]]$test_index))
      stopifnot((partition$p2p[[temp.run_ind]]$training_index == partition$p2p[[temp.run_ind + 1]]$training_index) != length(partition$p2p[[temp.run_ind + 1]]$training_index))
      
      if (length(temp.cell_lines$slope) > 0) {
        stopifnot((partition$c2p.slope[[temp.run_ind]]$training_index == partition$c2p.slope[[temp.run_ind + 1]]$training_index) != length(partition$c2p.slope[[temp.run_ind + 1]]$training_index))
        stopifnot((partition$cp2p.slope[[temp.run_ind]]$training_index == partition$cp2p.slope[[temp.run_ind + 1]]$training_index) != length(partition$cp2p.slope[[temp.run_ind + 1]]$training_index))
      }
      if (length(temp.cell_lines$AUC) > 0) {
        stopifnot((partition$c2p.AUC[[temp.run_ind]]$training_index == partition$c2p.AUC[[temp.run_ind + 1]]$training_index) != length(partition$c2p.AUC[[temp.run_ind + 1]]$training_index))
        stopifnot((partition$cp2p.AUC[[temp.run_ind]]$training_index == partition$cp2p.AUC[[temp.run_ind + 1]]$training_index) != length(partition$cp2p.AUC[[temp.run_ind + 1]]$training_index))
      }
      if (length(temp.cell_lines$IC50) > 0) {
        stopifnot((partition$c2p.IC50[[temp.run_ind]]$training_index == partition$c2p.IC50[[temp.run_ind + 1]]$training_index) != length(partition$c2p.IC50[[temp.run_ind + 1]]$training_index))
        stopifnot((partition$cp2p.IC50[[temp.run_ind]]$training_index == partition$cp2p.IC50[[temp.run_ind + 1]]$training_index) != length(partition$cp2p.IC50[[temp.run_ind + 1]]$training_index))
      }
    }
  
  } else {
    warning("Warning doing leave-one-out partitioning")
    for (temp.run_ind in 1:length(labels_patient)) {
      
      # if the training and test set does not contain at least 5 samples of each label, then resample
      temp.loop_count = 0
      temp.training_index.cp2p <- list()
      temp.training_index.c2p <- list()
      temp.training_index.p2p <- vector()
      
      temp.test_index <- temp.patients[temp.run_ind]
      temp.training_index.p2p <- setdiff(temp.patients, temp.test_index)
      stopifnot(length(temp.training_index.p2p) == training_amount.p)
      
      if (length(temp.cell_lines$slope) > 0) {
        temp.training_index.c2p$slope = temp.cell_lines$slope[cell_line_order$slope][1:training_amount.c]
        temp.training_index.cp2p$slope <- c(temp.training_index.c2p$slope, temp.training_index.p2p)
        stopifnot(length(temp.training_index.c2p$slope) == training_amount.c)
      }
      if (length(temp.cell_lines$AUC) > 0) {
        temp.training_index.c2p$AUC = temp.cell_lines$AUC[cell_line_order$AUC][1:training_amount.c]
        temp.training_index.cp2p$AUC <- c(temp.training_index.c2p$AUC, temp.training_index.p2p)
        stopifnot(length(temp.training_index.c2p$AUC) == training_amount.c)
      }
      if (length(temp.cell_lines$IC50) > 0) {
        temp.training_index.c2p$IC50 = temp.cell_lines$IC50[cell_line_order$IC50][1:training_amount.c]
        temp.training_index.cp2p$IC50 <- c(temp.training_index.c2p$IC50, temp.training_index.p2p)
        stopifnot(length(temp.training_index.c2p$IC50) == training_amount.c)
      }

      stopifnot(length(temp.training_index.p2p) == training_amount.p)
      stopifnot(length(intersect(temp.test_index, temp.training_index.p2p)) == 0)
      stopifnot(temp.test_index %in% temp.patients)
      
      if (length(temp.cell_lines$slope) > 0) {
        stopifnot(length(temp.training_index.cp2p$slope) == training_amount.p + training_amount.c)
        stopifnot(length(intersect(temp.test_index, temp.training_index.cp2p$slope)) == 0)
        stopifnot(length(intersect(temp.test_index, temp.training_index.c2p$slope)) == 0)
      }
      if (length(temp.cell_lines$AUC) > 0) {
        stopifnot(length(temp.training_index.cp2p$AUC) == training_amount.p + training_amount.c)      
        stopifnot(length(intersect(temp.test_index, temp.training_index.cp2p$AUC)) == 0)      
        stopifnot(length(intersect(temp.test_index, temp.training_index.c2p$AUC)) == 0)
      }
      if (length(temp.cell_lines$IC50) > 0) {
        stopifnot(length(temp.training_index.cp2p$IC50) == training_amount.p + training_amount.c)
        stopifnot(length(intersect(temp.test_index, temp.training_index.cp2p$IC50)) == 0)
        stopifnot(length(intersect(temp.test_index, temp.training_index.c2p$IC50)) == 0)
      }
      
      partition$p2p[[temp.run_ind]] <- list(test_index = temp.test_index
                                            , training_index = temp.training_index.p2p)
      
      if (length(temp.cell_lines$slope) > 0) {
        partition$c2p.slope[[temp.run_ind]] <- list(test_index = temp.test_index
                                                    , training_index = temp.training_index.c2p$slope)
        partition$cp2p.slope[[temp.run_ind]] <- list(test_index = temp.test_index
                                                     , training_index = temp.training_index.cp2p$slope)
      }
      if (length(temp.cell_lines$AUC) > 0) {
        partition$c2p.AUC[[temp.run_ind]] <- list(test_index = temp.test_index
                                                  , training_index = temp.training_index.c2p$AUC)
        partition$cp2p.AUC[[temp.run_ind]] <- list(test_index = temp.test_index
                                                   , training_index = temp.training_index.cp2p$AUC)
      }
      if (length(temp.cell_lines$IC50) > 0) {
        partition$c2p.IC50[[temp.run_ind]] <- list(test_index = temp.test_index
                                                   , training_index = temp.training_index.c2p$IC50)
        partition$cp2p.IC50[[temp.run_ind]] <- list(test_index = temp.test_index
                                                    , training_index = temp.training_index.cp2p$IC50)
      }
    }
    
    for (temp.run_ind in 1:(length(labels_patient) - 1)) {
      stopifnot((partition$p2p[[temp.run_ind]]$test_index == partition$p2p[[temp.run_ind + 1]]$test_index) != length(partition$p2p[[temp.run_ind + 1]]$test_index))
      stopifnot((partition$p2p[[temp.run_ind]]$training_index == partition$p2p[[temp.run_ind + 1]]$training_index) != length(partition$p2p[[temp.run_ind + 1]]$training_index))
      
      if (length(temp.cell_lines$slope) > 0) {
        stopifnot((partition$c2p.slope[[temp.run_ind]]$training_index == partition$c2p.slope[[temp.run_ind + 1]]$training_index) != length(partition$c2p.slope[[temp.run_ind + 1]]$training_index))
        stopifnot((partition$cp2p.slope[[temp.run_ind]]$training_index == partition$cp2p.slope[[temp.run_ind + 1]]$training_index) != length(partition$cp2p.slope[[temp.run_ind + 1]]$training_index))
      }
      if (length(temp.cell_lines$AUC) > 0) {
        stopifnot((partition$c2p.AUC[[temp.run_ind]]$training_index == partition$c2p.AUC[[temp.run_ind + 1]]$training_index) != length(partition$c2p.AUC[[temp.run_ind + 1]]$training_index))
        stopifnot((partition$cp2p.AUC[[temp.run_ind]]$training_index == partition$cp2p.AUC[[temp.run_ind + 1]]$training_index) != length(partition$cp2p.AUC[[temp.run_ind + 1]]$training_index))
      }
      if (length(temp.cell_lines$IC50) > 0) {
        stopifnot((partition$c2p.IC50[[temp.run_ind]]$training_index == partition$c2p.IC50[[temp.run_ind + 1]]$training_index) != length(partition$c2p.IC50[[temp.run_ind + 1]]$training_index))
        stopifnot((partition$cp2p.IC50[[temp.run_ind]]$training_index == partition$cp2p.IC50[[temp.run_ind + 1]]$training_index) != length(partition$cp2p.IC50[[temp.run_ind + 1]]$training_index))
      }
    }
  }
      
  return(partition)
}
