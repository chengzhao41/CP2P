# generates the random partitions with the test set common to all approaches
generate_random_partition <- function(
  input_labels_cell_lines, 
  input_labels_patient,
  training_amount.p,
  leave_one_out,
  cell_line_order,
  training_amount.c, 
  input_partition = NULL
  ) {
  
  stopifnot(length(input_labels_cell_lines) > 0)
  stopifnot(length(input_labels_cell_lines) == length(cell_line_order))
  stopifnot(length(input_labels_patient) >= 10)
  stopifnot(training_amount.p >= 10)
  stopifnot(!is.null(leave_one_out))
  stopifnot(length(input_labels_patient) > training_amount.p)
  stopifnot(length(input_labels_cell_lines$slope) >= training_amount.c)
  stopifnot(length(input_labels_cell_lines$AUC) >= training_amount.c)
  stopifnot(length(input_labels_cell_lines$IC50) >= training_amount.c)
  
  require("doParallel")

  temp.patients <- 1:length(input_labels_patient)
  
  temp.cell_lines <- list()
  temp.cell_lines$slope <- length(temp.patients) + 1:length(input_labels_cell_lines$slope)
  temp.cell_lines$AUC <- length(temp.patients) + 1:length(input_labels_cell_lines$AUC)
  temp.cell_lines$IC50 <- length(temp.patients) + 1:length(input_labels_cell_lines$IC50)
  
  partition <- list()
  
  for (temp.run_ind in 1:100) {
    
    # if the training and test set does not contain at least 5 samples of each label, then resample
    temp.loop_count = 0
    temp.training_index.cp2p <- list()
    temp.training_index.c2p <- list()
    temp.training_index.p2p <- vector()
    repeat {
      
      if (is.null(input_partition)) {
      temp.training_index.p2p <- sample(temp.patients, training_amount.p)
      stopifnot(length(temp.training_index.p2p) == training_amount.p)
      temp.test_index <- setdiff(temp.patients, temp.training_index.p2p)
      } else {
        temp.training_index.p2p <- input_partition[[temp.run_ind]]$training_index.p2p
        stopifnot(length(temp.training_index.p2p) == training_amount.p)
        temp.test_index <- input_partition[[temp.run_ind]]$test_index
      }
      
      if (length(temp.cell_lines$slope) > 0) {
        temp.training_index.c2p$slope = temp.cell_lines$slope[cell_line_order$slope][1:training_amount.c]
        temp.training_index.cp2p$slope <- c(temp.training_index.c2p$slope, temp.training_index.p2p)
      }
      if (length(temp.cell_lines$AUC) > 0) {
        temp.training_index.c2p$AUC = temp.cell_lines$AUC[cell_line_order$AUC][1:training_amount.c]
        temp.training_index.cp2p$AUC <- c(temp.training_index.c2p$AUC, temp.training_index.p2p)
      }
      if (length(temp.cell_lines$IC50) > 0) {
        temp.training_index.c2p$IC50 = temp.cell_lines$IC50[cell_line_order$IC50][1:training_amount.c]
        temp.training_index.cp2p$IC50 <- c(temp.training_index.c2p$IC50, temp.training_index.p2p)
      }
      
      stopifnot(length(temp.training_index.c2p$slope) == training_amount.c)
      stopifnot(length(temp.training_index.c2p$AUC) == training_amount.c)
      stopifnot(length(temp.training_index.c2p$IC50) == training_amount.c)
      
      if ((length(table(input_labels_patient[temp.test_index])) == 2 && min(table(input_labels_patient[temp.test_index])) >= 5
           && length(table(input_labels_patient[temp.training_index.p2p])) == 2 && min(table(input_labels_patient[temp.training_index.p2p])) >= 5)) {
        break
      }
      temp.loop_count = temp.loop_count + 1
      stopifnot(temp.loop_count < 100)
    }
    
    stopifnot(length(temp.training_index.p2p) == training_amount.p)
    stopifnot(length(temp.training_index.cp2p$slope) == training_amount.p + training_amount.c)
    stopifnot(length(temp.training_index.cp2p$AUC) == training_amount.p + training_amount.c)
    stopifnot(length(temp.training_index.cp2p$IC50) == training_amount.p + training_amount.c)
    stopifnot(length(intersect(temp.test_index, temp.training_index.p2p)) == 0)
    stopifnot(length(intersect(temp.test_index, temp.training_index.cp2p$slope)) == 0)
    stopifnot(length(intersect(temp.test_index, temp.training_index.cp2p$AUC)) == 0)
    stopifnot(length(intersect(temp.test_index, temp.training_index.cp2p$IC50)) == 0)
    stopifnot(length(intersect(temp.test_index, temp.training_index.c2p$slope)) == 0)
    stopifnot(length(intersect(temp.test_index, temp.training_index.c2p$AUC)) == 0)
    stopifnot(length(intersect(temp.test_index, temp.training_index.c2p$IC50)) == 0)
    stopifnot(temp.test_index %in% temp.patients)
    
    partition$p2p[[temp.run_ind]] <- list(test_index = temp.test_index
                                         , training_index = temp.training_index.p2p)
    partition$c2p.slope[[temp.run_ind]] <- list(test_index = temp.test_index
                                         , training_index = temp.training_index.c2p$slope)
    partition$c2p.AUC[[temp.run_ind]] <- list(test_index = temp.test_index
                                                , training_index = temp.training_index.c2p$AUC)
    partition$c2p.IC50[[temp.run_ind]] <- list(test_index = temp.test_index
                                                , training_index = temp.training_index.c2p$IC50)
    partition$cp2p.slope[[temp.run_ind]] <- list(test_index = temp.test_index
                                         , training_index = temp.training_index.cp2p$slope)
    partition$cp2p.AUC[[temp.run_ind]] <- list(test_index = temp.test_index
                                                 , training_index = temp.training_index.cp2p$AUC)
    partition$cp2p.IC50[[temp.run_ind]] <- list(test_index = temp.test_index
                                                 , training_index = temp.training_index.cp2p$IC50)
  }
  
  for (temp.run_ind in 1:99) {
    stopifnot((partition$p2p[[temp.run_ind]]$test_index == partition$p2p[[temp.run_ind + 1]]$test_index) != length(partition$p2p[[temp.run_ind + 1]]$test_index))
    stopifnot((partition$p2p[[temp.run_ind]]$training_index == partition$p2p[[temp.run_ind + 1]]$training_index) != length(partition$p2p[[temp.run_ind + 1]]$training_index))
    stopifnot((partition$c2p.slope[[temp.run_ind]]$training_index == partition$c2p.slope[[temp.run_ind + 1]]$training_index) != length(partition$c2p.slope[[temp.run_ind + 1]]$training_index))
    stopifnot((partition$c2p.AUC[[temp.run_ind]]$training_index == partition$c2p.AUC[[temp.run_ind + 1]]$training_index) != length(partition$c2p.AUC[[temp.run_ind + 1]]$training_index))
    stopifnot((partition$c2p.IC50[[temp.run_ind]]$training_index == partition$c2p.IC50[[temp.run_ind + 1]]$training_index) != length(partition$c2p.IC50[[temp.run_ind + 1]]$training_index))
    stopifnot((partition$cp2p.slope[[temp.run_ind]]$training_index == partition$cp2p.slope[[temp.run_ind + 1]]$training_index) != length(partition$cp2p.slope[[temp.run_ind + 1]]$training_index))
    stopifnot((partition$cp2p.AUC[[temp.run_ind]]$training_index == partition$cp2p.AUC[[temp.run_ind + 1]]$training_index) != length(partition$cp2p.AUC[[temp.run_ind + 1]]$training_index))
    stopifnot((partition$cp2p.IC50[[temp.run_ind]]$training_index == partition$cp2p.IC50[[temp.run_ind + 1]]$training_index) != length(partition$cp2p.IC50[[temp.run_ind + 1]]$training_index))
  }
  
  #partition$patients = temp.patients
  #partition$cell_lines = temp.cell_lines
  
  return(partition)
}
