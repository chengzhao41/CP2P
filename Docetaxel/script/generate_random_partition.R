generate_random_partition <- function(input_labels_cell_lines, input_labels_patient, test_amount) {

  require("doParallel")
  
  partition <- list()

  temp.cell_lines <- 1:length(input_labels_cell_lines)
  temp.patients <- length(input_labels_cell_lines) + 1:length(input_labels_patient)
  
  partition$patients = temp.patients
  partition$cell_lines = temp.cell_lines

  ## cpp
  cpp.partition <- foreach (temp.run_ind = 1:length(input_labels_patient), .errorhandling="stop") %dopar% {        
    
    temp.test_index <- temp.patients[temp.run_ind]
    temp.training_index.single <- c(temp.cell_lines, setdiff(temp.patients, temp.test_index))          
    stopifnot(length(temp.test_index) > 0)
    stopifnot(length(temp.training_index.single) > 0)    
    stopifnot(length(intersect(temp.test_index, temp.training_index.single)) == 0)
    stopifnot(temp.test_index %in% temp.patients)
    
    list(
      test_index = temp.test_index
      , training_index.single = temp.training_index.single    
    )
  }  
  partition$cpp = cpp.partition  
  
  for (temp.run_ind in 1:(length(input_labels_patient) - 1)) {
    stopifnot(sum(cpp.partition[[temp.run_ind]]$test_index == cpp.partition[[temp.run_ind + 1]]$test_index) != length(cpp.partition[[temp.run_ind + 1]]$test_index))
    stopifnot(sum(cpp.partition[[temp.run_ind]]$training_index.single == cpp.partition[[temp.run_ind + 1]]$training_index.single) != length(cpp.partition[[temp.run_ind + 1]]$training_index.single))
  }
  
  ## cp
  cp.partition <- foreach (temp.run_ind = 1:100, .errorhandling="stop") %dopar% {    
    
    temp.test_index <- temp.patients
    temp.training_index.single <- sample(temp.cell_lines, length(temp.cell_lines) * 0.9)
    stopifnot(length(temp.test_index) > 0)
    stopifnot(length(temp.training_index.single) > 0)
    stopifnot(length(intersect(temp.test_index, temp.training_index.single)) == 0)    
    stopifnot(temp.test_index %in% temp.patients)
    
    list(
      test_index = temp.test_index
      , training_index.single = temp.training_index.single    
    )
  }
  partition$cp = cp.partition  

  for (temp.run_ind in 1:99) {    
    stopifnot(sum(cp.partition[[temp.run_ind]]$training_index.single == cp.partition[[temp.run_ind + 1]]$training_index.single) != length(cp.partition[[temp.run_ind + 1]]$training_index.single))
  }  
  
  ## pp
  temp.patients2 = temp.patients - length(input_labels_cell_lines)  
  pp.partition <- foreach (temp.run_ind = 1:length(input_labels_patient), .errorhandling="stop") %dopar% {    
    
    temp.test_index <- cpp.partition[[temp.run_ind]]$test_index - length(input_labels_cell_lines) 
    temp.training_index.single <- setdiff(temp.patients2, temp.test_index)
    
    stopifnot(length(temp.test_index) > 0)
    stopifnot(length(temp.training_index.single) > 0)    
    stopifnot(length(intersect(temp.test_index, temp.training_index.single)) == 0)    
    stopifnot(temp.test_index %in% temp.patients2)
    
    list(
      test_index = temp.test_index
      , training_index.single = temp.training_index.single    
    )
  }
  partition$pp = pp.partition  
  
  for (temp.run_ind in 1:(length(input_labels_patient) - 1)) {
    stopifnot(sum(pp.partition[[temp.run_ind]]$test_index == pp.partition[[temp.run_ind + 1]]$test_index) != length(pp.partition[[temp.run_ind + 1]]$test_index))
    stopifnot(sum(pp.partition[[temp.run_ind]]$training_index.single == pp.partition[[temp.run_ind + 1]]$training_index.single) != length(pp.partition[[temp.run_ind + 1]]$training_index.single))
  }   
  
  ## cc
  cc.partition <- foreach (temp.run_ind = 1:100, .errorhandling="stop") %dopar% {    
    
    temp.test_index <- sample(temp.cell_lines, test_amount$cc)
    temp.training_index.single <- setdiff(temp.cell_lines, temp.test_index)
    
    stopifnot(length(temp.test_index) > 0)
    stopifnot(length(temp.training_index.single) > 0)    
    stopifnot(length(intersect(temp.test_index, temp.training_index.single)) == 0)    
    stopifnot(temp.test_index %in% temp.cell_lines)
    
    list(
      test_index = temp.test_index
      , training_index.single = temp.training_index.single    
    )
  }
  partition$cc = cc.partition  
  
  for (temp.run_ind in 1:99) {
    stopifnot(sum(cc.partition[[temp.run_ind]]$test_index == cc.partition[[temp.run_ind + 1]]$test_index) != length(cc.partition[[temp.run_ind + 1]]$test_index))
    stopifnot(sum(cc.partition[[temp.run_ind]]$training_index.single == cc.partition[[temp.run_ind + 1]]$training_index.single) != length(cc.partition[[temp.run_ind + 1]]$training_index.single))
  }   
  
  
  return(partition)
}

generate_random_partition.cp_var <- function(input_labels_cell_lines, input_labels_patient, training_amount) {
  
  require("doParallel")
  stopifnot(input_labels_cell_lines < training_amount)  
  
  temp.cell_lines <- 1:length(input_labels_cell_lines)
  temp.patients <- length(input_labels_cell_lines) + 1:length(input_labels_patient)    
  input_combined_labels <- c(input_labels_cell_lines, input_labels_patient)
  
  cp.partition <- foreach (temp.run_ind = 1:100, .errorhandling="stop") %dopar% {        
    temp.test_index <- temp.patients
    temp.training_index.single <- sample(temp.cell_lines, training_amount)    
    while ((length(table(input_combined_labels[temp.training_index.single])) != 2)) {
      temp.test_index <- temp.patients
      temp.training_index.single <- sample(temp.cell_lines, training_amount)
    }
    
    stopifnot(length(temp.test_index) > 0)
    stopifnot(length(temp.training_index.single) > 0)
    stopifnot(length(intersect(temp.test_index, temp.training_index.single)) == 0)    
    stopifnot(temp.test_index %in% temp.patients)
    
    list(
      test_index = temp.test_index
      , training_index.single = temp.training_index.single    
    )
  }  
  
  for (temp.run_ind in 1:99) {    
    stopifnot(sum(cp.partition[[temp.run_ind]]$training_index.single == cp.partition[[temp.run_ind + 1]]$training_index.single) < length(cp.partition[[temp.run_ind + 1]]$training_index.single))
  }      
  
  return(cp.partition)
}

generate_random_partition.cpp_var <- function(input_labels_cell_lines, input_labels_patient, training_amount) {
  
  require("doParallel")
  stopifnot(length(input_labels_patient) > training_amount)  
  
  temp.cell_lines <- 1:length(input_labels_cell_lines)
  temp.patients <- length(input_labels_cell_lines) + 1:length(input_labels_patient)    
  input_combined_labels <- c(input_labels_cell_lines, input_labels_patient)
  
  cpp.partition <- foreach (temp.run_ind = 1:100, .errorhandling="stop") %dopar% {        
    
    temp.patients_in_training <- sample(temp.patients, training_amount)
    temp.training_index.single <- c(temp.cell_lines, temp.patients_in_training)
    temp.test_index <- setdiff(temp.patients, temp.patients_in_training)
    
    while (length(table(input_combined_labels[temp.training_index.single])) != 2) {
      temp.patients_in_training <- sample(temp.patients, training_amount)
      temp.training_index.single <- c(temp.cell_lines, temp.patients_in_training)
      temp.test_index <- setdiff(temp.patients, temp.patients_in_training)      
    }
    
    stopifnot(length(temp.test_index) > 0)
    stopifnot(length(temp.training_index.single) > 0)    
    stopifnot(length(intersect(temp.test_index, temp.training_index.single)) == 0)
    stopifnot(temp.test_index %in% temp.patients)
    
    list(
      test_index = temp.test_index
      , training_index.single = temp.training_index.single    
    )
  }    
  
  return(cpp.partition)
}

generate_random_partition.pp_var <- function(input_labels_patient, training_amount) {
  
  require("doParallel")
  stopifnot(length(input_labels_patient) > training_amount)   
  
  temp.patients <- 1:length(input_labels_patient)  
  
  pp.partition <- foreach (temp.run_ind = 1:100, .errorhandling="stop") %dopar% {        
    
    temp.training_index.single <- sample(temp.patients, training_amount)    
    temp.test_index <- setdiff(temp.patients, temp.training_index.single)
    
    while ((length(table(input_labels_patient[temp.training_index.single])) != 2)) {
      temp.training_index.single <- sample(temp.patients, training_amount)    
      temp.test_index <- setdiff(temp.patients, temp.training_index.single)
    }
    
    stopifnot(length(temp.test_index) > 0)
    stopifnot(length(temp.training_index.single) > 0)    
    stopifnot(length(intersect(temp.test_index, temp.training_index.single)) == 0)
    stopifnot(temp.test_index %in% temp.patients)
    
    list(
      test_index = temp.test_index
      , training_index.single = temp.training_index.single    
    )
  }    
  
  return(pp.partition)
}
