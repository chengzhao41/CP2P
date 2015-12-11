# generate random partitions for a given a number of cell lines and patients
generate_random_partition.cpp_var2 <- function(
  input_labels_cell_lines = NULL, 
  input_labels_patient = NULL, 
  cell_line_training_amount = NULL, 
  patient_training_amount = NULL) {

  stopifnot(input_labels_cell_lines != NULL)
  stopifnot(input_labels_patient != NULL)
  stopifnot(cell_line_training_amount != NULL)
  stopifnot(patient_training_amount != NULL)
  
  require("doParallel")
  stopifnot(length(input_labels_patient) > patient_training_amount)
  stopifnot(length(input_labels_cell_lines) > cell_line_training_amount)
  
  temp.cell_lines <- 1:length(input_labels_cell_lines)
  temp.patients <- length(input_labels_cell_lines) + 1:length(input_labels_patient)    
  input_combined_labels <- c(input_labels_cell_lines, input_labels_patient)
  
  cpp.partition <- foreach (temp.run_ind = 1:100, .errorhandling="stop") %dopar% {        
    
    temp.patients_in_training <- sample(temp.patients, patient_training_amount)
    temp.cell_lines_in_training <- sample(temp.cell_lines, cell_line_training_amount)
    temp.training_index.single <- c(temp.cell_lines_in_training, temp.patients_in_training)
    temp.test_index <- setdiff(temp.patients, temp.patients_in_training)
    
    while ((length(table(input_combined_labels[temp.training_index.single])) != 2) || (length(table(input_combined_labels[temp.test_index])) != 2)) {
      temp.patients_in_training <- sample(temp.patients, patient_training_amount)
      temp.cell_lines_in_training <- sample(temp.cell_lines, cell_line_training_amount)
      temp.training_index.single <- c(temp.cell_lines_in_training, temp.patients_in_training)
      temp.test_index <- setdiff(temp.patients, temp.patients_in_training)     
    }
    
    stopifnot(length(temp.test_index) > 0)
    stopifnot(length(temp.training_index.single) == cell_line_training_amount + patient_training_amount)
    stopifnot(length(intersect(temp.test_index, temp.training_index.single)) == 0)
    stopifnot(temp.test_index %in% temp.patients)
    
    list(
      test_index = temp.test_index
      , training_index.single = temp.training_index.single    
    )
  }
  
  for (temp.run_ind in 1:99) {
    stopifnot(sum(cpp.partition[[temp.run_ind]]$test_index == cpp.partition[[temp.run_ind + 1]]$test_index) != length(cpp.partition[[temp.run_ind + 1]]$test_index))
    stopifnot(sum(cpp.partition[[temp.run_ind]]$training_index.single == cpp.partition[[temp.run_ind + 1]]$training_index.single) != length(cpp.partition[[temp.run_ind + 1]]$training_index.single))
  }
  
  return(cpp.partition)
}


# generate random partitions for a given a number of cell lines and patients
# cell line training data is added based on provided order
generate_random_partition.cpp_var3 <- function(
  input_labels_cell_lines = NULL, 
  input_labels_patient = NULL, 
  cell_line_training_amount = NULL, 
  patient_training_amount = NULL,
  cell_line_order = NULL) {

  stopifnot(input_labels_cell_lines != NULL)
  stopifnot(input_labels_patient != NULL)
  stopifnot(cell_line_training_amount != NULL)
  stopifnot(patient_training_amount != NULL)
  stopifnot(cell_line_order != NULL)  
  
  require("doParallel")
  stopifnot(length(input_labels_patient) > patient_training_amount)
  stopifnot(length(input_labels_cell_lines) > cell_line_training_amount)
  
  temp.cell_lines <- 1:length(input_labels_cell_lines)
  temp.patients <- length(input_labels_cell_lines) + 1:length(input_labels_patient)    
  input_combined_labels <- c(input_labels_cell_lines, input_labels_patient)
  
  cpp.partition <- foreach (temp.run_ind = 1:100, .errorhandling="stop") %dopar% {        
    
    temp.patients_in_training <- sample(temp.patients, patient_training_amount)
    temp.cell_lines_in_training <- temp.cell_lines[cell_line_order][1:cell_line_training_amount]
    temp.training_index.single <- c(temp.cell_lines_in_training, temp.patients_in_training)
    temp.test_index <- setdiff(temp.patients, temp.patients_in_training)
    
    while ((length(table(input_combined_labels[temp.training_index.single])) != 2) || (length(table(input_combined_labels[temp.test_index])) != 2)) {
      temp.patients_in_training <- sample(temp.patients, patient_training_amount)
      temp.training_index.single <- c(temp.cell_lines_in_training, temp.patients_in_training)
      temp.test_index <- setdiff(temp.patients, temp.patients_in_training)     
    }
    
    stopifnot(length(temp.test_index) > 0)
    stopifnot(length(temp.training_index.single) == cell_line_training_amount + patient_training_amount)    
    stopifnot(length(intersect(temp.test_index, temp.training_index.single)) == 0)
    stopifnot(temp.test_index %in% temp.patients)
    
    list(
      test_index = temp.test_index
      , training_index.single = temp.training_index.single    
    )
  }
  
  for (temp.run_ind in 1:99) {
    stopifnot(sum(cpp.partition[[temp.run_ind]]$test_index == cpp.partition[[temp.run_ind + 1]]$test_index) != length(cpp.partition[[temp.run_ind + 1]]$test_index))
    stopifnot(sum(cpp.partition[[temp.run_ind]]$training_index.single == cpp.partition[[temp.run_ind + 1]]$training_index.single) != length(cpp.partition[[temp.run_ind + 1]]$training_index.single))
  }
  
  return(cpp.partition)
}
