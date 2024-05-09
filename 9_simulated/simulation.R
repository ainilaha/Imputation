# Load necessary library
library(MASS)


simulateData = function(nrows=1500){
  # Define the base data parameters
  mu = rep(0, 15)  # 15-dimensional normal distribution with mean 0
  Sigma = diag(15)  # Identity matrix as the covariance matrix (independent columns)
  
  # Generate initial data
  set.seed(123)
  data = mvrnorm(n = nrows, mu = mu, Sigma = Sigma)
  
  # Introduce non-linear relationships
  # Quadratic relationship
  data[, 2] = data[, 1]^2 + rnorm(nrows)
  
  # Exponential relationship
  data[, 3] = exp(data[, 2]/3) + rnorm(nrows)
  
  # Sine transformation
  data[, 4] = sin(data[, 3]) + rnorm(nrows)
  
  # Interaction term (product of two columns)
  data[, 5] = data[, 1] * data[, 2] + rnorm(nrows)
  
  # Logarithmic transformation (adding a constant to avoid log(0) issue)
  data[, 6] = log(abs(data[, 5]) + 1) + rnorm(nrows)
  
  # Additional non-linear transformations for more complexity
  data[, 7] = sqrt(abs(data[, 6])) + rnorm(nrows)
  data[, 8] = data[, 7]^3 + rnorm(nrows)
  data[, 9] = tanh(data[, 8]) + rnorm(nrows)
  data[, 10] = (data[, 9] / data[, 4])^2 + rnorm(nrows)
  data[, 11] = cos(data[, 10]) + rnorm(nrows)
  data[, 12] = abs(data[, 11] - data[, 5]) + rnorm(nrows)
  data[, 13] = data[, 12] * data[, 3] + rnorm(nrows)
  data[, 14] = (data[, 13] + data[, 6] / data[, 1])^2 + rnorm(nrows)
  data[, 15] = data[, 14] / data[, 11] + rnorm(nrows)
  
  # Create a dataframe and return
  data = as.data.frame(data)
  colnames(data) = paste0(rep("v",15),1:15)
  data
}


addMissingness = function(dataframe, ratio=0.1) {
  # Calculate the total number of data points in the dataframe
  total_elements <- nrow(dataframe) * ncol(dataframe)
  
  # Calculate the number of elements to set as NA
  num_missing <- ceiling(total_elements * ratio)
  
  # Generate random indices to set as NA
  missing_indices <- sample(total_elements, num_missing)
  
  # Convert matrix indices to array indices
  array_indices <- arrayInd(missing_indices, dim(dataframe))
  
  # Assign NA to random elements in the dataframe
  dataframe[array_indices] <- NA
  
  # Return a list containing both the modified dataframe and the random indices
  # return(list("missData" = dataframe, "missIndex" = array_indices))
  dataframe
}
