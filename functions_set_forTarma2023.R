Gram_K<-function(x,Kernel){
  if (Kernel=='Gauss'){
    sigma <- sqrt(0.5*median(dist(x)^2))
    X_gram <- dnorm(as.matrix(dist(x,diag=T,upper=T)), mean=0,sd=sigma)
  }
  if (Kernel=='Laplace'){
    sigma <- sqrt(0.5*median(dist(x)^2))
    X_gram <- exp(-as.matrix(dist(x,diag=T,upper=T))/sigma )
  }
  return(X_gram)
}

wild_boot_logrank<-function(x, z, d, kernel_x, kernel_z, num_bootstrap=1999){
  n <- nrow(x)
  p <- ncol(x)  
  
  # Sort the data in order of increasing time.
  sorted_indices <- order(z)
  z <- z[sorted_indices]
  d <- d[sorted_indices]
  x <- x[sorted_indices, ]
  
  # Define Y_matrix[i,:] to be the vector of indicators who are at risk at the i-th event time.
  Y_matrix <- matrix(0, n, n)
  Y_matrix[upper.tri(Y_matrix, diag = TRUE)] <- 1
  
  # Define Y[i] count the number of individuals at risk at the i-th event time.
  Y <- n - seq_len(n) + 1
  
  # Define A[i,:] to be a normalized (each row sums to 1) indicator of being at risk at time i. (note this is the transpose of A in our paper).
  scale_by_Y <- diag(1/Y)
  A <- scale_by_Y %*% Y_matrix
  
  # Define censoring_matrix[i,j] to be d[i]d[j]
  censoring_matrix <- outer(d, d)
  
  # Subtract A from the identity matrix
  I_minus_A <- diag(n) - A
  # Define the kernel matrix on X and Z
  K_X_Gram<-Gram_K(x,kernel_x)
  K_Z_Gram<-1#Gram_K(z,kernel_z)
  
  # Define Lz to be the kernel matrix on Z, with elementwise multiplication of the censoring matrix.
  Lz <- K_Z_Gram * censoring_matrix
  
  # Define the first_product matrix that we can re-use for computation in the wilde bootstrap.
  first_product <- I_minus_A%*%K_X_Gram%*%t(I_minus_A)
  original_statistic <- sum(first_product*Lz)
  #set.seed(1000)
  statistic_list <- rep(0,num_bootstrap)
  for (b in 1:num_bootstrap){
    W <- rbinom(n, 1, 1/2)*2-1
    WM <- outer(W, W)
    bootstrapLz <- WM * Lz
    multmatrix <- first_product * bootstrapLz
    bootstrap_statistic <- sum(multmatrix)
    statistic_list[b] <- bootstrap_statistic
  }
  pval <-   (1 + length(which(statistic_list>original_statistic))) / (1 + num_bootstrap)
  return(list('KLR'=original_statistic, 'pvalue'=pval))
}


SCMI<-function(x, z, d, kernel_x, kernel_z, num_bootstrap=1999){
  n <- nrow(x)
  p <- ncol(x)  
  
  # Sort the data in order of increasing time.
  sorted_indices <- order(z)
  z <- z[sorted_indices]
  d <- d[sorted_indices]
  x <- x[sorted_indices, ]
  
  # Define Y_matrix[i,:] to be the vector of indicators who are at risk at the i-th event time.
  Y_matrix <- matrix(0, n, n)
  Y_matrix[upper.tri(Y_matrix, diag = TRUE)] <- 1
  
  # Define Y[i] count the number of individuals at risk at the i-th event time.
  Y <- n - seq_len(n) + 1
  
  # Define A[i,:] to be a normalized (each row sums to 1) indicator of being at risk at time i. (note this is the transpose of A in our paper).
  scale_by_Y <- diag(1/Y)
  A <- scale_by_Y %*% Y_matrix
  
  # Define censoring_matrix[i,j] to be d[i]d[j]
  censoring_matrix <- outer(d, d)
  
  # Subtract A from the identity matrix
  I_minus_A <- diag(n) - A
  # Define the kernel matrix on X and Z
  K_Z_Gram<-1#Gram_K(z,kernel_z)
  
  # Define Lz to be the kernel matrix on Z, with elementwise multiplication of the censoring matrix.
  Lz <- K_Z_Gram * censoring_matrix
  
  # Define the first_product matrix that we can re-use for computation in the wilde bootstrap.
  first_product<-matrix(0,n,n)
  for(s in 1:p){
    K_X_Gram<-Gram_K(x[,s],kernel_x)
    first_product <- first_product+I_minus_A%*%K_X_Gram%*%t(I_minus_A)
  }
  
  original_statistic <- sum(first_product*Lz)
  #set.seed(1000)
  statistic_list <- rep(0,num_bootstrap)
  for (b in 1:num_bootstrap){
    W <- rbinom(n, 1, 1/2)*2-1
    WM <- outer(W, W)
    bootstrapLz <- WM * Lz
    multmatrix <- first_product * bootstrapLz
    bootstrap_statistic <- sum(multmatrix)
    statistic_list[b] <- bootstrap_statistic
  }
  pval <-   (1 + length(which(statistic_list>original_statistic))) / (1 + num_bootstrap)
  return(list('KLR'=original_statistic, 'pvalue'=pval))
}