
# Squared-exponential covariance kernel.
#
# Builds the covariance matrix used to simulate GP components and fit
# MultiAddGP components. When `Identity_row` is supplied, the kernel is masked
# by the corresponding block structure before adding diagonal jitter.
SE <- function(X, sigma = 1, rho = median(as.matrix(dist(t(X)))), 
              specific_value=NULL, Identity_row = NULL, X_row, jitter = 1e-10){

  if (is.null(specific_value) & is.null(Identity_row)){
      dist <- as.matrix(dist(t(X[X_row,,drop=FALSE])))
      G <- sigma^2 * exp(-dist^2/(2 * rho^2)) + jitter * diag(ncol(dist))
  }
  else {
      block <- t(t(X[Identity_row,]))%*%t(X[Identity_row,]) # block matrix
      dist <- as.matrix(dist(t(X[X_row,,drop=FALSE])))
      G <- sigma^2 * exp(-dist^2/(2 * rho^2)) 

      if (all(dim(block) != dim(G))){
        stop("The dimension of block matrix and G should be the same")
      }
      G <- block*G + jitter * diag(ncol(dist))
  }
  return(G)
}

# Center component functions over samples.
#
# Used before comparing true and estimated component functions. For posterior
# arrays, centering is applied separately for each coordinate and iteration.
standardize_function <-function(X,iteration = TRUE){
  if (iteration == TRUE){
    lambda <- array(dim =c(dim(X)[1],dim(X)[2],dim(X)[3]))
    for (j in 1:dim(X)[1]){
      for (i in 1:dim(X)[3]){
        lambda[j,,i] <- (X[j,,i]- mean(X[j,,i]))
        }
    }
  }
  else{
    lambda <- array(dim =c(dim(X)[1],dim(X)[2]))
    for (j in 1:dim(X)[1]){
        lambda[j,] <- X[j,]- mean(X[j,])
        }
  }
  return(lambda)
}


# Intervention time transformation.
#
# Used internally by `intervention_kernel()`.
# The output is normalized to the interval [-1, 1].
warp <- function(x, mode = mode, scale = scale, skew = skew) {
  # Shift x so that mode is centered at zero
  x_shifted <- (x - mode) / scale
  
  # Apply a skewed Gaussian function with heavy right tail
  gaussian_core <- exp(-0.5 * x_shifted^2)  # Standard Gaussian curve
  skewed_tail <- 1 / (1 + exp(-skew * x_shifted))  # Logistic function for right-heavy tail
  transformed_x <- gaussian_core * skewed_tail
  
  # Normalize to output btw [-1, 1]
  transformed_x <- 2 * (transformed_x - min(transformed_x)) / (max(transformed_x) - min(transformed_x)) - 1
  
  return(transformed_x)
}

# Intervention covariance from ordinary elapsed time.
intervention_kernel <- function(X, sigma = 1, rho = median(as.matrix(dist(t(X)))),
                                specific_value = NULL, Identity_row = NULL,
                                X_row, jitter = 1e-10,
                                mode = 20, scale = 3, skew = 1) {
  if (length(X_row) != 1L) stop("X_row must select one raw-time row")
  z <- warp(X[X_row, ], mode = mode, scale = scale, skew = skew)
  # Apply the covariate's decimal-precision conversion.
  X[X_row, ] <- as.numeric(as.character(z))

  SE(
    X, sigma = sigma, rho = rho, specific_value = specific_value,
    Identity_row = Identity_row, X_row = X_row, jitter = jitter
  )
}

# Convert ALR-scale component functions to CLR scale.
alr_to_clr <- function(f,clr_D){
  f <- alrInv_array(f,clr_D,1)
  f <- clr_array(f,1)
  return(f)
}
