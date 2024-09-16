## This files contains all kernel functions, functions created for summarizing the results, and coverage ratio function

## warp function
warp.fun <- function(t,a=2.5,b=-1.2,c=2.5){
  omega <- 2*c*(-0.5 + 1/(1+ exp(-a*(t-b))))
  
  # Calculate the minimum and maximum values of omega
  omega_min <- min(omega)
  omega_max <- max(omega)

  # Desired interval
  desired_min <- -4.3
  desired_max <- 1.6

  # Linearly scale omega to the desired interval
  omega_scaled <- ((omega - omega_min) / (omega_max - omega_min)) * (desired_max - desired_min) + (desired_min)

  return(omega_scaled)
}


## Squared Exponential Kernel
SE <- function(X, sigma = 1, rho = median(as.matrix(dist(t(X)))), 
              specific_value=NULL,Identity_row = NULL,X_row,jitter = 1e-10){

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


# Period Kernel
Period <- function(X, sigma = 1, rho = median(as.matrix(dist(t(X)))), p = 24,
                  specific_value = NULL, Identity_row = NULL,X_row, jitter = 1e-10){

  if (is.null(specific_value) & is.null(Identity_row)){
      dist <- as.matrix(dist(t(X[X_row,,drop=FALSE])))
      G <- sigma^2 * exp(-2*sin(pi*dist/p)^2/(rho^2)) + jitter * diag(ncol(dist))
  }
  else {
      block <- t(t(X[Identity_row,]))%*%t(X[Identity_row,]) # block matrix
      dist <- as.matrix(dist(t(X[X_row,,drop=FALSE])))
      G <- sigma^2 * exp(-2*sin(pi*dist/p)^2/(rho^2))

      if (all(dim(block) != dim(G))){
        stop("The dimension of I_Matrix and G should be the same")
      }
      G <- block*G + jitter * diag(ncol(dist))
  }
  return(G)
}

# Linear Kernel
Linear <- function(X,sigma_b=1,sigma_v = 1,
                   specific_value = NULL,Identity_row = NULL,X_row, jitter = 1e-1) {
  if (is.null(specific_value) & is.null(Identity_row)){
   G <- sigma_v + sigma_b * (t(X[X_row,,drop=FALSE]) %*% X[X_row,,drop=FALSE]) 
   G <- G +  jitter * diag(ncol(G))
  }
  else {
    block <- t(t(X[Identity_row,]))%*%t(X[Identity_row,]) # block matrix
    G <- sigma_v + sigma_b * (t(X[X_row,,drop=FALSE]) %*% X[X_row,,drop=FALSE])
    if (all(dim(block) != dim(G))){
      stop("The dimension of I_Matrix and G should be the same")
    }
    G <- block*G + jitter * diag(ncol(G))
  }
  return(G)
}

# Rational Quadratic Kernel
# when a -> \infty, the SE kernel becomes the exponential kernel
RQ <- function(X,sigma =1, rho = median(as.matrix(dist(t(X)))), a = 1,
              specific_value=NULL,Identity_row = NULL,X_row,jitter = 1e-10){

  if (is.null(specific_value) & is.null(Identity_row)){
      dist <- as.matrix(dist(t(X[X_row,,drop=F])))
      G <- sigma^2 *(1+ (dist^2/(2 * a*rho^2)))^(-a) + jitter * diag(ncol(dist))
  }
  else {
      block <- t(t(X[Identity_row,]))%*%t(X[Identity_row,]) # block matrix
      dist <- as.matrix(dist(t(X[X_row,,drop=F])))
      G <- sigma^2 *(1+ (dist^2/(2 * a*rho^2)))^(-a) 

      if (all(dim(block) != dim(G))){
        stop("The dimension of I_Matrix and G should be the same")
      }
      G <- block*G + jitter * diag(ncol(dist))
  }
  return(G)
}



# Squared Exponential Kernel with warp function
# params: X, sigma, rho, specific_value, Identity_row, X_row, warp.fun, jitter
# warp.fun: function to warp the X, return Z = f(X)
SE.warp <- function(X,sigma =1, rho = median(as.matrix(dist(t(X)))),
                    specific_value = NULL, Identity_row=NULL,X_row, warp.fun = NULL,jitter = 1e-10){

  if (is.null(specific_value) & is.null(Identity_row)){
      dist <- as.matrix(dist(t(X[X_row,,drop=F])))
      G <- sigma^2 * exp(-dist^2/(2 * rho^2)) + jitter * diag(ncol(dist))
  }
  else {
  block <- t(t(X[Identity_row,]))%*%t(X[Identity_row,]) # block matrix
      if (is.null(warp.fun)){
        dist <- as.matrix(dist(t(X[X_row,,drop=F])))
        G <- sigma^2 * exp(-dist^2/(2 * rho^2))
        G <- block*G + + jitter * diag(ncol(dist))       
      }
      else {
        Z <- warp.fun(t(X[X_row,]))
        # x <- Z * X[Identity_row,]
        dist <- as.matrix(dist(t(Z)))
        G <- sigma^2 * exp(-dist^2/(2 * rho^2))
        G <- block*G  + jitter*diag(ncol(dist))
      }    
  }
  return(G)
}

# Rational Quadratic Kernel with warp function
RQ.warp <- function(X,sigma =1, rho = median(as.matrix(dist(t(X)))), a = 1,
              specific_value=NULL,Identity_row = NULL,X_row, warp.func = NULL,jitter = 1e-10){
  
  if (is.null(specific_value) & is.null(Identity_row)){
      dist <- as.matrix(dist(t(X[X_row,,drop=F])))
      G <- sigma^2 *(1+ (dist^2/(2 * a*rho^2)))^(-a) + jitter * diag(ncol(dist))
  }
  else {
    block <- t(t(X[Identity_row,]))%*%t(X[Identity_row,]) # block matrix
    if (is.null(warp.func)){
      dist <- as.matrix(dist(t(X[X_row,,drop=F])))
      G <- sigma^2 *(1+ (dist^2/(2 * a*rho^2)))^(-a) 
      G <- block* G + jitter * diag(ncol(dist))
    }
    else {
      Z <- warp.fun(t(X[X_row,]))
      dist <- as.matrix(dist(t(Z)))
      G <- sigma^2 *(1+ (dist^2/(2 * a*rho^2)))^(-a)
      G <- block* G + jitter * diag(ncol(dist))
    }
  }
  return(G)
}




# Function for summary lambda (CLR)
# 0. X: mod$Lambda[[num]]
# 1.filter_lambda: filter the lambda based on the filter_lambda, eg. filter_lambda = row of dummy lambda, must be defined as dummy
# 2.time_column: time column, time culumn of X
# 3.filter_coord: filter the coord based on the filter_coord, eg. filter_coord = 2, filter the lambda based on the coord = 2
# output: summary of specific lambda
# 4. CLR: if TRUE, then CLR, if FALSE, then ALR
Summary_function <- function(X,iteration = TRUE,filter_lambda,time_column,filter_coord, CLR = TRUE){
      if (CLR == TRUE){
          if (iteration == TRUE){
            summary <- gather_array(X,val,coord,sample,iter)%>%
              arrange(iter,coord) %>%
              mutate(time = rep(time_column,D*max(iter)))%>%
              mutate(filter = rep(filter_lambda,D*max(iter))) %>%
              filter(!is.na(val)) %>%
              filter(coord == filter_coord) %>%
              filter(filter == 1) %>%
              group_by(time) %>%
              summarise_posterior(val, na.rm=TRUE) %>%
              ungroup()
          }
          else{
            summary <- gather_array(X,val,coord,sample)%>%
              arrange(coord) %>%
              mutate(time = rep(time_column,D))%>%
              mutate(filter = rep(filter_lambda,D)) %>%
              filter(!is.na(val)) %>%
              filter(coord ==filter_coord) %>%
              filter(filter == 1) %>%
              group_by(time) %>%
              summarise_posterior(val, na.rm=TRUE) %>%
              ungroup()
          }
      }
      else if (CLR == FALSE){ #ALR instead
          if (iteration == TRUE){
            summary <- gather_array(X,val,coord,sample,iter)%>%
              arrange(iter,coord) %>%
              mutate(time = rep(time_column,(D-1)*max(iter)))%>%
              mutate(filter = rep(filter_lambda,(D-1)*max(iter))) %>%
              filter(!is.na(val)) %>%
              filter(coord == filter_coord) %>%
              filter(filter == 1) %>%
              group_by(time) %>%
              summarise_posterior(val, na.rm=TRUE) %>%
              ungroup()
          }
          else{
            summary <- gather_array(X,val,coord,sample)%>%
              arrange(coord) %>%
              mutate(time = rep(time_column,(D-1)))%>%
              mutate(filter = rep(filter_lambda,(D-1))) %>%
              filter(!is.na(val)) %>%
              filter(coord ==filter_coord) %>%
              filter(filter == 1) %>%
              group_by(time) %>%
              summarise_posterior(val, na.rm=TRUE) %>%
              ungroup()
          }
    }
return(summary)
}

standardize_function <-function(X,iteration = TRUE){
  if (iteration == TRUE){
    lambda <- array(dim =c(D-1,N,samplesize))
    for (j in 1:dim(X)[1]){
      for (i in 1:dim(X)[3]){
        lambda[j,,i] <- (X[j,,i]- mean(X[j,,i]))
        }
    }
  }
  else{
    lambda <- array(dim =c(D-1,N))
    for (j in 1:dim(X)[1]){
        lambda[j,] <- X[j,]- mean(X[j,])
        }
  }
  return(lambda)
}


## Function to calculate the ratio of coverage percentages between predicted and pseudo datasets.
## The coverage percentage is the proportion of true mean values that fall within the prediction intervals.
##
## Parameters:
## - true:   Data frame containing the true mean values (true$mean).
## - pred:   Data frame with predicted values, including p2.5 and p97.5 percentiles (pred$p2.5, pred$p97.5).
## - pesudo: Data frame similar to pred, representing NaddGP model predictions (pesudo$p2.5, pesudo$p97.5).
coverage.ratio <- function(true, pred,pesudo){
  pred.upper <- pred$p97.5
  pred.lower <- pred$p2.5
  pred$true <- true$mean
  pred$coverage <- ifelse(true$mean > pred.lower & true$mean < pred.upper, 1, 0)
  pred.coverage_percent <- sum(pred$coverage)/length(pred$coverage) 

  pesudo.upper <- pesudo$p97.5
  pesudo.lower <- pesudo$p2.5
  pesudo$true <- true$mean
  pesudo$coverage <- ifelse(true$mean > pesudo.lower & true$mean < pesudo.upper, 1, 0)
  pesudo.coverage_percent <- sum(pesudo$coverage)/length(pesudo$coverage) 

  if (pesudo.coverage_percent == 0) {
    warning("Pesudo coverage is 0. To compute the ratio, Pesudo coverage percent is setting as 0.01")
    pesudo.coverage_percent <- 0.01
    coverage.ratio <- pred.coverage_percent/pesudo.coverage_percent 
  }
  else {
     coverage.ratio <- pred.coverage_percent/pesudo.coverage_percent 
  }
  return(coverage.ratio) 
}