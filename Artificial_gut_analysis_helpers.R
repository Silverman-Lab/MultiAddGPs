# Helpers used by artificial gut analysis.R

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
