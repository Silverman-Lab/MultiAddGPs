source('~/artificial_gut/01_data_preprocessing.R', chdir = TRUE)

D <- nrow(Y)
N <- nrow(mallard_family$sample_data)
C <- nrow(X)
samplesize <- 2000

# plot Y in clr space
Y_clr_tidy <- clr_array(Y+0.65, parts = 1) %>% 
  gather_array(mean, coord, sample) %>% 
  arrange(coord) %>%
  mutate(Hour_diff = rep(Hour_diff,D))%>%
  mutate(source = "Y")
  # filter(coord == 2)

# ggplot(Y_clr_tidy, aes(x = Hour_diff, y = mean)) +  
#     # geom_point()+
#     geom_point(alpha=0.4)+
#     geom_line(alpha=0.4)+
#   facet_wrap(~coord, scales="free_y",ncol=2)


## Standardize the input data


upsilon <- D+5
Model <- function(Y,X,upsilon,
                  sigma1,l1,
                  # nu,l2,
                  sigma2,l2,a2,
                  noise
                  ) {
    set.seed(893)
    Xi <- matrix(0.4,D-1,D-1)
    diag(Xi) <- noise 
    # Xi <- Xi + 1e-10
    # linear_component_theta <- matrix(0,D-1 ,nrow(X)-2)
    # linear_component_gamma <- diag(nrow(X)-2)

  # stationary kernel implementation
    Gamma_vessel1 <-function(X) {SE(X,sigma = sigma1, rho = l1, specific_value = 1,Identity_row = 1,X_row = C)}
    Gamma_vessel2 <-function(X) {SE(X,sigma = sigma1, rho = l1, specific_value = 1,Identity_row = 2,X_row = C)}
    Gamma_vessel3 <-function(X) {SE(X,sigma = sigma1, rho = l1, specific_value = 1,Identity_row = 3,X_row = C)}
    Gamma_vessel4 <-function(X) {SE(X,sigma = sigma1, rho = l1, specific_value = 1,Identity_row = 4,X_row = C)}
    Gamma_vessel1Spike <-function(X) {RQ(X,sigma = sigma2, rho = l2,a=a2, specific_value = 1,Identity_row = 7,X_row = C)}
    Gamma_vessel2Spike <-function(X) {RQ(X,sigma = sigma2, rho = l2,a=a2, specific_value = 1,Identity_row = 8,X_row = C)}

    # Gamma_vessel1Spike <-function(X) {RQ.warp(X,sigma =sigma2,rho = l2, a= a2,specific_value = 1,Identity_row = 7,X_row = C)}
    # Gamma_vessel2Spike <-function(X) {RQ.warp(X,sigma =sigma2,rho = l2, a= a2,specific_value = 1,Identity_row = 8,X_row = C)}

    Theta_kernel <- function(X) matrix(0, D-1, ncol(X))

    Theta <- list(Theta_kernel,Theta_kernel,Theta_kernel,Theta_kernel,Theta_kernel,Theta_kernel)
    Gamma <- list(Gamma_vessel1, Gamma_vessel2,Gamma_vessel3,Gamma_vessel4,Gamma_vessel1Spike,Gamma_vessel2Spike)

    mod <- fido::basset(Y, X, upsilon, Theta, Gamma, Xi, verbose = TRUE, seed = 893)
}


MLL <- function(Y,X,upsilon,
               sigma1,l1,
              #  nu,l2,
               sigma2,l2,a2,
               noise
               ){

  set.seed(893)
  sigma <- c(sigma1, sigma2)
  sigma <- sigma[order(sigma)]

  rho <- c(l1, l2)
  rho <- rho[order(rho)]

  # prior for length/sigma ~ InverseGamma(alpha,beta) 
  alpha1 <- 10
  beta1  <- 20 
  
  alpha2 <- 10
  beta2 <- 10 

  log_prior_l1 <- (alpha1 * log(beta1) - lgamma(alpha1) - (alpha1 + 1) * log(l1) - beta1/ l1)
  log_prior_l2 <- (alpha2 * log(beta2) - lgamma(alpha2) - (alpha2 + 1) * log(l2) - beta2/ l2)

  log_prior_s1 <- (alpha2 * log(beta2) - lgamma(alpha2) - (alpha2 + 1) * log(sigma1) - beta2/ sigma1)
  log_prior_s2 <- (alpha1 * log(beta1) - lgamma(alpha1) - (alpha1 + 1) * log(sigma2) - beta1/ sigma2)

  prior <- log_prior_l1 + log_prior_l2 + log_prior_s1 + log_prior_s2
  prior_lambda <- 120
# 

  mod <- Model(Y,X,upsilon,
               sigma1= sigma[1],l1 = rho[2],
              #  nu,l2,
               sigma2 = sigma[2],l2 = rho[1],a2,
               noise
               )

  list(Score = mod$logMarginalLikelihood + prior*prior_lambda)
}

set.seed(893)
optObj <- bayesOpt(
  FUN = function(sigma1,l1,
                #  nu,l2,
                 sigma2,l2
                # noise 
                 )
    MLL(Y,X,upsilon,
        sigma1,l1,
        # nu,l2,
        sigma2,l2,a2 =2,
        noise = 1
        ),
  bounds = list(
    sigma1 = c(0.1,2),
    l1 = c(1,3),
    sigma2 = c(1,4),
    l2 = c(0.1,2)
  ),
  initPoints = 10,
  iters.n = 20,
  verbose = 1
)

optObj$scoreSummary
getBestPars(optObj)
 