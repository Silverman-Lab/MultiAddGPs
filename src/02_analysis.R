
source('/home/thcodelia/Odelia/Project/MultiAddGPs/src/01_data_preprocessing.R', chdir = TRUE)

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

    # linear_component_theta <- matrix(beta_0,D-1,1)
    # linear_component_gamma <- diag(1)

  # stationary kernel implementation
    Gamma_vessel1 <-function(X) {SE(X,sigma = sigma1, rho = l1, specific_value = 1,Identity_row = 2,X_row = C)}
    Gamma_vessel2 <-function(X) {SE(X,sigma = sigma1, rho = l1, specific_value = 1,Identity_row = 3,X_row = C)}
    Gamma_vessel3 <-function(X) {SE(X,sigma = sigma1, rho = l1, specific_value = 1,Identity_row = 4,X_row = C)}
    Gamma_vessel4 <-function(X) {SE(X,sigma = sigma1, rho = l1, specific_value = 1,Identity_row = 5,X_row = C)}
    # Gamma_vessel1Spike <-function(X) {RQ(X,sigma = sigma2, rho = l2,a=a2, specific_value = 1,Identity_row = 5,X_row = C)}
    # Gamma_vessel2Spike <-function(X) {RQ(X,sigma = sigma2, rho = l2,a=a2, specific_value = 1,Identity_row = 6,X_row = C)}

    Gamma_vessel1Spike <-function(X) {RQ.warp(X,sigma =sigma2,rho = l2, a= a2,specific_value = 1,Identity_row = 8,X_row = C)}
    Gamma_vessel2Spike <-function(X) {RQ.warp(X,sigma =sigma2,rho = l2, a= a2,specific_value = 1,Identity_row = 9,X_row = C)}

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

    # prior for length ~ InverseGamma(alpha,beta) 
    # alpha <- 6.407
    # beta  <- 12.24
    alpha1 <- 20
    beta1  <- 42  # with model around 2
    
    alpha2 <- 20
    beta2 <- 21 # with mode around 1

    log_prior_l1 <- (alpha1 * log(beta1) - lgamma(alpha1) - (alpha1 + 1) * log(l1) - beta1/ l1)
    log_prior_l2 <- (alpha2 * log(beta2) - lgamma(alpha2) - (alpha2 + 1) * log(l2) - beta2/ l2)

    # # prior for sigma ~ InverseGamma(alpha,beta) 
    alpha_s1 <- 20
    beta_s1  <- 14.7 # with mode around 0.7

    # alpha_s2 <- 5
    # beta_s2 <- 18

    log_prior_s1 <- (alpha_s1 * log(beta_s1) - lgamma(alpha_s1) - (alpha_s1 + 1) * log(sigma1) - beta_s1/ sigma1)
    log_prior_s2 <- (alpha1 * log(beta1) - lgamma(alpha1) - (alpha1 + 1) * log(sigma2) - beta1/ sigma2)

    prior <- log_prior_l1 + log_prior_l2 + log_prior_s1 + log_prior_s2
    prior_lambda <- 100
    # 

    mod <- Model(Y,X,upsilon,
                sigma1= sigma[1],l1 = rho[2],
                #  nu,l2,
                sigma2 = sigma[2],l2 = rho[1],a2,
                noise
                )

    list(Score = mod$logMarginalLikelihood + prior*prior_lambda)
    }

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
        sigma1 = c(0.1,1),
        l1 = c(1,3),
        # a1 = c(1,20),
        # nu = c (0.1,1),
        # l2 = c(1,3),
        sigma2 = c(1,4),
        l2 = c(0.1,2)
        # a2 = c(10,20)
        # noise = c(1,20)
    ),
    initPoints = 10,
    iters.n = 20,
    verbose = 1
    )

optObj$scoreSummary
getBestPars(optObj)
 