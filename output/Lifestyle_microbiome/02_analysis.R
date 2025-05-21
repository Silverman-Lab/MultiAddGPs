
source('~/Lifestyle_microbiome/01_data_preprocessing.R', chdir = TRUE)


D <- nrow(Y)
N <- nrow(train_metadata_st)
C <- nrow(X)
samplesize <- 2000



upsilon <- D+ 100
Omega <- diag(D)
G <- cbind(diag(D-1),-1)
Xi <- (upsilon-D)*G%*%Omega%*%t(G)
# Xi <- matrix(50,D-1,D-1)
# diag(Xi) <- 100 



model <- function(Y,X,upsilon,Xi,D,
                  sigma_travel,rho_travel,
                  sigma_general,rho_general
                  ){
        
        set.seed(893)

        Gamma_travel <- function(X){SE(X,sigma = sigma_travel, rho = rho_travel,  X_row=3,jitter = 1e-6)}
        Gamma_general <- function(X){SE(X,sigma = sigma_general, rho = rho_general, X_row=2,jitter = 1e-6)}
        Theta_kernel <- function(X){matrix(0,D-1,ncol(X))}
        

        linear_component_theta <- matrix(0,D-1,1) 
        linear_component_gamma <- diag(1) + 1e-6 # a matrix with (row and col) dim of the length of linear components

        # Remember the order of putting the linear component is important, always put the linear component first
        Theta <- list(linear_component_theta,Theta_kernel,Theta_kernel)
        Gamma <- list(linear_component_gamma,Gamma_travel,Gamma_general)
        # Note: the order of output lambda depends on the order of the input lambda list

        fit <- fido::basset(Y,X,upsilon,Theta=Theta,Gamma=Gamma,Xi=Xi, linear =c(1),
                            samplesize = 2000,verbose = TRUE,seed = 893,
                            jitter = 1e-4, max_iter = 1000,
                            init = alr_array(Y+0.65, parts = 1))
}

MLL <- function(Y,X,upsilon, Xi,D,
                sigma_travel,rho_travel,
                sigma_general,rho_general
                ){
                    set.seed(893)
                    sigma_1 <- min(sigma_travel, sigma_general)
                    sigma_2 <- max(sigma_travel, sigma_general)

                    rho_1 <- min(rho_travel, rho_general)
                    rho_2 <- max(rho_travel, rho_general)


                    #prior for sigma and rho

                    prior_sigma_travel <- dnorm(sigma_travel, mean = 5, sd = 1, log=TRUE)
                    prior_sigma_general <- dnorm(sigma_general, mean = 5, sd = 1,log=TRUE)
                    prior_rho_travel <- dnorm(rho_travel, mean = 5, sd = 1,log=TRUE)
                    prior_rho_general <- dgamma(1 / rho_general, shape = 10, rate = 8.8, log=TRUE) - log(rho_general^2)

                    log_prior <- prior_sigma_travel + prior_sigma_general +                                                                                                                 
                             prior_rho_travel + prior_rho_general
                    prior_lambda <- 100

                    mod <- model(Y,X,upsilon,Xi,D,
                                sigma_travel = sigma_1,rho_travel = rho_2,
                                sigma_general = sigma_2,rho_general = rho_1
                                )
                    list(Score = mod$logMarginalLikelihood + log_prior*prior_lambda)
                }

# MLL(Y,X,upsilon,2,3,4,5)


set.seed(893)
optObj <- bayesOpt(
    FUN = function(sigma_travel,rho_travel,
                   sigma_general,rho_general
                  ){
                    MLL(Y,X,upsilon, Xi,D,
                        sigma_travel,rho_travel,
                        sigma_general,rho_general)},
                    bounds = list(
                        sigma_travel = c(2,10),
                        rho_travel = c(2,6),
                        sigma_general = c(1,5),
                        rho_general = c(0.5,2)
                        ),
                    initPoints = 10,
                    iters.n = 20,
                    verbose =1
)

optObj$scoreSummary
getBestPars(optObj)