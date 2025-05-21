
source('/home/thcodelia/Odelia/Project/MultiAddGPs/src/helper_functions.R', chdir = TRUE)
source('/home/thcodelia/Odelia/Project/MultiAddGPs/src/package_loading.R', chdir = TRUE)



Simulation_linear <- function(D=4,N=500,hyper_params=c(3,5,1,20),multi.depth=200,seed=534536){
        set.seed(seed)
        D <- D # number of dimensions
        C <- 3 # number of covariates
        N <- N # number of samples

        time <- seq(0,N-1,by = 1)
        Date <- time
        time.s <- (time-mean(time))/sd(time)
        batch <- sample(rep(c("Batch1","Batch2"),each = N/2))
        X1 <- as.data.frame(cbind(batch, time.s)) %>% mutate(time = as.numeric(time.s))
        X1['Date'] = Date

        X <- t(model.matrix(~batch+time,data = X1)) # Transpose to get the right dimensions

        hyper_params <- hyper_params


        ## Step 1: Sigma ~ IW(upsilon,Xi)
        Sigma <- matrix(0.9,D-1,D-1)
        diag(Sigma) <- 1.5

        ## Step 2:  f_3(X_3) ~ GP(0,\Sigma,\Gamma_3(X_3))  # time effect_daily
        Gamma3 <- Period(X,sigma = hyper_params[1], rho = hyper_params[2], p = hyper_params[3], X_row =3)
        Theta3 <- matrix(0,D-1,N)
        x <- matrix(rnorm(N*(D-1)),D-1,byrow = TRUE)
        f3 <- Theta3 +t(chol(Sigma))%*%x%*%(chol(Gamma3))

        
        ## Linear Kernel instead
        ## Step 3: f_2(X_2) ~ GP(0,\Sigma,\Gamma_2(X_2))  # time effect_hourly\
        Gamma2 <- Linear(X,sigma_b =hyper_params[4], sigma_v = 1, X_row =3)
        Theta2 <- matrix(0,D-1,N)
        x <- matrix(rnorm(N*(D-1)),D-1,byrow = TRUE)
        f2 <- Theta2 +t(chol(Sigma))%*%x%*%(chol(Gamma2))


        ## Step 4: f_1 ~ MN(\theta_1,\Sigma,\Gamma_1)  # treatment effect
        beta_0 <- 2.7
        beta_1 <- 3
        f1 <- cbind(matrix(beta_0,D-1,1), matrix(beta_1,D-1,1))
        f1X1 <- f1%*%X[1:2,]


        ## Step 5: eta_.j ~ N(F(X),Sigma,I_N)
        F <- f1X1 + f2 + f3

        eta <- matrix(rnorm(N*(D-1)),D-1,byrow = TRUE)
        eta <-  F + t(chol(Sigma))%*%eta%*%(diag(N)*1)

        ## Step 6: p_.j = alr(eta_.j)
        pai <- t(alrInv(t(eta)))

        ## Step 5: y ~ Multinomial(500,\p)
        Y <- matrix(0, D, N)
        for (i in 1:N) Y[,i] <- rmultinom(1,multi.depth,prob = pai[,i])

        zerocount <- sum(Y==0)/length(Y)

        return(list(Y=Y,
                    X=X,
                    eta=eta,
                    pai= pai,
                    f1X1=f1X1,
                    f2=f2,
                    f3=f3,
                    F=F, 
                    Sigma=Sigma,
                    hyper_params=hyper_params,
                    beta_0=beta_0,
                    beta_1=beta_1,
                    X1=X1,
                    D = D,
                    N = N,
                    percent.zero = zerocount))
}

Simulation <- function( D=4, 
                        P=20, 
                        C=4, 
                        samplesize=2000,
                        hyper_params=c(4,25,30,1,30),
                        beta_0=2.7,
                        beta_1=1,
                        multi_depth=20,
                        seed=534536) {
    set.seed(seed)
    
    # Generate time sequence
    start_time <- as.POSIXct("2024-01-10 00:00:00", tz = "UTC")
    end_time <- as.POSIXct("2024-01-25 23:00:00", tz = "UTC")
    time_interval <- 60*60 # 1 hour
    Hour_time <- seq(start_time, end_time, by = time_interval)
    
    time <- Hour_time
    o <- order(time)
    
    Date <- rep(time[o], each = 1)
    time <- as.numeric(time[o])
    time <- (time - min(time))/3600
    
    # Generate batch assignments
    batch <- sample(rep(c("batch1", "batch2"), each = 72))
    X1 <- as.data.frame(cbind(batch, time)) %>% 
          mutate(time = as.numeric(time))
    X1['Date'] = Date
    
    # Create design matrix
    X <- t(model.matrix(~batch+time, data = X1))
    N <- ncol(X)
    
    ## Step 1: Generate Sigma
    Sigma <- matrix(0.9, D-1, D-1)
    diag(Sigma) <- 1.5
    
    ## Step 2: Generate f_3 (daily time effect)
    Gamma3 <- Period(X, sigma = hyper_params[1], rho = hyper_params[2], 
                     p = hyper_params[3], X_row = 3)
    Theta3 <- matrix(0, D-1, N)
    x <- matrix(rnorm(N*(D-1)), D-1, byrow = TRUE)
    f3 <- Theta3 + t(chol(Sigma)) %*% x %*% (chol(Gamma3))
    
    ## Step 3: Generate f_2 (hourly time effect)
    Gamma2 <- SE(X, sigma = hyper_params[4], rho = hyper_params[5], X_row = 3)
    Theta2 <- matrix(0, D-1, N)
    x <- matrix(rnorm(N*(D-1)), D-1, byrow = TRUE)
    f2 <- Theta2 + t(chol(Sigma)) %*% x %*% (chol(Gamma2))
    
    ## Step 4: Generate f_1 (treatment effect)
    f1 <- cbind(matrix(beta_0, D-1, 1), matrix(beta_1, D-1, 1))
    f1X1 <- f1 %*% X[1:2,]
    
    ## Step 5: Generate eta
    eta <- matrix(rnorm(N*(D-1)), D-1, byrow = TRUE)
    eta <- f1X1 + f2 + f3 + t(chol(Sigma)) %*% eta %*% (diag(N)*1)
    
    ## Step 6: Generate probabilities
    pai <- t(alrInv(t(eta)))
    
    ## Step 7: Generate counts
    Y <- matrix(0, D, N)
    for (i in 1:N) {
        Y[,i] <- rmultinom(1, multi_depth, prob = pai[,i])
    }
    
    # Calculate zero count percentage
    zerocount <- sum(Y==0)/length(Y)
    
    # Return all generated data and parameters
    return(list(
        Y = Y,
        X = X,
        eta = eta,
        pai = pai,
        f1X1 = f1X1,
        f2 = f2,
        f3 = f3,
        F = f1X1 + f2 + f3,
        Sigma = Sigma,
        hyper_params = hyper_params,
        beta_0 = beta_0,
        beta_1 = beta_1,
        X1 = X1,
        D = D,
        N = N,
        samplesize = samplesize,
        percent_zero = zerocount
    ))
}


# Simulation Function
#' @Description: this is a function that simulates data for a MultiAddGP model with three 
#' covariates including one fixed effect batch and two time effects (General and Intervention). 
#'
#' @param D  Number of taxa
#' @param N  Number of samples
#' @param hyper_params  A vector of hyperparameters for the GP model
#' @return function return a list of simulated data
#' @details 
#' @examples
Simulation_warp <- function(D=3,N=384,hyper_params=c(3,4,1,0.5),multi.depth=200,seed=5336){
        set.seed(seed)
        D <- D # number of dimensions
        C <- 3 # number of covariates
        N <- N # number of samples

        time <- seq(0,N-1,by = 1)
        Date <- (time-mean(time))/sd(time)
        batch <- sample(rep(c("Batch1","Batch2"),each = N/2))
        treatment <- warp(time, mode = 20, scale = 3, skew = 1)

        X1 <- as.data.frame(cbind(batch, time,treatment))%>% 
              mutate(time = as.numeric(time), treatment = as.numeric(treatment))
        X1['Date'] = Date

        X <- t(model.matrix(~batch+Date+treatment,data = X1)) # Transpose to get the right dimensions


        hyper_params <- hyper_params

        ## Step 1: Sigma ~ IW(upsilon,Xi)
        Sigma <- matrix(0.5,D-1,D-1)
        diag(Sigma) <- 1

        ## Step 2:  f_3(X_3) ~ GP(0,\Sigma,\Gamma_3(X_3))  # treatment_effect
        set.seed(5336)
        Gamma3 <- SE(X,sigma = hyper_params[1], rho = hyper_params[2],X_row = 4)
        Theta3 <- matrix(0,D-1,N)
        x <- matrix(rnorm(N*(D-1)),D-1,byrow = TRUE)
        f3 <- Theta3 +t(chol(Sigma))%*%x%*%(chol(Gamma3))
      #   par(mfrow=c(7,1))
      #   plot(X1$Date,f3[2,])
      #   points(f3[2,], type = "l", col = "red")
    
        ## Step 3: f_2(X_2) ~ GP(0,\Sigma,\Gamma_2(X_2))  # general effect
        set.seed(5336)
        Gamma2 <- SE(X,sigma =hyper_params[3] , rho = hyper_params[4],X_row =3)
        Theta2 <- matrix(0,D-1,N)
        x <- matrix(rnorm(N*(D-1)),D-1,byrow = TRUE)
        f2 <- Theta2 +t(chol(Sigma))%*%x%*%(chol(Gamma2))
      #   plot(X1$Date,f2[2,])
      #   plot(f2[2,]+f3[2,])

        ## Step 4: f_1 ~ MN(\theta_1,\Sigma,\Gamma_1)  # Batch effect
        set.seed(5336)
        beta_0 <- 0.4
        beta_1 <- 1
        f1 <- cbind(matrix(beta_0,D-1,1), matrix(beta_1,D-1,1))
        f1X1 <- f1%*%X[1:2,]
      #   plot(X1$Date,f1X1[2,])
        # plot(f1X1[2,]+f3[2,]+f2[2,])

        ## Step 5: eta_.j ~ N(F(X),Sigma,I_N)
        #F(X) = f_1*X_1 + f_2(X_2) + f_3(X_3)
        # noise <- 1
        F <- f1X1 + f2 + f3
      #   plot(F[2,])
        
        set.seed(5336)
        eta <- matrix(rnorm(N*(D-1)),D-1,byrow = TRUE)
        eta <-  F + t(chol(Sigma))%*%eta%*%(diag(N)*1)
      #   plot(X1$Date,eta[2,])


        ## Step 6: p_.j = alrINV(eta_.j)
        pai <- t(alrInv(t(eta)))
        # pai <- softmax(eta)

        ## Step 5: y ~ Multinomial(500,\p)
        Y <- matrix(0, D, N)
        for (i in 1:N) Y[,i] <- rmultinom(1,multi.depth,prob = pai[,i])
      #   plot(X1$Date,Y[1,])


        zerocount <- sum(Y==0)/length(Y)
      #   print(paste("The percentage of zeros in the data is",zerocount))

        return(list(Y=Y,
                    X=X,
                    eta=eta,
                    pai= pai,
                    f1X1=f1X1,
                    f2=f2,
                    f3=f3,
                    F=F, 
                    Sigma=Sigma,
                    hyper_params=hyper_params,
                    beta_0=beta_0,
                    beta_1=beta_1,
                    X1=X1,
                    D = D,
                    N = N,
                    depth = multi.depth,
                    percent.zero = zerocount))
}

