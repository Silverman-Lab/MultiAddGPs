
source('/home/thcodelia/Odelia/Project/MultiAddGPs/src/helper_functions.R', chdir = TRUE)
source('/home/thcodelia/Odelia/Project/MultiAddGPs/src/00_package_loading.R', chdir = TRUE)



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

