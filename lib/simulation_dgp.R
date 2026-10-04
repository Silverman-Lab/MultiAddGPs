# Simulation 1 data-generating process.
#
# Defines the GP-simulated DGP used to generate multinomial count data.
# The design matrix contains the intercept, batch indicator, standardized
# Date, and raw time. The intervention kernel transforms raw time internally;
# no warped covariate is returned.

sim <- function(D=3,N=384,hyper_params=c(3,4,1,0.5),multi.depth=200,seed){
        set.seed(seed)
        D <- D # number of dimensions
        C <- 3 # number of covariates
        N <- N # number of samples

        time <- seq(0,N-1,by = 1)
        Date <- (time-mean(time))/sd(time)
        batch <- sample(rep(c("Batch1","Batch2"),each = N/2))

        X1 <- data.frame(batch = batch, time = time, Date = Date,
                         stringsAsFactors = FALSE)

        X <- t(model.matrix(~batch+Date+time,data = X1)) # Samples in columns

        hyper_params <- hyper_params

        ## Step 1: Sigma ~ IW(upsilon,Xi)
        Sigma <- matrix(0.5,D-1,D-1)
        diag(Sigma) <- 1

        ## Step 2:  f_3(X_3) ~ GP(0,\Sigma,\Gamma_3(X_3))  # intervention effect
        Gamma3 <- intervention_kernel(
          X = X,
          sigma = hyper_params[1],
          rho = hyper_params[2],
          X_row = "time",
          jitter = 1e-10
        )
        Theta3 <- matrix(0,D-1,N)
        x <- matrix(rnorm(N*(D-1)),D-1,byrow = TRUE)
        f3 <- Theta3 +t(chol(Sigma))%*%x%*%(chol(Gamma3))

        ## Step 3: f_2(X_2) ~ GP(0,\Sigma,\Gamma_2(X_2))  # general effect
        Gamma2 <- SE(X,sigma =hyper_params[3] , rho = hyper_params[4],X_row =3)
        Theta2 <- matrix(0,D-1,N)
        x <- matrix(rnorm(N*(D-1)),D-1,byrow = TRUE)
        f2 <- Theta2 +t(chol(Sigma))%*%x%*%(chol(Gamma2))

        ## Step 4: f_1 ~ MN(\theta_1,\Sigma,\Gamma_1)  # Batch effect
        beta_0 <- 0.4
        beta_1 <- 1
        f1 <- cbind(matrix(beta_0,D-1,1), matrix(beta_1,D-1,1))
        f1X1 <- f1%*%X[1:2,]

        ## Step 5: eta_.j ~ N(F(X),Sigma,I_N)
        #F(X) = f_1*X_1 + f_2(X_2) + f_3(X_3)
        F <- f1X1 + f2 + f3
        eta <- matrix(rnorm(N*(D-1)),D-1,byrow = TRUE)
        eta <-  F + t(chol(Sigma))%*%eta%*%(diag(N)*1)

        ## Step 6: p_.j = alrINV(eta_.j)
        pai <- t(alrInv(t(eta)))

        ## Step 5: y ~ Multinomial(500,\p)
        Y <- matrix(0, D, N)
        for (i in 1:N) Y[,i] <- rmultinom(1,multi.depth,prob = pai[,i])

        zerocount <- sum(Y==0)/length(Y)
        print(paste("The percentage of zeros in the data is",zerocount))

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
