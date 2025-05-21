

source("Project/MultiAddGPs/src/helper_functions.R")
source("Project/MultiAddGPs/src/package_loading.R")
source("Project/MultiAddGPs/src/simulation_function.R")

simulation_data <- Simulation(D=4,P=20,C=4,
                              samplesize=2000,
                              hyper_params=c(4,25,30,1,30),
                              beta_0=2.7,
                              beta_1=1,
                              multi_depth=20,
                              seed=534536)

Y <- simulation_data$Y
X <- simulation_data$X
eta <- simulation_data$eta
pai <- simulation_data$pai
f1X1 <- simulation_data$f1X1
f2 <- simulation_data$f2
f3 <- simulation_data$f3
F <- simulation_data$F
Sigma <- simulation_data$Sigma
hyper_params <- simulation_data$hyper_params
beta_0 <- simulation_data$beta_0
beta_1 <- simulation_data$beta_1
X1 <- simulation_data$X1
D <- simulation_data$D
N <- simulation_data$N
samplesize <- simulation_data$samplesize
percent_zero <- simulation_data$percent_zero
Lambda_time <- X1$time



# NaddGPs count model: 
                              # Y_j+ 0.5  = alr(eta_j)
                              # eta ~ N(F(X),Sigma,I_n)
                              # F(X) ~  GP(beta_0 + \theta_1X_1 + Theta_2 + Theta_3,Sigma,X^T*Gamma1*X + Gamma2(X_2) + Gamma3(X_2)+I_n)  
                              # Sigma ~ IW(upsilon,Xi)

###################################### start here #############################

NaddGPs <- function(Y=Y,X=X,D=D,sigma1,l1,p,sigma2,l2,b0,b1,N){
      noise_eta <- 1
      Sigma <- diag(D-1)*2
      Theta1 <-  cbind(matrix(b0,D-1,1), matrix(b1,D-1,1))
      Gamma1 <- diag(nrow(X)-1)*0.1
      Gamma2 <- Period(X,sigma = sigma1, rho = l1, p = p, X_row = 3)
      Gamma3 <- SE(X,sigma = sigma2, rho = l2, X_row = 3)
      Theta_kernel <- matrix(0, D-1, ncol(X))
      # combine all theta and gamma
      Theta <-  Theta1%*%X[1:2,] + Theta_kernel + Theta_kernel
      Gamma <- t(X[1:2,])%*%Gamma1%*%X[1:2,] + Gamma2 + Gamma3
      # Transform Y to alr
      eta <- alr_array(Y+0.65, parts = 1)
      return(list(eta,Theta,Gamma,Sigma,noise_eta,N,Theta_kernel,Theta1,Gamma1,Gamma2,Gamma3))
}

############## Optimize the hyperparameters############################################ 
MLL_NaddGPs<- function(Y,X,Sigma,D,sigma1,l1,p,sigma2,l2,b0 = beta_0,b1 = beta_1,N=N){

      NaddGPs <- NaddGPs(Y,X,D,sigma1,l1,p,sigma2,l2,b0 = beta_0,b1 = beta_1,N)
      eta <- NaddGPs[[1]]
      Theta <- NaddGPs[[2]]
      Gamma <- NaddGPs[[3]]
      Sigma <- NaddGPs[[4]]
      noise_eta <- NaddGPs[[5]]
      N <- NaddGPs[[6]]
      # Analytical form Marginal likelihood
      set.seed(893)

      prior_sigma1 <- dnorm(sigma1,mean = 4,sd = 1,log = TRUE)
      prior_sigma2 <- dnorm(sigma2,mean = 1,sd = 1,log = TRUE)

      prior_l1 <- dnorm(l1,mean = 25,sd = 1,log = TRUE)
      prior_l2 <- dnorm(l2,mean = 30,sd = 1,log = TRUE) 
      
      prior_log <- prior_sigma1 + prior_sigma2 + prior_l1 + prior_l2
      prior_lambda <- 100

      logMarginalLikelihood <- dmvnorm(as.vector(eta),mean=as.vector(Theta),sigma = kronecker(Sigma,(Gamma+diag(N)*noise_eta)),log = TRUE)
      list(Score=logMarginalLikelihood + prior_log * prior_lambda)
}

set.seed(893)
optObj_NaddGPs <- bayesOpt(
    FUN= function(sigma1,l1,
                  sigma2,l2
                 ){
                  MLL_NaddGPs(Y,X,Sigma,D,
                           sigma1,l1,p = 30,
                           sigma2,l2,
                           b0 = beta_0,
                           b1 = beta_1,
                           N=N)},
                  bounds = list( 
                        sigma1 = c(1,hyper_params[1]+5),
                        l1 = c(hyper_params[2]-10,hyper_params[2]+10),
                        sigma2 = c(1,hyper_params[4]+5),
                        l2 = c(hyper_params[5]-10,hyper_params[5]+10)
                        ),
                  initPoints = 10,
                  iters.n = 40,
                  verbose = 1
)
optObj_NaddGPs$scoreSummary
getBestPars(optObj_NaddGPs)

########################################################################################

# Fit the hyperparameters back to the NaddGPs_model and Sample F

params <- getBestPars(optObj_NaddGPs)
# params <-hyper_params


NaddGPs_count_model <- function(Y,X,D,beta_0,beta_1,N,params,samplesize){
    p <- 3 # number of component in the model
    l <- 2 # number of linear component
    
    # Pre-allocate arrays
    sample_posterior <- lapply(1:(p+1), function(g) array(dim=c(D-1,N,samplesize)))
    
    # Get NaddGPs results
    NaddGPs <- NaddGPs(Y,X,D,
                 sigma1 = params$sigma1,
                 l1 = params$l1,
                 p = 30,
                 sigma2 = params$sigma2,
                 l2 = params$l2,
                 b0 = beta_0,
                 b1 = beta_1,
                 N = N)
    
    # Extract components
    eta <- NaddGPs[[1]]
    Theta <- NaddGPs[[2]]
    Gamma <- NaddGPs[[3]]
    Sigma <- NaddGPs[[4]]
    noise_eta <- NaddGPs[[5]]
    N <- NaddGPs[[6]]
    Theta_kernel <- NaddGPs[[7]]
    Theta1 <- NaddGPs[[8]]
    Gamma1 <- NaddGPs[[9]]
    Gamma2 <- NaddGPs[[10]]
    Gamma3 <- NaddGPs[[11]]
    
    # Pre-compute frequently used matrices
    noise_eta_inv <- solve(diag(N) * noise_eta)
    Gamma_inv <- solve(Gamma)
    Gamma1_inv <- solve(Gamma1)
    Gamma2_inv <- solve(Gamma2)
    Gamma3_inv <- solve(Gamma3)
    
    # Compute mean and covariance matrices
    mean_F <- (eta %*% noise_eta_inv + Theta %*% Gamma_inv) %*% solve(noise_eta_inv + Gamma_inv)
    row_cov <- Sigma
    col_cov <- solve(noise_eta_inv + Gamma_inv)
    
    # Transform Matrix normal to multivariate normal for sampling
    product <- kronecker(row_cov, col_cov)
    vec_mean <- as.vector(mean_F)
    
    # Generate all random samples at once
    stand <- matrix(rnorm(N * (D-1) * samplesize), (D-1) * N, samplesize)
    
    # Eigenvalue decomposition for kronecker product
    L_product <- Eigenvalue_decomp(product)
    
    # Sampling from F - vectorized
    Sam_F <- L_product %*% stand
    Sam_F <- sweep(Sam_F, 1, vec_mean, FUN=`+`)
    
    # Reorganize the sample F - vectorized
    Sample_F <- array(Sam_F, dim = c(D-1, N, samplesize))
    sample_posterior[[4]] <- Sample_F
    
    # Back sampling 
    for(i in 1:p) {
        if(i == 1) {

            K <- Gamma2 + Gamma3
            K_inv <- solve(K)
            X_sub <- X[1:2,]
            XKXt <- X_sub %*% K_inv %*% t(X_sub)
            solve_term <- solve(XKXt + Gamma1_inv)
            
            for(s in 1:samplesize) {
                Lambda_star <- Sample_F[,,s] - 2 * Theta_kernel
                mean <- (Lambda_star %*% K_inv %*% t(X_sub) + Theta1 %*% Gamma1_inv) %*% solve_term
                L <- Eigenvalue_decomp(solve_term)
                sample_i <- mean + t(chol(Sigma)) %*% matrix(rnorm(l * (D-1)), D-1, l) %*% L
                sample_posterior[[i]][,,s] <- sample_i %*% X_sub
            }
        } else if(i > 1 && i < p) {

            K_inv <- Gamma3_inv
            solve_term <- solve(K_inv + Gamma2_inv)
            
            for(s in 1:samplesize) {
                Lambda_star <- Sample_F[,,s] - sample_posterior[[i-1]][,,s] - Theta_kernel
                mean <- (Lambda_star %*% K_inv + Theta_kernel %*% Gamma2_inv) %*% solve_term
                L <- Eigenvalue_decomp(solve_term)
                sample_i <- mean + t(chol(Sigma)) %*% matrix(rnorm(N * (D-1)), D-1, N) %*% L
                sample_posterior[[i]][,,s] <- sample_i
            }
        } else if(i == p) {
            sample_posterior[[i]] <- Sample_F - sample_posterior[[i-1]] - sample_posterior[[i-2]]
        }
    }
    
    return(sample_posterior)
}

start_time <- Sys.time()
sample_posterior <- NaddGPs_count_model(Y,X,D,beta_0,beta_1,N,params,samplesize)
end_time <- Sys.time()
end_time - start_time

############################################### end NaddGPs count model ###############################


############################################## Start AGP Model #####################################
### prior specification
D <- nrow(Y)
N <- N
C <- nrow(X)
upsilon <- D+1000
beta_0 <- beta_0
f1 <- beta_1
samplesize <- samplesize


AGP_model <- function(Y,X,upsilon,D,
                      sigma1,l1,p,
                      sigma2,l2,
                      beta_0 = beta_0,
                      f1=beta_1,
                      n_samples = samplesize){
    
    s <- matrix(0.9,D-1,D-1)
    diag(s) <- 1
    Xi <- (upsilon-D-1)*(s)
    linear_component_theta <-  cbind(matrix(beta_0,D-1,1), matrix(f1,D-1,1))
    linear_component_gamma <- diag(nrow(X)-1)*0.1 + 1e-6

    # Gamma.time1 <- function(X) RQ(X[1,,drop=F], sigma=sigma1, rho=l1,a=a1)
    Gamma.time1 <- function(X) Period(X,sigma = sigma1, rho = l1, p = p, X_row = 3)
    # Gamma.time2 <- function(X) SE(X[3,,drop=F], sigma=sigma2, rho=l2)
    Gamma.time2 <- function(X) SE(X, sigma=sigma2, rho=l2, X_row = 3)

    # Gamma_burnout <- function(X) SE(X[2,,drop=F], sigma=sigma3, rho=l3)
    Theta_kernel <- function(X) matrix(0, D-1, ncol(X))

    Theta <- list(linear_component_theta,Theta_kernel,Theta_kernel)
    Gamma <- list(linear_component_gamma,Gamma.time1,Gamma.time2)

    # Theta <- list(Theta_kernel,Theta_kernel)
    # Gamma <- list(Gamma.time1,Gamma.time2)

    # Now fit the model
    mod <- fido::basset(Y, X, upsilon, Theta, Gamma, Xi,linear=c(1:2), 
                        verbose = TRUE, n_samples = n_samples,
                        seed = 893,jitter = 1e-4, max_iter = 1000,
                        init = alr_array(Y+0.65, parts = 1))
    return(mod)
}

######################################### Start Bayesian Optimization ############################
MLL_AGP <- function(Y,X,upsilon,D,
                    sigma1,l1,p,
                    sigma2,l2,
                    beta_0,
                    f1,
                    n_samples){
    set.seed(893)
    
    sigma1 <- max(sigma1,sigma2)
    sigma2 <- min(sigma1,sigma2)

    l1 <- min(l1,l2)
    l2 <- max(l1,l2)

    prior_sigma1 <- dnorm(sigma1,mean = 4,sd = 1,log = TRUE)
    prior_sigma2 <- dnorm(sigma2,mean = 1,sd = 1,log = TRUE)

    prior_l1 <- dnorm(l1,mean = 25,sd = 1,log = TRUE)
    prior_l2 <- dnorm(l2,mean = 30,sd = 1,log = TRUE)
    
    prior_log <- prior_sigma1 + prior_sigma2 + prior_l1 + prior_l2
    prior_lambda <- 100

    mod <- AGP_model(Y,X,upsilon,D,
                     sigma1,l1,p,
                     sigma2,l2,
                     beta_0,f1,
                     n_samples = samplesize)
    list(Score=mod$logMarginalLikelihood + prior_log * prior_lambda)
}

set.seed(893)
optObj_AGP <- bayesOpt(
    FUN= function(sigma1,l1,
                  sigma2,l2
                 ){
                  MLL_AGP(Y,X,upsilon,D,
                          sigma1,l1,p = 30,
                          sigma2,l2,
                          beta_0,f1,
                          n_samples = samplesize)},
                  bounds = list( 
                          sigma1 = c(1,hyper_params[1]+5),
                          l1 = c(hyper_params[2]-10,hyper_params[2]+10),
                          sigma2 = c(1,hyper_params[4]+5),
                          l2 = c(hyper_params[5]-10,hyper_params[5]+10)
                          ),
                  initPoints = 10,
                  iters.n = 20,
                  verbose = 1
)

optObj_AGP$scoreSummary
getBestPars(optObj_AGP)
############################################ End Bayesian Optimization ############################
   
##################################### Fit hyperparameters back to the model #######################
# hyper <- hyper_params

hyper <- getBestPars(optObj_AGP)

mod <- AGP_model(Y,X,upsilon,D,
                sigma1 = hyper$sigma1,
                l1 = hyper$l1,
                p = 30,
                sigma2 = hyper$sigma2,
                l2 = hyper$l2,
                beta_0 = beta_0,
                f1=beta_1,
                n_samples = samplesize)

predicted <- predict(mod,response = "Lambda")

#################################################################################


#################### standrization and summarize True \F ##############
true_F <- Summary_taxa2(f1X1+ f2 + f3,iteration = FALSE)
true_Eta <- Summary_taxa2(eta,iteration = FALSE)

f1X1_standarized <- standardize_function(f1X1,iteration = FALSE)
true_f1x1 <- Summary_taxa2(f1X1_standarized,iteration = FALSE)

f2_standarized <- standardize_function(f2,iteration = FALSE)
true_f2 <-Summary_taxa2(f2_standarized,iteration = FALSE)

f3_standarized <- standardize_function(f3,iteration = FALSE)
true_f3 <-Summary_taxa2(f3_standarized,iteration = FALSE)

f1f2_standarized <- standardize_function(f1X1+f3,iteration = FALSE)
true_f1_f2 <- Summary_taxa2(f1f2_standarized,iteration = FALSE)


######################### standrization and summarize MultiAddGPs \F ###################
Sample_F <- sample_posterior[[4]]
Lambda1  <- sample_posterior[[1]]
Lambda2  <- sample_posterior[[2]]
Lambda3  <- sample_posterior[[3]]

sample_tidy <- Summary_taxa2(Sample_F)

## for linear term
f1_sample<- array(dim=c(D-1,N,samplesize))
for (i in 1:dim(mod$Lambda[[1]])[3]){
    f1_sample[,,i] <- mod$Lambda[[1]][,,i] %*% X[1:2,]
}

# for GP term
lambda1 <- standardize_function(f1_sample)
f1_summary <- Summary_taxa2(lambda1)

lambda1_NaddGPs <- standardize_function(Lambda1)
f1_NaddGPs_summary <- Summary_taxa2(lambda1_NaddGPs)

lambda2 <- standardize_function(mod$Lambda[[2]])
f2_summary <- Summary_taxa2(lambda2)

lambda2_NaddGPs <- standardize_function(Lambda2)
f2_NaddGPs_summary <- Summary_taxa2(lambda2_NaddGPs)

lambda3 <- standardize_function(mod$Lambda[[3]])
f3_summary <- Summary_taxa2(lambda3)

lambda3_NaddGPs <- standardize_function(Lambda3)
f3_NaddGPs_summary <- Summary_taxa2(lambda3_NaddGPs)

f12_summary <- Summary_taxa2(lambda1+lambda2)
f12_NaddGPs_summary <- Summary_taxa2(lambda1_NaddGPs+lambda2_NaddGPs)

f123_summary <- Summary_taxa2(lambda1+lambda2+lambda3)
f123_NaddGPs_summary <- Summary_taxa2(lambda1_NaddGPs+lambda2_NaddGPs+lambda3_NaddGPs)
    
# ################################# Predicted F ########################################

# predicted <- predict(mod,response = "Lambda")
predicted_tidy <- Summary_taxa2(predicted)

################################# Plot ########################################
custom_colors <- c('#E43D40','#0D698B','#746C70','#658EA9')

sample_tidy$source <- "NAddGP"
predicted_tidy$source <- "MultiAddGP"
true_F$source <- "True"
combine_NaddGPs_pred <- rbind(sample_tidy,predicted_tidy, true_F)


F_plot <- ggplot(data = combine_NaddGPs_pred, aes(x = Date, y = mean)) +
  # Add ribbons for intervals, ensure differentiation by both batch and source
  geom_ribbon(aes(ymin = p2.5, ymax = p97.5, fill = source), alpha = 0.2) +
  geom_ribbon(aes(ymin = p25, ymax = p75, fill = source), alpha = 0.5) +
  facet_wrap(~batch,ncol=1)+
  # Add mean lines, differentiate by batch with color and source with linetype
  geom_line(aes(color = source,linetype = batch), size = 0.7) +
  # geom_point(data = Y_alr_tidy,aes(x=Date,y=mean), alpha=0.4) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  # scale_color_brewer(palette = "Set1") + # Custom color scheme for batch
  # Additional customizations
  theme_minimal() +
  theme(axis.title.y =element_text(size = 12),
        axis.title.x =element_text(size = 12),
        axis.text.x = element_text(angle = 0, hjust = 1),
        plot.title = element_text(size = 18),
        plot.subtitle = element_text(size =12),
        legend.title = element_blank()) +
  labs(y = "Lambda",x= "Hour")+
  labs(title = "D",
  subtitle = expression(F(X) == beta[0]+f[1]* " \u00D7 " * batch  + f[2](Hourly)+f[3](Daily)))


f1_NaddGPs_summary$source <- "NAddGP_f1"
f1_summary$source <- "MultiAddGP_f1"
true_f1x1$source <- "True_f1"
combine_NaddGPs_pred_f1<- rbind(f1_NaddGPs_summary,f1_summary, true_f1x1)


f1_plot <- ggplot(data = combine_NaddGPs_pred_f1, aes(x = Date, y = mean, group = interaction(batch, source))) +
  # Add ribbons for intervals, ensure differentiation by both batch and source
  geom_ribbon(aes(ymin = p2.5, ymax = p97.5, fill = source, group = interaction(batch, source)), alpha = 0.2) +
  geom_ribbon(aes(ymin = p25, ymax = p75, fill = source, group = interaction(batch, source)), alpha = 0.5) +
  # Add mean lines, differentiate by batch with color and source with linetype
  geom_line(aes(color = source, linetype = batch), size = 1) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  # scale_color_brewer(palette = "Set1") + # Custom color scheme for batch
  # Additional customizations
  theme_minimal() +
  theme(legend.position = "none")+
  theme(axis.title.y =element_text(size = 12),
        axis.title.x =element_text(size = 12),
        axis.text.x = element_text(angle = 0, hjust = 1),
        plot.title = element_text(size = 18),
        plot.subtitle = element_text(size =12)) +
  labs(y = "Lambda",x= "Hour")+
  labs(title = "A",
       subtitle = expression(F(X) == beta[0]+f[1]* " \u00D7 " * Batch))

f2_NaddGPs_summary$source <- "NAddGP_f2"
f2_summary$source <- "MultiAddGP_f2"
true_f3$source <- "True_f2"
combine_NaddGPs_pred_f2<- rbind(f2_NaddGPs_summary,f2_summary, true_f3)[-1,]


    
f2_plot <- ggplot(data = combine_NaddGPs_pred_f2, aes(x = Date, y = mean, group =  source)) +
  # Add ribbons for intervals, ensure differentiation by both batch and source
  geom_ribbon(aes(ymin = p2.5, ymax = p97.5,fill= source, group = source), alpha = 0.2) +
  geom_ribbon(aes(ymin = p25, ymax = p75,fill = source, group = source), alpha = 0.5) +
  # Add mean lines, differentiate by batch with color and source with linetype
  geom_line(aes(color = source), size = 0.7) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  # scale_color_brewer(palette = "Set1") + # Custom color scheme for batch  
  # Additional customizations
  theme_minimal() +
  theme(legend.position = "none")+
  theme(axis.title.y =element_text(size = 12),
        axis.title.x =element_text(size = 12),
        axis.text.x = element_text(angle = 0, hjust = 1),
        plot.title = element_text(size = 18),
        plot.subtitle = element_text(size =12)) +
  labs(y = "Lambda",x= "Hour")+
  labs(title = "B",
       subtitle = expression(F(X) == f[2](Hourly)))



f3_NaddGPs_summary$source <- "NAddGP_f3" 
f3_summary$source <- "MultiAddGP_f3"
true_f2$source <- "True_f3"
combine_NaddGPs_pred_f3<- rbind(f3_NaddGPs_summary,f3_summary, true_f2)[-1,]
    

f3_plot <- ggplot(data = combine_NaddGPs_pred_f3, aes(x = Date, y = mean, group =  source)) +
  # Add ribbons for intervals, ensure differentiation by both batch and source
  geom_ribbon(aes(ymin = p2.5, ymax = p97.5,fill= source,group = source), alpha = 0.2) +
  geom_ribbon(aes(ymin = p25, ymax = p75,   fill= source,group = source), alpha = 0.5) +
  # Add mean lines, differentiate by batch with color and source with linetype
  geom_line(aes(color = source), size = 0.7) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  # scale_color_brewer(palette = "Set1") + # Custom color scheme for batch
  # Additional customizations
  theme_minimal() +
  theme(legend.position = "none")+
  theme(axis.title.y =element_text(size = 12),
        axis.title.x =element_text(size = 12),
        axis.text.x = element_text(angle = 0, hjust = 1),
        plot.title = element_text(size = 18),
        plot.subtitle = element_text(size =12)) +
  labs(y = "Lambda",x= "Hour")+
  labs(title = "C",
       subtitle = expression(F(X) == f[3](Daily)))

layout_matrix <- rbind(c(1, 2, 4),
                       c(3, 3, 4))

p <- grid.arrange(f1_plot,f3_plot,f2_plot,F_plot,layout_matrix=layout_matrix, widths = c(1.1, 1.3, 1.8))
# 



