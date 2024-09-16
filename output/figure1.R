source('/home/thcodelia/Odelia/Project/MultiAddGPs/src/00_package_loading.R', chdir = TRUE)
source('/home/thcodelia/Odelia/Project/MultiAddGPs/src/helper_functions.R', chdir = TRUE)
load("/home/thcodelia/Odelia/Project/MultiAddGPs/data/simulated_data2.RData")
Simulation_result <- readRDS("/home/thcodelia/Odelia/Project/MultiAddGPs/data/Simulation_result2.rds")


Summary_function_sim <- function(X,iteration = TRUE){
  if (iteration == TRUE){
    summary <- gather_array(X,val,coord,sample,iter)%>%
      arrange(iter,coord) %>%
      mutate(time = rep(Lambda_time,samplesize*(D-1))) %>%
      mutate(Date = rep(X1$time,samplesize*(D-1)))%>%
      mutate(batch = rep(X1$batch,samplesize*(D-1)))%>%
      filter(!is.na(val)) %>%
      filter(coord == 2) %>%
      group_by(Date,batch) %>%
      summarise_posterior(val, na.rm=TRUE) %>%
      ungroup()
  }
  else{
    summary <- gather_array(X,val,coord,sample)%>%
      arrange(coord) %>%
      mutate(time = rep(Lambda_time,(D-1))) %>%
      mutate(Date = rep(X1$time,(D-1)))%>%
      mutate(batch = rep(X1$batch,(D-1)))%>%
      filter(!is.na(val)) %>%
      filter(coord == 2) %>%
      group_by(Date,batch) %>%
      summarise_posterior(val, na.rm=TRUE) %>%
      ungroup()
  }
  return(summary)
}

Y <- simulated_data2$Y
X <- simulated_data2$X
eta <- simulated_data2$eta
pai <- simulated_data2$pai
f1X1 <- simulated_data2$f1X1
f2 <- simulated_data2$f2
f3 <- simulated_data2$f3
Sigma <- simulated_data2$Sigma
hyper_params <- simulated_data2$hyper_params
beta_0 <- simulated_data2$beta_0
beta_1 <- simulated_data2$beta_1
X1 <- simulated_data2$X1
D <- simulated_data2$D
N <- simulated_data2$N
samplesize <- simulated_data2$samplesize

sample_posterior <- Simulation_result$sample_posterior
Y_alr_tidy <- Simulation_result$Y_alr_tidy
mod <- Simulation_result$mod
Bayesian_result <- Simulation_result$Bayesian_result
Bayesian_result_sudo <- Simulation_result$Bayesian_result_sudo
predicted <- Simulation_result$predicted
true_F <- Simulation_result$true_F
true_Eta <- Simulation_result$true_Eta
f1_summary <- Simulation_result$f1_summary
f2_summary <- Simulation_result$f2_summary
f3_summary <- Simulation_result$f3_summary
f12_summary <- Simulation_result$f12_summary
f123_summary <- Simulation_result$f123_summary
f123_summary <- Simulation_result$f123_summary


Lambda_time <- X1$time

#################### standrization and summarize True \F ##############
true_F <- Summary_function_sim(f1X1+ f2 + f3,iteration = FALSE)
true_Eta <- Summary_function_sim(eta,iteration = FALSE)

f1X1_standarized <- standardize_function(f1X1,iteration = FALSE)
true_f1x1 <- Summary_function_sim(f1X1_standarized,iteration = FALSE)

f2_standarized <- standardize_function(f2,iteration = FALSE)
true_f2 <-Summary_function_sim(f2_standarized,iteration = FALSE)

f3_standarized <- standardize_function(f3,iteration = FALSE)
true_f3 <-Summary_function_sim(f3_standarized,iteration = FALSE)

f1f2_standarized <- standardize_function(f1X1+f3,iteration = FALSE)
true_f1_f2 <- Summary_function_sim(f1f2_standarized,iteration = FALSE)


######################### standrization and summarize MultiAddGPs \F ###################
Sample_F <- sample_posterior[[4]]
Lambda1  <- sample_posterior[[1]]
Lambda2  <- sample_posterior[[2]]
Lambda3  <- sample_posterior[[3]]

sample_tidy <- Summary_function_sim(Sample_F)

## for linear term
f1_sample<- array(dim=c(D-1,N,samplesize))
for (i in 1:dim(mod$Lambda[[1]])[3]){
    f1_sample[,,i] <- mod$Lambda[[1]][,,i] %*% X[1:2,]
}

# for GP term
lambda1 <- standardize_function(f1_sample)
f1_summary <- Summary_function_sim(lambda1)

lambda1_sudo <- standardize_function(Lambda1)
f1_sudo_summary <- Summary_function_sim(lambda1_sudo)

lambda2 <- standardize_function(mod$Lambda[[2]])
f2_summary <- Summary_function_sim(lambda2)

lambda2_sudo <- standardize_function(Lambda2)
f2_sudo_summary <- Summary_function_sim(lambda2_sudo)

lambda3 <- standardize_function(mod$Lambda[[3]])
f3_summary <- Summary_function_sim(lambda3)

lambda3_sudo <- standardize_function(Lambda3)
f3_sudo_summary <- Summary_function_sim(lambda3_sudo)

f12_summary <- Summary_function_sim(lambda1+lambda2)
f12_sudo_summary <- Summary_function_sim(lambda1_sudo+lambda2_sudo)

f123_summary <- Summary_function_sim(lambda1+lambda2+lambda3)
f123_sudo_summary <- Summary_function_sim(lambda1_sudo+lambda2_sudo+lambda3_sudo)
    
# ################################# Predicted F ########################################

# predicted <- predict(mod,response = "Lambda")
predicted_tidy <- Summary_function_sim(predicted)

################################# Plot ########################################
custom_colors <- c('#E43D40','#0D698B','#746C70','#658EA9')

sample_tidy$source <- "NAddGP"
predicted_tidy$source <- "MultiAddGP"
true_F$source <- "True"
combine_sudo_pred <- rbind(sample_tidy,predicted_tidy, true_F)


F_plot <- ggplot(data = combine_sudo_pred, aes(x = Date, y = mean)) +
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


f1_sudo_summary$source <- "NAddGP_f1"
f1_summary$source <- "MultiAddGP_f1"
true_f1x1$source <- "True_f1"
combine_sudo_pred_f1<- rbind(f1_sudo_summary,f1_summary, true_f1x1)


f1_plot <- ggplot(data = combine_sudo_pred_f1, aes(x = Date, y = mean, group = interaction(batch, source))) +
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

f2_sudo_summary$source <- "NAddGP_f2"
f2_summary$source <- "MultiAddGP_f2"
true_f3$source <- "True_f2"
combine_sudo_pred_f2<- rbind(f2_sudo_summary,f2_summary, true_f3)[-1,]


    
f2_plot <- ggplot(data = combine_sudo_pred_f2, aes(x = Date, y = mean, group =  source)) +
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



f3_sudo_summary$source <- "NAddGP_f3" 
f3_summary$source <- "MultiAddGP_f3"
true_f2$source <- "True_f3"
combine_sudo_pred_f3<- rbind(f3_sudo_summary,f3_summary, true_f2)[-1,]
    

f3_plot <- ggplot(data = combine_sudo_pred_f3, aes(x = Date, y = mean, group =  source)) +
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






