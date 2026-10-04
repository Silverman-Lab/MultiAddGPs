## Artificial gut analysis: four vessel baselines and two post-starvation effects

library(fido)
library(dplyr)
library(tidyr)
library(Biostrings)
library(phyloseq)
library(ggplot2)
library(purrr)
library(gridExtra)
library(grid)
library(gtable)
library(patchwork)
library(ggtext)

set.seed(59292)



# Data preparation --------------------------------------------------------
analysis_dir <- "."
figure_dir <- file.path(analysis_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(analysis_dir, "Artificial_gut_analysis_helpers.R"), chdir = TRUE)
load(file.path(analysis_dir, "data", "Data2_Mallard", "mallard_family.RData"))

mallard_family$sample_data <- mallard_family$sample_data %>% 
    select(X.SampleID,time,Vessel,SampleType,batch)%>%
    group_by(Vessel)%>%
    arrange(time)%>%
    filter(!duplicated(time))%>%
    mutate(time_num = as.numeric(time),
           Hour = (time_num - min(time_num))/3600,
           SampleType = ifelse(time_num == "1448031600","Hourly",SampleType),
           SampleType = ifelse(time_num %in% c("1448550000","1448636400","1448722800","1448809200"),"Daily",SampleType),
           Hour_diff = (Hour-mean(Hour))/sd(Hour), 
           Hour_sample = ifelse(SampleType == "Hourly",1,0),
           Daily_sample = ifelse(SampleType == "Daily",1,0),
           burnout  = ifelse(time < "2015-11-11 15:00:00",1,0),
           Vessel1 = ifelse(Vessel == 1,1,0),
           Vessel2 = ifelse(Vessel == 2,1,0),
           Vessel3 = ifelse(Vessel == 3,1,0),
           Vessel4 = ifelse(Vessel == 4,1,0),
           Vessel1Burn = ifelse(Vessel == 1 & burnout != 1,1,0),
           Vessel2Burn = ifelse(Vessel == 2 & burnout != 1,1,0)
           )%>%
    arrange(Vessel)


#match sample_data to the otu_table
mallard_family$otu_table <- mallard_family$otu_table[mallard_family$sample_data$X.SampleID,] 

# Design matrix for vessel baselines, post-starvation effects, and time
X1 <- mallard_family$sample_data %>%
               dplyr::select(X.SampleID,Vessel,Vessel1,Vessel2,Vessel3,Vessel4,Vessel1Burn,Vessel2Burn,Hour_diff,Hour)

Hour <- X1$Hour
Hour_diff <- X1$Hour_diff
Vessel <- X1$Vessel
Vessel1 <- X1$Vessel1
Vessel2 <- X1$Vessel2
Vessel3 <- X1$Vessel3
Vessel4 <- X1$Vessel4
Vessel1Burn <- X1$Vessel1Burn
Vessel2Burn <- X1$Vessel2Burn

X<- t(model.matrix(~ Vessel1+Vessel2+Vessel3+Vessel4+Vessel1Burn+Vessel2Burn+Hour_diff+0,data = X1)) # Transpose to get the right dimensions

#  Extract Data / dimensions from Phyloseq object
Y <- t(as(mallard_family$otu_table, "matrix"))
rownames(Y) <- mallard_family$tax_table$Family

if (identical(colnames(Y), X1$X.SampleID)) {
  colnames(Y) <- colnames(X)
} else {
  warning("The order of the samples in the otu_table and sample_data are not the same")
}
D <- nrow(Y)
N <- nrow(mallard_family$sample_data)
C <- nrow(X)

# Model specification -----------------------------------------------------
upsilon <- D+5
Model <- function(Y,X,upsilon,
                  sigma1,l1,
                  sigma2,l2,a2,
                  noise
                  ) {
    set.seed(893)
    Xi <- matrix(0.4,D-1,D-1)
    diag(Xi) <- noise 
  # stationary kernel implementation
    Gamma_vessel1 <-function(X) {SE(X,sigma = sigma1, rho = l1, specific_value = 1,Identity_row = 1,X_row = C)}
    Gamma_vessel2 <-function(X) {SE(X,sigma = sigma1, rho = l1, specific_value = 1,Identity_row = 2,X_row = C)}
    Gamma_vessel3 <-function(X) {SE(X,sigma = sigma1, rho = l1, specific_value = 1,Identity_row = 3,X_row = C)}
    Gamma_vessel4 <-function(X) {SE(X,sigma = sigma1, rho = l1, specific_value = 1,Identity_row = 4,X_row = C)}
    # Rows 5 and 6 are post-starvation indicators for Vessels 1 and 2.
    Gamma_vessel1Disrupt <-function(X) {RQ(X,sigma = sigma2, rho = l2,a=a2, specific_value = 1,Identity_row = 5,X_row = C)}
    Gamma_vessel2Disrupt <-function(X) {RQ(X,sigma = sigma2, rho = l2,a=a2, specific_value = 1,Identity_row = 6,X_row = C)}
    Theta_kernel <- function(X) matrix(0, D-1, ncol(X))
    Theta <- list(Theta_kernel,Theta_kernel,Theta_kernel,Theta_kernel,Theta_kernel,Theta_kernel)
    Gamma <- list(Gamma_vessel1, Gamma_vessel2,Gamma_vessel3,Gamma_vessel4,Gamma_vessel1Disrupt,Gamma_vessel2Disrupt)
    mod <- fido::basset(Y, X, upsilon, Theta, Gamma, Xi, verbose = TRUE, seed = 893)
}


# Hyperparameters ---------------------------------------------------------
# Set to TRUE to re-run Bayesian hyperparameter optimization.
use_optimization <- FALSE

if (use_optimization) {
  MLL <- function(Y, X, upsilon, sigma1, l1, sigma2, l2, a2, noise) {
    set.seed(893)
    sigma <- sort(c(sigma1, sigma2)) # smaller magnitude: baseline
    rho <- sort(c(l1, l2))           # larger length scale: baseline
    alpha1 <- 10
    beta1 <- 16
    alpha2 <- 20
    beta2 <- 16
    log_prior_l1 <- alpha1 * log(beta1) - lgamma(alpha1) - (alpha1 + 1) * log(rho[2]) - beta1 / rho[2]
    log_prior_l2 <- alpha2 * log(beta2) - lgamma(alpha2) - (alpha2 + 1) * log(rho[1]) - beta2 / rho[1]
    log_prior_s1 <- alpha2 * log(beta2) - lgamma(alpha2) - (alpha2 + 1) * log(sigma[1]) - beta2 / sigma[1]
    log_prior_s2 <- alpha1 * log(beta1) - lgamma(alpha1) - (alpha1 + 1) * log(sigma[2]) - beta1 / sigma[2]
    prior <- log_prior_l1 + log_prior_l2 + log_prior_s1 + log_prior_s2
    prior_lambda <- 120
    mod <- Model(Y, X, upsilon, sigma[1], rho[2], sigma[2], rho[1], a2, noise)
    list(Score = mod$logMarginalLikelihood + prior * prior_lambda)
  }

  set.seed(893)
  optObj <- ParBayesianOptimization::bayesOpt(
    FUN = function(sigma1, l1, sigma2, l2)
      MLL(Y, X, upsilon, sigma1, l1, sigma2, l2, a2 = 2, noise = 1),
    bounds = list(sigma1 = c(0.1, 2), l1 = c(1, 3),
                  sigma2 = c(1, 4), l2 = c(0.1, 2)),
    initPoints = 10,
    iters.n = 20,
    verbose = 1
  )
  optObj$scoreSummary
  ParBayesianOptimization::getBestPars(optObj)
  # Diagnostic view of the optimization scores from the original analysis.
  plotly::plot_ly(
    x = ~optObj$scoreSummary$sigma2,
    y = ~optObj$scoreSummary$l2,
    z = ~optObj$scoreSummary$Score,
    type = "scatter3d", mode = "markers",
    marker = list(size = 5, colorscale = c("Blues"))
  )
  best <- ParBayesianOptimization::getBestPars(optObj)
  hyper <- c(min(best$sigma1, best$sigma2), max(best$l1, best$l2),
             max(best$sigma1, best$sigma2), min(best$l1, best$l2), 2, 1)
} else {
  hyper <- c(1.03183062716154, 1.62553177010268,
             1.61533599942923, 0.849494453410152, 2, 1) # optimized hyperparameters
}
mod <- Model(Y,X,upsilon,sigma1=hyper[1],l1=hyper[2],
                         sigma2=hyper[3],l2=hyper[4], 
                         a2 = hyper[5],
                         noise = hyper[6])



# Posterior summaries in CLR coordinates -------------------------------
mod_clr <- to_clr(mod)
predicted <- predict(mod_clr,response = "Lambda")
rownames(predicted) <- rownames(Y)
taxa_names <- as(mallard_family$tax_table$Family,"vector")


# summarize prediction in clr space 
predicted_clr_tidy <- gather_array(predicted,val,coord,sample,iter)%>%
    arrange(coord,iter)%>%
    mutate(coord = paste0("CLR(",taxa_names[coord], ")"))%>%
    mutate(coord = factor(coord,levels = unique(coord)))%>%
    mutate(Hour = rep(Hour,max(iter)*D))%>%
    mutate(Vessel = as.factor(rep(Vessel,max(iter)*D)))%>%
    filter(!is.na(val))%>%
    group_by(Vessel,coord,Hour)%>%
    summarise_posterior(val,na.rm = TRUE)%>%
    ungroup()
    # filter(coord == 5)

Y_clr_tidy <- clr_array(Y+0.65, parts = 1) %>% 
  gather_array(mean, coord, sample) %>% 
  mutate(coord = paste0("CLR(",taxa_names[coord], ")"))%>%
  mutate(coord = factor(coord,levels = unique(coord)))%>%  
  arrange(coord) %>%
  mutate(Hour = rep(Hour,D))%>%
  mutate(Vessel = as.factor(rep(Vessel,D)))%>%
  mutate(source = "Y")
  # filter(coord == 5)

custom_colors <- c('#5ba300','#f57600','#8babf1','#0073e6')
custom_colors_v12 <- c('#5ba300','#f57600')

# Calculate ymin and ymax for each facet
facet_lims <- Y_clr_tidy %>%
  group_by(coord) %>%
  summarise(ymin = min(mean), ymax = max(mean))

# Merge the calculated limits back to the original data
predicted_clr_all <- predicted_clr_tidy %>%
  left_join(facet_lims, by = "coord")
predicted_clr_tidy <- predicted_clr_all %>%
  filter(coord == "CLR(Rikenellaceae)" | coord == "CLR(Fusobacteriaceae)")



# Starvation-component summaries ---------------------------------------
combine_vessel_summaries <- function(mod_clr, lambda_component,filter_lambda_vessel1, filter_lambda_vessel2, time_column, coords, CLR = TRUE) {
  combined_summary <- data.frame() 
  
  component1<-lambda_component[1]
  component2<-lambda_component[2]

  for (coord in coords) {
    vessel1_summary <- Summary_function(
      mod_clr$Lambda[[component1]],
      iteration = TRUE,
      filter_lambda = filter_lambda_vessel1,
      time_column = time_column,
      filter_coord = coord,
      CLR = CLR
    )
    vessel2_summary <- Summary_function(
      mod_clr$Lambda[[component2]],
      iteration = TRUE,
      filter_lambda = filter_lambda_vessel2,
      time_column = time_column,
      filter_coord = coord,
      CLR = CLR
    )
    # Add the Vessel and Coord identifiers
    vessel1_summary$Vessel <- "1"
    vessel2_summary$Vessel <- "2"
    vessel1_summary$Coord <- coord
    vessel2_summary$Coord <- coord
    # Combine the summaries
    combined_summary <- rbind(combined_summary, vessel1_summary, vessel2_summary)
  }
  return(combined_summary)
}


dis_allcoord <- combine_vessel_summaries(mod_clr,lambda_component = c(5,6), Vessel1, Vessel2, Hour, c(1:10))%>%
                arrange(Coord)%>%
                mutate(Coord = paste0("CLR(",taxa_names[Coord], ")"))%>%
                mutate(Coord = factor(Coord,levels = unique(Coord)))

dis_allcoord_all <- dis_allcoord
dis_allcoord <- dis_allcoord_all %>%
                filter(Coord == "CLR(Rikenellaceae)" | Coord == "CLR(Fusobacteriaceae)")

p <- ggplot(predicted_clr_tidy,aes(x=Hour,y=mean))+
      facet_wrap(~coord, scales="free_y", ncol=1) +
      geom_ribbon(aes(ymin = p2.5, ymax = p97.5,fill= Vessel), alpha = 0.2) +
      geom_ribbon(aes(ymin = p25,  ymax = p75,  fill= Vessel), alpha = 0.5) +
      geom_line(aes(color=Vessel),alpha = 0.7,size = 0.7) +
      # geom_point(data = Y_clr_tidy, aes(x = Hour, y = mean,color=Vessel),alpha =0.3)+
      scale_fill_manual(values = custom_colors) +        
      scale_color_manual(values = custom_colors) +
      # Additional customizations
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(color = "black", size =20,face = "italic"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.y = element_text(size = 20, margin = margin(t = 15)),
        axis.title.x = element_text(size = 20, margin = margin(t = 15)),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 20),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 20),
        plot.title = element_text(size = 25, margin = margin(b = 20), hjust = 0.5),
        legend.spacing.y = unit(16, "pt"),
        legend.key.height = unit(1.5, "cm"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))+
        # panel.spacing = unit(2, "lines")) +
      labs(y = "F", x= "Hour", title = expression(F == f^{" longterm"}+f^{" starvation"}))+
      annotate("rect", xmin = 240, xmax = 264, ymin = -Inf, ymax = Inf, fill = "#717D7E", alpha = 0.25) 


disruption_plot <- ggplot(data = dis_allcoord, aes(x = time, y = mean)) +
    facet_wrap(~Coord, scales = "free_y", ncol = 1) +
    geom_ribbon(aes(ymin = p2.5, ymax = p97.5, fill = Vessel), alpha = 0.2) +
    geom_ribbon(aes(ymin = p25, ymax = p75, fill = Vessel), alpha = 0.5) +
    geom_line(aes(color = Vessel), alpha = 0.4, size = 0.7) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    scale_fill_manual(values = custom_colors_v12) +
    scale_color_manual(values = custom_colors_v12) +
    annotate("rect", xmin = 240, xmax = 264, ymin = -Inf, ymax = Inf,
            fill = "#717D7E", alpha = 0.25) +
    labs(y = "", x = "Hour", title = expression(f^{" starvation"})) +
    guides(fill = "none", color = "none") +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(color = "black", size = 20, face = "italic"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.title.y = element_text(size = 20, margin = margin(t = 15)),
      axis.title.x = element_text(size = 20, margin = margin(t = 15)),
      axis.text.x = element_text(angle = 0, hjust = 1, size = 20),
      axis.text.y = element_text(angle = 0, hjust = 1, size = 20),
      plot.title = element_text(size = 25, margin = margin(b = 20), hjust = 0.5),
      legend.position = "none"
    )

two_family_plot <- p + disruption_plot
two_family_plot <- two_family_plot + plot_layout(guides = "collect") &
  theme(legend.position = "right")
ggsave(file.path(figure_dir, "TwoFamilyCombined.png"), plot = two_family_plot, width = 18, height = 12)

# All-family figures from the 2024 analysis ------------------------------
p_all <- ggplot(predicted_clr_all, aes(x = Hour, y = mean)) +
  facet_wrap(~coord, scales = "free_y", ncol = 5) +
  geom_ribbon(aes(ymin = p2.5, ymax = p97.5, fill = Vessel), alpha = 0.4) +
  geom_line(aes(color = Vessel), alpha = 0.7, linewidth = 0.7) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(color = "black", size = 30),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title.y = element_text(size = 35, margin = margin(t = 15)),
    axis.title.x = element_text(size = 35, margin = margin(t = 15)),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    plot.title = element_text(size = 27, margin = margin(b = 20)),
    legend.spacing.y = unit(16, "pt"),
    legend.key.height = unit(1.5, "cm"),
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 30),
    panel.spacing = unit(2, "lines")
  ) +
  labs(y = "F", x = "Hour") +
  annotate("rect", xmin = 240, xmax = 264, ymin = -Inf, ymax = Inf,
           fill = "#717D7E", alpha = 0.25)

ggsave(file.path(figure_dir, "AllVessels.png"), plot = p_all, width = 35, height = 18)

disruption_plot_all <- ggplot(dis_allcoord_all, aes(x = time, y = mean)) +
  facet_wrap(~Coord, scales = "free_y", ncol = 5) +
  geom_ribbon(aes(ymin = p2.5, ymax = p97.5, fill = Vessel), alpha = 0.2) +
  geom_ribbon(aes(ymin = p25, ymax = p75, fill = Vessel), alpha = 0.5) +
  geom_line(aes(color = Vessel), alpha = 0.4, linewidth = 0.7) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  scale_fill_manual(values = custom_colors_v12) +
  scale_color_manual(values = custom_colors_v12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(color = "black", size = 30),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title.y = element_text(size = 35, margin = margin(t = 15)),
    axis.title.x = element_text(size = 35, margin = margin(t = 15)),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    plot.title = element_text(size = 27, margin = margin(b = 20)),
    legend.spacing.y = unit(16, "pt"),
    legend.key.height = unit(1.5, "cm"),
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 30),
    panel.spacing = unit(2, "lines")
  ) +
  labs(y = "F", x = "Hour") +
  annotate("rect", xmin = 240, xmax = 264, ymin = -Inf, ymax = Inf,
           fill = "#717D7E", alpha = 0.25)

ggsave(file.path(figure_dir, "Disruption.png"), plot = disruption_plot_all, width = 35, height = 18)


# Inverse-gamma priors used below: baseline amplitude/disruption length scale
# (20, 16), and baseline length scale/disruption amplitude (10, 16).
dir.create(figure_dir, showWarnings = FALSE)
png(file.path(figure_dir, "Web_Figure6.png"), width = 2100, height = 1200, res = 300)
par(mfrow = c(1, 2))
prior_x <- seq(0.001, 4, length.out = 1000)
for (prior_shape in c(20, 10)) {
  plot(prior_x, dgamma(1 / prior_x, shape = prior_shape, rate = 16) / prior_x^2,
       type = "l", col = "darkred", lwd = 1.5, xlim = c(0, 4), yaxt = "n",
       xlab = sprintf("InverseGamma(%d, 16)", prior_shape),
       ylab = if (prior_shape == 20) "Prior Density" else "")
}
dev.off()
