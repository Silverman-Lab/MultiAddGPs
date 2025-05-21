
source("/home/thcodelia/Odelia/Project/MultiAddGPs/src/helper_functions.R")
source("/home/thcodelia/Odelia/Project/MultiAddGPs/src/package_loading.R")
source("/home/thcodelia/Odelia/Project/MultiAddGPs/src/simulation_function.R")

library(lgpr)

custom_colors <- c('#E43D40','#0D698B','#746C70','#658EA9')



simulation_warp <- Simulation_warp(D=3,N=50,hyper_params=c(2,3,1,0.5),multi.depth=2000,seed=5332)

X <- simulation_warp$X
Y <- simulation_warp$Y
D <- nrow(Y)
N <- ncol(X)
C <- nrow(X)
time <- seq(0,N-1,by = 1)


################################################## fit lgpr model ###################################################
X1 <- as.data.frame(t(X))
data <- cbind(X1, Y=Y[2,]) %>% 
        mutate(Batch = as.factor(ifelse(batchBatch2 == 1, "Batch1", "Batch2")),
               ID = as.factor("001"),
               time = seq(0,N-1,by = 1))%>%
        dplyr::select(-batchBatch2)
lgpr <- readRDS("~/Simulation/lgpr.rds")
p <- pred(lgpr,data,draws = sample.int(5000,2000))

##################################################### fit MultiAddGP model #####################################################



hyper_params <- simulation_warp$hyper_params
depth <- simulation_warp$depth


upsilon <- D+100

Gamma_general <- function(X){SE(X,sigma = hyper_params[3], rho = hyper_params[4], X_row=3,jitter = 1e-6)}
Gamma_treatment <- function(X){SE(X,sigma =hyper_params[1], rho = hyper_params[2],  X_row=4,jitter = 1e-6)}
Theta_kernel <- function(X){matrix(0,D-1,ncol(X))}

linear_component_theta <- matrix(0,D-1,2) 
linear_component_gamma <- diag(2) + 1e-6 # a matrix with (row and col) dim of the length of linear components

# Remember the order of putting the linear component is important, always put the linear component first
Theta <- list(linear_component_theta,Theta_kernel,Theta_kernel)
Gamma <- list(linear_component_gamma,Gamma_general,Gamma_treatment)
# Note: the order of output lambda depends on the order of the input lambda list

fit <- fido::basset(Y,X,upsilon,Theta=Theta,Gamma=Gamma,linear =c(1:2),
                    samplesize = 2000,verbose = TRUE,seed = 893,
                    jitter = 1e-4, max_iter = 1000,
                    init = alr_array(Y+0.65, parts = 1))


str(fit)
fit <- to_clr(fit)
predictX <- X
predictY <- Y
Predict <- predict(fit,newdata = X, reponse="Lambda")

########################################### fit back to F ###########################################
Predict_tidy <- gather_array(Predict,val,coord,sample,iter)%>%
                mutate(Day = time[sample],
                      Batch = as.factor(X[2,sample]),
                      Batch = ifelse(Batch == 1, "Batch1","Batch2")) %>%
                filter(!is.na(val))%>%
                group_by(Day,coord,Batch)%>%
                summarise_posterior(val,na.rm = TRUE)%>%
                ungroup()%>%
                filter(coord == "2")%>%
                dplyr::select(-coord)


F <- alr_to_clr(simulation_warp$F,D)
True_F <-  standardize_function(F,iteration = FALSE)

True_F_tidy <- gather_array(F, mean, coord, sample) %>% 
               mutate(Day = time[sample], 
                      Batch = as.factor(X[2,sample]),
                      Batch = ifelse(Batch == 1, "Batch1","Batch2")) %>%
               filter(!is.na(mean))%>%
               group_by(Day,coord,Batch)%>%
               summarise_posterior(mean,na.rm = TRUE)%>%
               ungroup() %>%
               filter(coord == "2")%>%
               dplyr::select(-coord)
        
F_tidy <- p@f %>%
       t() %>%
       data.frame() %>%
       rownames_to_column("Sample") %>%
       pivot_longer(cols = -Sample, names_to = "Iteration", values_to = "val")%>%
       mutate(Iteration = as.numeric(gsub("X","",Iteration)),
              Sample = as.numeric(Sample))%>%
       mutate(Day = time[Sample],
              Batch = data$Batch[Sample])%>%
       group_by(Day,Batch)%>%
       summarise_posterior(val,na.rm = TRUE)%>%
       ungroup()


True_F_tidy$source <- "True"
Predict_tidy$source <- "MultiAddGPs"
F_tidy$source <- "Lgpr"

combined_tidy_F <- rbind(True_F_tidy, Predict_tidy,F_tidy)

F_plot <- ggplot(combined_tidy_F, aes(x = Day, y = mean)) +  
      facet_wrap(~Batch,ncol=1)+
        # facet_wrap(~coord, scales="free_y",ncol=4)+
        geom_line(aes(color = source,linetype = Batch),size = 0.7)+
        geom_ribbon(aes(ymin = p2.5, ymax=p97.5,fill=source),alpha=0.3)+
        geom_ribbon(aes(ymin =p25,ymax=p75,fill=source),alpha=0.8)+
        theme_minimal() +
        scale_fill_manual(values = custom_colors) +
        scale_color_manual(values = custom_colors) +
        # theme(legend.position = "none")+
        theme(axis.title.y =element_text(size = 20,margin = margin(t = 15)),
              axis.title.x =element_text(size = 20,margin = margin(t = 15)),
              axis.text.x = element_text(angle = 0, hjust = 1,size =17),
              axis.text.y = element_text(angle = 0, hjust = 1,size =17),
              plot.title = element_text(size = 30,margin = margin(b = 20)),
              plot.subtitle = element_text(size =27,margin = margin(b = 25)),
              legend.title = element_blank(),
              strip.text = element_text(size = 20),
              legend.text = element_text(size = 17),
              legend.spacing.y = unit(16, "pt"),
              legend.key.height = unit(1.5, "cm")) +
        labs(y = "F",x= "Time")+
        labs(subtitle = expression(F(X)==b[0]+b[1]*x^(Batch)+f^(Non-Stationary)+f^(Trend)))


####################################### fit each Lambda ########################################


Lambda_linear <- array(0,dim = dim(fit$Lambda[2][[1]]))
Lambda_General <-  array(0,dim = dim(fit$Lambda[3][[1]]))
Lambda_Treatment <- array(0,dim = dim(fit$Lambda[2][[1]]))

## seperate the Lambda

for (i in 1:dim(fit$Lambda[2][[1]])[3]){
    Lambda_linear[,,i] <- fit$Lambda[1][[1]][,,i]%*%X[1:2,]}

for(i in 1:dim(fit$Lambda[2][[1]])[3]){
     Lambda_General[,,i] <- fit$Lambda[2][[1]][,,i]}

for(i in 1:dim(fit$Lambda[2][[1]])[3]){
     Lambda_Treatment[,,i] <- fit$Lambda[3][[1]][,,i] }


# Center the Lambda and Y
Lambda_linear <- standardize_function(Lambda_linear,iteration = TRUE)
Lambda_General<-standardize_function(Lambda_General,iteration = TRUE)
Lambda_Treatment<-standardize_function(Lambda_Treatment,iteration = TRUE)
Y_st <- standardize_function(clr_array(Y+0.65, parts = 1),iteration = FALSE)


Lambda_linear_tidy <- gather_array(Lambda_linear,val,coord,sample,iter)%>%
                    mutate(Day = time[sample],
                           Batch = as.factor(X[2,sample]),
                           Batch = ifelse(Batch == 1, "Batch1","Batch2")) %>%
                    filter(!is.na(val))%>%
                    group_by(Day,coord,Batch)%>%
                    summarise_posterior(val,na.rm = TRUE)%>%
                    ungroup() %>%
                    filter(coord == "2")%>%
                    dplyr::select(-coord)


Lambda_General_tidy <- gather_array(Lambda_General,val,coord,sample,iter)%>%
                    mutate(Day = time[sample],
                           Batch = as.factor(X[2,sample]),
                           Batch = ifelse(Batch == 1, "Batch1","Batch2")) %>%
                    filter(!is.na(val))%>%
                    group_by(Day,coord,Batch)%>%
                    summarise_posterior(val,na.rm = TRUE)%>%
                    ungroup()%>%
                    filter(coord == "2")%>%
                    dplyr::select(-coord)

Lambda_Treatment_tidy <- gather_array(Lambda_Treatment,val,coord,sample,iter)%>%
                    mutate(Day = time[sample],
                            Batch = as.factor(X[2,sample]),
                            Batch = ifelse(Batch == 1, "Batch1","Batch2")) %>%
                    filter(!is.na(val))%>%
                    group_by(Day,coord,Batch)%>%
                    summarise_posterior(val,na.rm = TRUE)%>%
                    ungroup()%>%
                    filter(coord == "2")%>%
                    dplyr::select(-coord)    


####################################### true lambda ########################################

f1 <- alr_to_clr(simulation_warp$f1X1,D)
f2 <- alr_to_clr(simulation_warp$f2,D)
f3 <- alr_to_clr(simulation_warp$f3,D)


True_f1X1 <- standardize_function(f1,iteration = FALSE)
True_f2 <- standardize_function(f2,iteration = FALSE)
True_f3 <- standardize_function(f3,iteration = FALSE)

# True_F_tidy <- gather_array(True_F, mean, coord, sample) %>% 
#                mutate(Day = time[sample], 
#                      Batch = as.factor(X[2,sample]),
#                      Batch = ifelse(Batch == 1, "Batch1","Batch2")) %>%
#                filter(coord == "2")%>%
#                dplyr::select(-coord)

True_f1X1_tidy <- gather_array(True_f1X1, mean, coord, sample) %>%
                    mutate(Day = time[sample],
                            Batch = as.factor(X[2,sample]),
                            Batch = ifelse(Batch == 1, "Batch1","Batch2")) %>%
                    filter(!is.na(mean))%>%
                    group_by(Day,coord,Batch)%>%
                    summarise_posterior(mean,na.rm = TRUE)%>%
                    ungroup() %>%
                    filter(coord == "2")%>%
                    dplyr::select(-coord)    


True_f2_tidy <- gather_array(True_f2, mean, coord, sample) %>%
                    mutate(Day = time[sample],
                            Batch = as.factor(X[2,sample]),
                            Batch = ifelse(Batch == 1, "Batch1","Batch2")) %>%
                    filter(!is.na(mean))%>%
                    group_by(Day,coord,Batch)%>%
                    summarise_posterior(mean,na.rm = TRUE)%>%
                    ungroup() %>%
                    filter(coord == "2")%>%
                    dplyr::select(-coord)    


True_f3_tidy <- gather_array(True_f3, mean, coord, sample) %>%
                    mutate(Day = time[sample],
                            Batch = as.factor(X[2,sample]),
                            Batch = ifelse(Batch == 1, "Batch1","Batch2")) %>%
                    filter(!is.na(mean))%>%
                    group_by(Day,coord,Batch)%>%
                    summarise_posterior(mean,na.rm = TRUE)%>%
                    ungroup() %>%
                    filter(coord == "2")%>%
                    dplyr::select(-coord)    

##################################### lgpr lambda ########################################


Linear_tidy_lgpr <- p@f_comp[3][[1]] %>%
                    t() %>%
                    data.frame() %>%
                    rownames_to_column("Sample") %>%
                    pivot_longer(cols = -Sample, names_to = "Iteration", values_to = "val")%>%
                    mutate(Iteration = as.numeric(gsub("X","",Iteration)),
                            Sample = as.numeric(Sample))%>%
                    mutate(Day = time[Sample],
                            Batch = data$Batch[Sample])%>%
                    group_by(Day,Batch)%>%
                    summarise_posterior(val,na.rm = TRUE)%>%
                    ungroup()

General_tidy_lgpr <- p@f_comp[2][[1]] %>%
                    t() %>%
                    data.frame() %>%
                    rownames_to_column("Sample") %>%
                    pivot_longer(cols = -Sample, names_to = "Iteration", values_to = "val")%>%
                    mutate(Iteration = as.numeric(gsub("X","",Iteration)),
                            Sample = as.numeric(Sample))%>%
                    mutate(Day = time[Sample],
                            Batch = data$Batch[Sample])%>%
                    group_by(Day,Batch)%>%
                    summarise_posterior(val,na.rm = TRUE)%>%
                    ungroup()

Treatment_tidy_lgpr <- p@f_comp[4][[1]] %>%
                    t() %>%
                    data.frame() %>%
                    rownames_to_column("Sample") %>%
                    pivot_longer(cols = -Sample, names_to = "Iteration", values_to = "val")%>%
                    mutate(Iteration = as.numeric(gsub("X","",Iteration)),
                            Sample = as.numeric(Sample))%>%
                    mutate(Day = time[Sample],
                            Batch = data$Batch[Sample])%>%
                    group_by(Day,Batch)%>%
                    summarise_posterior(val,na.rm = TRUE)%>%
                    ungroup()

########################################### plot ###########################################

True_f1X1_tidy$source <- "True"
Lambda_linear_tidy$source <- "Predicted"
Linear_tidy_lgpr$source <- "Lgpr"

combined_tidy_f1X1 <- rbind(True_f1X1_tidy, Lambda_linear_tidy,Linear_tidy_lgpr)


f1_plot <- ggplot(combined_tidy_f1X1, aes(x = Day, y = mean,group = interaction(Batch, source))) +  
            # facet_wrap(~coord, scales="free_y",ncol=4)+
            geom_line(aes(color = source, linetype = Batch),size = 0.7)+
            # geom_line() +
        #     geom_ribbon(aes(ymin = p2.5, ymax=p97.5,fill = source, group = interaction(Batch, source)),alpha=0.2)+
            geom_ribbon(aes(ymin =p25,ymax=p75,fill = source, group = interaction(Batch, source)),alpha=0.5)+
            scale_fill_manual(values = custom_colors) +        
            scale_color_manual(values = custom_colors)+
            theme_minimal() +
            theme(legend.position = "none")+
            theme(axis.title.y =element_text(size = 20,margin = margin(t = 15)),
                  axis.title.x =element_text(size = 20,margin = margin(t = 15)),
                  axis.text.x = element_text(angle = 0, hjust = 1,size =17),
                  axis.text.y = element_text(angle = 0, hjust = 1,size =17),
                  plot.title = element_text(size = 27,margin = margin(b = 20)),
                  plot.subtitle = element_text(size =20,margin = margin(b = 25))) +
            labs(y = "F",x= "Hour")+
            labs(title = "A",
                subtitle = expression(F(X) == b[0]+b[1]*x^(Batch)))            



True_f2_tidy$source <- "True"
Lambda_General_tidy$source <- "Predicted"
General_tidy_lgpr$source <- "Lgpr"

combined_tidy_f2 <- rbind(True_f2_tidy, Lambda_General_tidy,General_tidy_lgpr)



f2_plot <- ggplot(combined_tidy_f2, aes(x = Day, y = mean,group = source)) +  
            # facet_wrap(~coord, scales="free_y",ncol=4)+
            geom_line(aes(color = source),size = 0.7)+
            geom_ribbon(aes(ymin = p2.5, ymax=p97.5,fill= source, group = source),alpha=0.2)+
            geom_ribbon(aes(ymin =p25,ymax=p75,fill= source, group = source),alpha=0.5)+
            scale_fill_manual(values = custom_colors) +        
            scale_color_manual(values = custom_colors)+
            theme_minimal() +
            theme(legend.position = "none") + 
            theme(axis.title.y =element_text(size = 20,margin = margin(t = 15)),
                    axis.title.x =element_text(size = 20,margin = margin(t = 15)),
                    axis.text.x = element_text(angle = 0, hjust = 1,size =17),
                    axis.text.y = element_text(angle = 0, hjust = 1,size =17),
                    plot.title = element_text(size = 27,margin = margin(b = 20)),
                    plot.subtitle = element_text(size =20,margin = margin(b = 25))) +
            labs(y = "F",x= "time")+
            labs(subtitle = expression(F(X) == f^(General)))



True_f3_tidy$source <- "True"
Lambda_Treatment_tidy$source <- "Predicted"
Treatment_tidy_lgpr$source <- "Lgpr"

combined_tidy_f3 <- rbind(True_f3_tidy, Lambda_Treatment_tidy,Treatment_tidy_lgpr)


f3_plot <- ggplot(combined_tidy_f3, aes(x = Day, y = mean,group=source)) +  
            # facet_wrap(~coord, scales="free_y",ncol=4)+
            geom_line(aes(color = source),size = 0.7)+
            geom_ribbon(aes(ymin = p2.5, ymax=p97.5,fill= source, group = source),alpha=0.2)+
            geom_ribbon(aes(ymin =p25,ymax=p75,fill= source, group = source),alpha=0.5)+
            scale_fill_manual(values = custom_colors) +        
            scale_color_manual(values = custom_colors)+
            theme_minimal() +
            theme(legend.position = "none") + 
            theme(axis.title.y =element_text(size = 20,margin = margin(t = 15)),
                    axis.title.x =element_text(size = 20,margin = margin(t = 15)),
                    axis.text.x = element_text(angle = 0, hjust = 1,size =17),
                    axis.text.y = element_text(angle = 0, hjust = 1,size =17),
                    plot.title = element_text(size = 27,margin = margin(b = 20)),
                    plot.subtitle = element_text(size =20,margin = margin(b = 25))) +
            labs(y = "F",x= "time")+
            labs(subtitle = expression(F(X) == f^(Intervention)))

layout_matrix <- rbind(c(1, 2, 3, 4))

p <- grid.arrange(f1_plot,f3_plot,f2_plot,F_plot,layout_matrix=layout_matrix, widths = c(1.1, 1.3, 1.8,2))



