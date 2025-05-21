source('~/Lifestyle_microbiome/02_analysis.R', chdir = TRUE)

# Load required packages
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)

#################################### fit the model ########################################################

sigma_travel_o <- getBestPars(optObj)$sigma_travel
sigma_general_o <- getBestPars(optObj)$sigma_general
rho_travel_o <- getBestPars(optObj)$rho_travel
rho_general_o <- getBestPars(optObj)$rho_general

fit <- model(Y,X,upsilon,Xi,D,
                sigma_travel = sigma_travel_o,rho_travel = rho_travel_o,
                sigma_general = sigma_general_o,rho_general = rho_general_o)


str(fit)
fit.clr <- to_clr(fit)
str(fit.clr)

####################################### Validate the result using the validate data #########################################

fit = fit.clr
predictX = X
predictY = Y
predictTaxa_name = taxa_names
predictmetadata = train_metadata_st


Predict <- predict(fit,newdata = predictX, reponse="Lambda")

Predict_tidy <- gather_array(Predict, val, coord, sample, iter)%>%
              mutate(Day = predictX[2,sample],
                     Original_Day = predictmetadata$COLLECTION_DAY[sample])%>%
              filter(!is.na(val))%>%
              group_by(Original_Day,coord)%>%
              summarise_posterior(val,na.rm = TRUE)%>%
              ungroup()%>%
              mutate(coord = paste0("CLR(",predictTaxa_name[coord], ")"))%>%
              filter(coord %in% c("CLR(Fusobacteriaceae)",
                                   "CLR(Monoglobaceae)","CLR(Yersiniaceae)"))

predictY_tidy <- clr_array(predictY+0.65, parts = 1) %>% 
       gather_array(mean, coord, sample) %>% 
       mutate(Day = predictX[2,sample],
              Original_Day = predictmetadata$COLLECTION_DAY[sample],
              coord =paste0("CLR(",predictTaxa_name[coord], ")"))%>%
       filter(coord %in% c("CLR(Fusobacteriaceae)",
                            "CLR(Monoglobaceae)","CLR(Yersiniaceae)"))

       

custom_colors <- c('#5ba300','#f57600','#8babf1','#0073e6')
custom_colors_v12 <- c('#5ba300','#f57600')


F_plot <- ggplot(Predict_tidy, aes(x = Original_Day, y = mean)) +  
              facet_wrap(~coord, scales="free_y", ncol=1, strip.position = "right" ) +
              geom_line(color = "blue", size = 1) +
              # geom_point(data = predictY_tidy, 
              #         color = "#56B4E9", alpha = 0.4) +  # Light blue points
              geom_ribbon(aes(ymin = p2.5, ymax = p97.5), fill = "#377EB8", alpha = 0.3) +  # Darker blue for 95% CI
              geom_ribbon(aes(ymin = p25, ymax = p75), fill = "#377EB8", alpha = 0.8) +  # Darker blue for 50% CI
              theme(
                     panel.background = element_rect(fill = "white", color = NA),
                     plot.background = element_rect(fill = "white", color = NA),
                     panel.grid.major = element_line(color = "gray90"),
                     panel.grid.minor = element_line(color = "gray95"),
                     strip.background = element_rect(fill = "white", color = "black"),
                     strip.text = element_text(color = "black", size = 13),
                     panel.border = element_rect(color = "black", fill = NA, size = 1),
                     panel.spacing = unit(1, "lines"),  # Increase spacing between facets
                     axis.title.y = element_text(size = 15, margin = margin(t = 15)),
                     axis.title.x = element_text(size = 15, margin = margin(t = 15)),
                     axis.text.x = element_text(angle = 0, hjust = 1, size = 15),
                     axis.text.y = element_text(angle = 0, hjust = 1, size = 15),
                     plot.subtitle = element_text(size = 18, margin = margin(b = 20),hjust = 0.5),
                     legend.spacing.y = unit(16, "pt"),
                     legend.key.height = unit(1.5, "cm")) +
              labs(subtitle = expression(F(X)==f^(trend)+f^(travel)),
                     x = "Day", y = "F")

## plot each Lambda
intercept <- fit.clr$Lambda[1][[1]]

Lambda_Travel <- array(0,dim = dim(fit.clr$Lambda[2][[1]]))
Lambda_General <-  array(0,dim = dim(fit.clr$Lambda[2][[1]]))

## Add the intercept to the Lambda
for(i in 1:dim(fit.clr$Lambda[2][[1]])[3]){
     Lambda_Travel[,,i] <- fit.clr$Lambda[2][[1]][,,i] + intercept[,,i]}

for(i in 1:dim(fit.clr$Lambda[2][[1]])[3]){
     Lambda_General[,,i] <- fit.clr$Lambda[3][[1]][,,i] + intercept[,,i]}

# Center the Lambda and Y
Lambda_Travel<-standardize_function(Lambda_Travel,iteration = TRUE)
Lambda_General<-standardize_function(Lambda_General,iteration = TRUE)
Y_st <- standardize_function(clr_array(Y+0.65, parts = 1),iteration = FALSE)

Lambda_General_tidy <- gather_array(Lambda_General, val, coord, sample,iter) %>% 
                    mutate(Day = X[2,sample],
                           Original_Day = train_metadata_st$COLLECTION_DAY[sample])%>%
                    filter(!is.na(val))%>%
                    group_by(Original_Day,coord)%>%
                    summarise_posterior(val,na.rm = TRUE)%>%
                    ungroup()%>%
                    mutate(coord = paste0("CLR(",taxa_names[coord], ")"))%>%
                    filter(coord %in% c("CLR(Fusobacteriaceae)",
                                        "CLR(Monoglobaceae)","CLR(Yersiniaceae)"))

Lambda_Travel_tidy <- gather_array(Lambda_Travel, val, coord, sample,iter) %>%
                    mutate(Day = X[2,sample],
                           Original_Day = train_metadata_st$COLLECTION_DAY[sample])%>%
                    filter(!is.na(val))%>%
                    group_by(Original_Day,coord)%>%
                    summarise_posterior(val,na.rm = TRUE)%>%
                    ungroup()%>%
                    mutate(coord = paste0("CLR(",taxa_names[coord], ")"))%>%
                    filter(coord %in% c("CLR(Fusobacteriaceae)",
                                         "CLR(Monoglobaceae)","CLR(Yersiniaceae)"))

Y_training_tidy <- Y_st %>% 
            gather_array(val, coord, sample) %>% 
            mutate(Day = X[2,sample],
                   Original_Day = train_metadata_st$COLLECTION_DAY[sample], 
                   coord =paste0("CLR(",taxa_names[coord], ")"))%>%
            filter(coord %in% c("CLR(Fusobacteriaceae)",
                                "CLR(Monoglobaceae)","CLR(Yersiniaceae)"))


custom_colors <- c('#5ba300','#f57600','#8babf1','#0073e6')
custom_colors_v12 <- c('#5ba300','#f57600')


General_plot <- ggplot(Lambda_General_tidy, aes(x = Original_Day, y = mean)) +  
                        facet_wrap(~coord, scales="free_y", ncol=1) +
                        geom_line(color = "blue", size = 1) +
                        # geom_point(data = Y_training_tidy, aes(x = Day, y = val), 
                        #         color = "#56B4E9", alpha = 0.4) +  # Light blue points
                        geom_ribbon(aes(ymin = p2.5, ymax = p97.5), fill = "#377EB8", alpha = 0.3) +  # Darker blue for 95% CI
                        geom_ribbon(aes(ymin = p25, ymax = p75), fill = "#377EB8", alpha = 0.8) +  # Darker blue for 50% CI
                        theme(
                            panel.background = element_rect(fill = "white", color = NA),
                            plot.background = element_rect(fill = "white", color = NA),
                            panel.grid.major = element_line(color = "gray90"),
                            panel.grid.minor = element_line(color = "gray95"),
                            strip.background = element_blank(),
                            strip.text = element_blank(), 
                            panel.spacing = unit(1, "lines"),  # Increase spacing between facets
                            panel.border = element_rect(color = "black", fill = NA, size = 1),
                            axis.title.y = element_text(size = 15, margin = margin(t = 15)),
                            axis.title.x = element_text(size = 15, margin = margin(t = 15)),
                            axis.text.x = element_text(angle = 0, hjust = 1, size = 15),
                            axis.text.y = element_text(angle = 0, hjust = 1, size = 15),
                            plot.subtitle = element_text(size = 18, margin = margin(b = 20),hjust = 0.5),
                            legend.spacing.y = unit(16, "pt"),
                            legend.key.height = unit(1.5, "cm")
                            ) +
                        labs( y = "F",
                              x = "Day",
                             subtitle = expression(F(X)==f^(trend)))


Travel_plot <-ggplot(Lambda_Travel_tidy, aes(x = Original_Day, y = mean)) +  
                    facet_wrap(~coord, scales="free_y", ncol=1) +
                    geom_line(color = "blue", size = 1) +
                    # geom_point(data = Y_training_tidy, aes(x = Day, y = val), 
                    #         color = "#56B4E9", alpha = 0.2) +  # Light blue points
                    geom_ribbon(aes(ymin = p2.5, ymax = p97.5), fill = "#377EB8", alpha = 0.3) +  # Darker blue for 95% CI
                    geom_ribbon(aes(ymin = p25, ymax = p75), fill = "#377EB8", alpha = 0.8) +  # Darker blue for 50% CI
                    theme(
                        panel.background = element_rect(fill = "white", color = NA),
                        plot.background = element_rect(fill = "white", color = NA),
                        panel.grid.major = element_line(color = "gray90"),
                        panel.grid.minor = element_line(color = "gray95"),
                        strip.background = element_blank(),
                        strip.text = element_blank(), 
                        panel.spacing = unit(1, "lines"),  # Increase spacing between facets
                        panel.border = element_rect(color = "black", fill = NA, size = 1),
                        axis.title.y =element_text(size = 15, margin = margin(t = 15)),
                        axis.title.x = element_text(size = 15, margin = margin(t = 15)),
                        axis.text.x = element_text(angle = 0, hjust = 1, size = 15),
                        axis.text.y = element_text(angle = 0, hjust = 1, size = 15),
                        plot.subtitle = element_text(size = 18, margin = margin(b = 20),hjust = 0.5),
                        legend.spacing.y = unit(16, "pt"),
                        legend.key.height = unit(1.5, "cm")
                    ) +
                   labs(subtitle = expression(F(X)==f^(travel)),
                        x = "Day")

combined_plot <- General_plot + Travel_plot +F_plot   # Arrange in 1 row, 3

combined <- combined_plot + plot_layout(ncol = 1, nrow = 3)  # Ensure 4 rows, 3 columns
