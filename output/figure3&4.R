
source('/home/thcodelia/Odelia/Project/MultiAddGPs/src/02_analysis.R', chdir = TRUE)

hyper <- c(getBestPars(optObj)$sigma1,getBestPars(optObj)$l1,
           getBestPars(optObj)$sigma2,getBestPars(optObj)$l2,2,
           1)
           
mod <- Model(Y,X,upsilon,sigma1=hyper[1],l1=hyper[2],
                         sigma2=hyper[3],l2=hyper[4], 
                         a2 = hyper[5],
                         noise = hyper[6])




## Prediction: plot in CLR space  ##########
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
predicted_clr_tidy <- predicted_clr_tidy %>%
  left_join(facet_lims, by = "coord")

p <- ggplot(predicted_clr_tidy,aes(x=Hour,y=mean))+
      facet_wrap(~coord, scales="free_y", ncol=5) +
      geom_ribbon(aes(ymin = p2.5, ymax = p97.5,fill= Vessel), alpha = 0.4) +
      # geom_ribbon(aes(ymin = p25,  ymax = p75,  fill= Vessel), alpha = 0.4) +
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
        strip.text = element_text(color = "black",size = 10),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), 
        axis.text.x = element_text(angle = 0, hjust = 1),
        plot.title = element_text(size = 10),
        legend.title = element_blank()) +   
      labs(y = expression(Lambda), x= "Hour")+
      annotate("rect", xmin = 240, xmax = 264, ymin = -Inf, ymax = Inf, fill = "#717D7E", alpha = 0.15) 


# Lambda Decomposition
combine_vessel_summaries <- function(mod_clr, lambda_component,filter_lambda_vessel1, filter_lambda_vessel2, time_column, coords, CLR = TRUE) {
  combined_summary <- data.frame()  # Initialize an empty dataframe to store combined summaries
  
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
    vessel1_summary$Vessel <- "Vessel 1"
    vessel2_summary$Vessel <- "Vessel 2"
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


disruption_plot <- ggplot(data = dis_allcoord, aes(x = time, y = mean)) +
                      facet_wrap(~Coord, scales="free_y", ncol=5) +
                      geom_ribbon(aes(ymin = p2.5, ymax = p97.5,fill = Vessel), alpha = 0.2) +
                      geom_ribbon(aes(ymin = p25, ymax = p75,   fill = Vessel), alpha = 0.5) +
                      geom_line(aes(color=Vessel),alpha = 0.4, size = 0.7) +
                      geom_hline(yintercept = 0, color = "black", linetype = "dashed") +  
                      scale_fill_manual(values = custom_colors_v12) +          
                      scale_color_manual(values = custom_colors_v12) +  
                      # Additional customizations
                      theme(
                              legend.position = "none",
                              panel.background = element_rect(fill = "white", color = NA),
                              plot.background = element_rect(fill = "white", color = NA),
                              panel.grid.major = element_line(color = "gray90"),
                              panel.grid.minor = element_line(color = "gray95"),
                              axis.title.y = element_text(size = 8),  
                              axis.title.x = element_text(size = 8), 
                              axis.text.x = element_text(angle = 0, hjust = 1),
                              plot.title = element_text(size = 16),
                              panel.border = element_rect(color = "black", fill = NA, size = 1) 
                          )+
                      labs(x= "Hour")
