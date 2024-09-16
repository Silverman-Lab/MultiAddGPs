
source('/home/thcodelia/Odelia/Project/MultiAddGPs/src/00_package_loading.R', chdir = TRUE)
source('/home/thcodelia/Odelia/Project/MultiAddGPs/src/helper_functions.R', chdir = TRUE)
load("/home/thcodelia/Odelia/Project/MultiAddGPs/data/mallard_family.RData")

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
           spike = ifelse(time > "2015-11-11 15:00:00" & time < "2015-11-19 16:00:00 ",1,0),
           burnout  = ifelse(time < "2015-11-11 15:00:00",1,0),
           Vessel1Spike = ifelse(Vessel == 1 & spike == 1,1,0),
           Vessel2Spike = ifelse(Vessel == 2 & spike == 1,1,0),
           Vessel1 = ifelse(Vessel == 1,1,0),
           Vessel2 = ifelse(Vessel == 2,1,0),
           Vessel3 = ifelse(Vessel == 3,1,0),
           Vessel4 = ifelse(Vessel == 4,1,0),
           Vessel1Burn = ifelse(Vessel == 1 & burnout != 1,1,0),
           Vessel2Burn = ifelse(Vessel == 2 & burnout != 1,1,0)
           )%>%
    # filter(!duplicated(Hour_diff))%>% # 11-22 and 11-29 have replicates, currently not consider the replicates
    arrange(Vessel)
    # split(.$Vessel)


#match sample_data to the otu_table
mallard_family$otu_table <- mallard_family$otu_table[mallard_family$sample_data$X.SampleID,] 


# Onehot encoding for the Batch and select X
X1 <- mallard_family$sample_data %>%
               dplyr::select(X.SampleID,Vessel,Vessel1,Vessel2,Vessel3,Vessel4,Vessel1Spike,Vessel2Spike,Vessel1Burn,Vessel2Burn,Hour_diff,Hour)

# Hour<- X1$Hour_sample
# Daily <- X1$Daily_sample
Hour <- X1$Hour
Hour_diff <- X1$Hour_diff
Vessel <- paste0("Vessel ",X1$Vessel)
Vessel1 <- X1$Vessel1
Vessel2 <- X1$Vessel2
Vessel3 <- X1$Vessel3
Vessel4 <- X1$Vessel4
Vessel1Spike <- X1$Vessel1Spike
Vessel2Spike <- X1$Vessel2Spike
Veseel1Burn <- X1$Vessel1Burn
Vessel2Burn <- X1$Vessel2Burn

X<- t(model.matrix(~ Vessel1+Vessel2+Vessel3+Vessel4+Vessel1Spike+Vessel2Spike+Veseel1Burn+Vessel2Burn+Hour_diff,data = X1)) # Transpose to get the right dimensions

#  Extract Data / dimensions from Phyloseq object
Y <- t(as(mallard_family$otu_table, "matrix"))
rownames(Y) <- mallard_family$tax_table$Family

if (identical(colnames(Y), X1$X.SampleID)) {
  colnames(Y) <- colnames(X)
} else {
  warning("The order of the samples in the otu_table and sample_data are not the same")
}

