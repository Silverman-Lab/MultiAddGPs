source('~/src/loading_package.R', chdir = TRUE)
source('~/src/source.R')
# source('/home/thcodelia/Odelia/Project/Project1_AGP/results/daivd_timeseries/1_27_2025/Model training.R')


# load the filtered data, 98 taxa with 328 sample and 41% zero count
fil_data <- readRDS("~/data/lifestyle.rds")


####################################### Split train, validate, and test data #######################################

metadata <- data.frame(sample_data(fil_data)) %>%
            arrange(COLLECTION_DAY)

### split the data by equally seperate the time sequence to 33 groups and randomaly select 2 samples from each group
metadata$group <- cut(seq_len(nrow(metadata)),breaks = 33, labels = FALSE)

validate_indices <- c()
test_indices <- c()
set.seed(1040)
for (i in 1:33) {
    group_samples <- metadata %>% filter(group == i)
    sampled <- sample(rownames(group_samples),2)

    validate_indices <- c(validate_indices, sampled[1])
    test_indices <- c(test_indices, sampled[2])
}

# ### split the data by time sequence
validate_data <- prune_samples(c(validate_indices,test_indices), fil_data) %>% tax_glom("Family")
test_data <- prune_samples(test_indices, fil_data) %>% tax_glom("Family")
train_indices <- setdiff(rownames(metadata), c(validate_indices, test_indices))
train_data <- prune_samples(train_indices, fil_data)
train_data <- train_data %>% tax_glom("Family")


###################################### standardize the metadata #######################################
# Standardize the metadata
train_metadata_st <- sample_data(train_data) %>%
                    data.frame() %>%
                    # mutate(Travel = as.factor(Travel)) %>%
                    mutate(COLLECTION_DAY = as.numeric(as.character(COLLECTION_DAY))) %>%
                    mutate(Day = (COLLECTION_DAY - mean(COLLECTION_DAY, na.rm = TRUE)) / sd(COLLECTION_DAY, na.rm = TRUE)) %>%
                    # mutate(Day = as.numeric(COLLECTION_DAY)) %>%
                    mutate(Travel = warp(COLLECTION_DAY,mode=92,scale =20,skew=2)) %>%
                    mutate(across(
                        starts_with("NUTRITION_"),
                        ~ as.numeric(as.character(.))  # Convert to numeric first
                    )) %>%
                    mutate(across(
                        starts_with("NUTRITION_"),
                        ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)  # Standardize
                    ))%>%
                    dplyr:: select(Day,COLLECTION_DAY, ID, Travel, NUTRITION_CALCIUM_PRECEDING_DAY:NUTRITION_SUGAR_PRECEDING_DAY)
          
# plot(train_metadata_st$Day, is.na(train_metadata_st$NUTRITION_CALCIUM_PRECEDING_DAY))
# plot(train_metadata_st$Day,train_metadata_st$Travel)

validate_data_st <- sample_data(validate_data) %>%
                    data.frame() %>%
                    mutate(Travel = as.factor(Travel)) %>%
                    mutate(COLLECTION_DAY = as.numeric(as.character(COLLECTION_DAY))) %>%
                    mutate(Day = (COLLECTION_DAY - mean(COLLECTION_DAY, na.rm = TRUE)) / sd(COLLECTION_DAY, na.rm = TRUE)) %>%
                    # mutate(Day = as.numeric(COLLECTION_DAY)) %>%
                    mutate(Travel = warp(COLLECTION_DAY,mode=92,scale =20,skew=2)) %>%
                    mutate(across(
                        starts_with("NUTRITION_"),
                        ~ as.numeric(as.character(.))  # Convert to numeric first
                    )) %>%
                    mutate(across(
                        starts_with("NUTRITION_"),
                        ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)  # Standardize
                    ))%>%
                    dplyr:: select(Day, COLLECTION_DAY, ID, Travel, NUTRITION_CALCIUM_PRECEDING_DAY:NUTRITION_SUGAR_PRECEDING_DAY)
                    # arrange(Day)


X <- t(model.matrix(~ Day + Travel, data = train_metadata_st))
Y <- as(otu_table(train_data), "matrix")
rownames(Y) <- data.frame(tax_table(train_data))$Family   
taxa_names <- as(data.frame(tax_table(train_data))$Family,"vector")

if(any(colnames(Y) != colnames(X))){
    stop("The column names of the training metadata and the OTU table do not match")
}

X_validate <- t(model.matrix(~ Day + Travel , data = validate_data_st))
Y_validate <- as(otu_table(validate_data), "matrix")
rownames(Y_validate) <- data.frame(tax_table(validate_data))$Family
taxa_names_validate <- as(data.frame(tax_table(validate_data))$Family,"vector")

if(any(colnames(Y_validate) != colnames(X_validate))){
    stop("The column names of the validated metadata and the OTU table do not match")
}

optObj_AGP <- bayesOpt(
    FUN= function(sigma1,l1,sigma2,l2
                 ){
                  MLL_sudo(Y,X,Sigma,D,
                           sigma1,l1,p = 20,
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
                  init_points = 10,
                  n_iter = 20,
                  acq = "ucb",
                  kappa = 2.576,
                  verbose = 1
)


