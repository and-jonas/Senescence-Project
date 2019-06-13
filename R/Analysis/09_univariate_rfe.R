#====================================================================================== -

#HEADER ----

# AUTHOR: Jonas Anderegg

#====================================================================================== -

.libPaths("T:/R3UserLibs")

library(tidyverse)
library(caret)
library(mlbench)
library(Hmisc)
library(randomForest)
library(corrplot)

#set working directory
wd <- "O:/Projects/KP0011/1/Senescence-Project/"
setwd(wd)

#load required function
source("R/Utils/003_univar_rfe_utils.R")

#define source and sink directories
dirfrom <- "Analysis/ANOVA/"
# dirto <- "Analysis/ANOVA/"

#====================================================================================== -

#Data preparation ----

# Load BLUEs for GY, GPC and DynPars

data <- readRDS(paste(dirfrom, "Expwise.rds", sep = ""))

BLUE <- data %>% dplyr::select(Exp, Trait, BLUE) %>% 
  filter(!is.na(BLUE)) %>% unnest() %>% 
  spread(., Trait, BLUE)

#====================================================================================== -

#UNIVARIATE ANALYSIS ----

# read data
dynpars <- readRDS("T:/PhD/Data_Analysis_3/Trait_pred/Data/dynpars_decorr.rds")
dynpars <- c(dynpars, "onsen_lin_SI_NDVI_nb_ASD",
             "midsen_lin_SI_NDVI_nb_ASD",
             "endsen_lin_SI_NDVI_nb_ASD",
             "tsen_lin_SI_NDVI_nb_ASD")

#uni-variate analysis
xy <- BLUE %>%
  #select all index derived dynamics parameters
  dplyr::select(one_of(dynpars),
                #select all flag leaf based parameters
                contains("SnsFl0"), 
                #select all canopy based parameters
                contains("SnsCnp"),
                #select response variable
                GY, GPC, Exp)
  
do_cor_test <- function(data){
  cor <- cor.test(data$value, data$response) %>%
    broom::tidy()
}

# perform linear regression for each predictor
# extract r squared and p values of linear regressions
data <- xy %>% 
  #transform data to long format
  gather(sendynpar, value, contains("sen"), factor_key = TRUE) %>% 
  gather(Response, response, GY:GPC, factor_key = TRUE) %>%
  #exclude obs with missing data
  filter(complete.cases(.)) %>%
  #reorder columns
  dplyr::select(Exp, sendynpar, Response, value, response) %>%
  #create nested tibble
  nest(value, response) %>% as_tibble() %>% 
  mutate(cor_test = purrr::map(data, do_cor_test)) %>%
  mutate(estimate = purrr::map_dbl(cor_test, "estimate"),
         pval = purrr::map_dbl(cor_test, "p.value")) %>%
  arrange(Exp, Response, desc(abs(estimate))) %>% 
  group_by(Exp, Response) %>% 
  top_n(n = 3, wt = abs(estimate))

out <- data %>% select(1:3, estimate, pval) %>%
  rowwise() %>%
  mutate(sig = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "")))) %>%
  transmute(Exp, sendynpar, Response, Est = paste(round(estimate, 3), sig, sep = "")) 
    

    # #fit a linear model response~x where x is each individual sendynpar
    # #and response is either GY or GPC
    # mutate(fit_lin = purrr::map(data, lm)) %>%
    # #extract the R squared and p-values of the linear regressions
    # mutate(glance = purrr::map(fit_lin, broom::glance),
    #        rsq = purrr::map_dbl(glance, "r.squared"),
    #        pval = purrr::map_dbl(glance, "p.value")) %>%
    # arrange(Year, Response, desc(rsq))
  
  out <- data %>% select(1:3, rsq, pval) %>%
    rowwise() %>%
    mutate(sig = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "")))) %>%
    transmute(Year, sendynpar, Response, Rsq = paste(round(rsq, 3), sig, sep = "")) 
  
  #correct index names as possible
  out$sendynpar <- gsub(pattern = "SI_", x = out$sendynpar, replacement = "")
  out$sendynpar <- gsub(pattern = "_r", x = out$sendynpar, replacement = "")
  out$sendynpar <- gsub(pattern = "_lin", x = out$sendynpar, replacement = "")

  write.csv(out, "O:/Projects/KP0011/1/Analysis/RFE/univar_cor.csv", row.names = FALSE)  
  
  # perform linear regression for each predictor
  # extract r squared and p values of linear regressions
  
  # read data
  BLUEs <- readRDS("T:/PhD/DataAnalysis/FPWW012_FPWW018/BLUEs_3Y.rds")
  dynpars <- readRDS("T:/PhD/Data_Analysis_3/Trait_pred/Data/dynpars_decorr.rds")
  dynpars <- c(dynpars, "onsen_lin_SI_NDVI_nb_ASD",
               "midsen_lin_SI_NDVI_nb_ASD",
               "endsen_lin_SI_NDVI_nb_ASD",
               "tsen_lin_SI_NDVI_nb_ASD")
  
  SIs <- unique(sapply(strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', dynpars), ' '), "[[", 2))
  names <- names(BLUEs)[grepl(paste0(SIs, collapse = "|"), names(BLUEs))]
  
  #uni-variate analysis
  xy <- BLUEs %>%
    select(#select Gen and Year
           Gen, Year, 
           #select all index derived dynamics parameters
           one_of(names),
           #select all flag leaf based parameters
           contains("SnsFl0"), 
           #select all canopy based parameters
           contains("SnsCnp")) %>%
    filter(Year == 2016) %>%
    select(-contains("gom"))

  data <- xy %>% 
    #transform data to long format
    tidyr::gather(sendynpar, value, endsen_lin_SI_1200_r:tsen_lin_SnsCnp, factor_key = TRUE)
 
  data$index <- sapply(strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', data$sendynpar), ' '), "[[", 2)
  data$sendynpar <- sapply(strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', data$sendynpar), ' '), "[[", 1)
  
  data <- data %>% select(Gen, Year, index, sendynpar, value)
  
  do_cor_test <- function(data){
    cor <- cor.test(data$onsen_lin, data$tsen_lin) %>%
      broom::tidy()
  }

  dat <- data %>%
    tidyr::spread(key = sendynpar, value = value) %>%
    group_by(Year, index) %>%
    nest() %>%
    mutate(cor_test = purrr::map(data, do_cor_test)) %>%
    mutate(estimate = purrr::map_dbl(cor_test, "estimate"),
           pval = purrr::map_dbl(cor_test, "p.value")) %>%
    arrange(Year, desc(abs(estimate)))
  
  mod <- lm(GY ~ tsen_lin_SI_NPCI_r, data = xy)
  summary(mod)
  
  cor(xy$tsen_lin_SI_NPCI_r, xy$tsen_lin_SI_PSRI_r, use = "complete.obs")
  plot(xy$tsen_lin_SI_NPCI_r, xy$tsen_lin_SI_PSRI_r)
  
  
  d16 <- xy[xy$Year == 2016,]
  
  mod0 <- lm(GY ~ onsen_lin_SI_PSRI_r, data = d16)
  preds <- predict(mod0, d16)
  
  # function to calculate RMSE  
  rmse = function(actual, predicted) {
    sqrt(mean((actual - predicted) ^ 2, na.rm = TRUE))
  }
  
  rmse <- rmse(d16$GY, preds)
  
  sub <- d16 %>% select(contains("gom"))
  
  M <- cor(sub, use = "pairwise.complete.obs")
  corrplot(M, method = "color")
  
  onsen_pars <- names(cleaned)
  
  zdf <- as.data.frame(as.table(M))
  zdf
  
  subset(zdf, abs(Freq) > 0.96)
  
  #assess coefficients of variation in the BLUEs
  d <- BLUEs %>%           
    #select response variable
    select(GY, GPC, Year) %>%
    #transform to long
    tidyr::gather(Trait, value, GY:GPC, factor_key = TRUE) %>%
    arrange(Year) %>%
    group_by(Year, Trait) %>%
    nest(-c(Year, Trait)) %>%
    mutate(cv = purrr::map(data, function(x) sd(as.numeric(unlist(x))), na.rm = TRUE)) %>%
    unnest(cv)
    
  data <- xy %>% 
    tibble::as.tibble()%>%
    filter(Year == 2016) 
  
  ggplot(data, aes(x = midsen_lin_SI_NDWI1, y = GY)) +
    geom_point()

#====================================================================================== -

#RF-RFE ----
  
## An initial split of the full data set is created 
## using stratified sampling based on percentiles of the outcome;
## This results in a training set of 260 samples and an independant validation set of 73 samples;
## Then, the training dataset is randomly resampled 25 times (80:20) to produce 25 different train-test sets.
## For each resample, recursive feature elimination is carried out

# read data
  BLUEs <- readRDS("T:/PhD/DataAnalysis/FPWW012_FPWW018/BLUEs_3Y.rds")
  dynpars <- readRDS("T:/PhD/Data_Analysis_3/Trait_pred/Data/dynpars_decorr.rds")
  dynpars <- c(dynpars, "onsen_lin_SI_NDVI_nb_ASD",
               "midsen_lin_SI_NDVI_nb_ASD",
               "endsen_lin_SI_NDVI_nb_ASD",
               "tsen_lin_SI_NDVI_nb_ASD")
  


# select only 1 year data to start with
  xy <- BLUEs %>%
    filter(Year == "2017") %>%
    #select all index derived dynamics parameters
    select(one_of(dynpars),
           #select all flag leaf based parameters
           contains("SnsFl0"),
           #select all canopy based parameters
           contains("SnsCnp"),
           #select response variable
           GPC) %>%
    #exclude obs with missing data
    filter(complete.cases(.))

# exclude extreme values, obvious outliers
  xy <- xy[!xy$GY < 3, ]

# perform linear regression for each predictor
# extract r squared and p values of linear regressions
  data <- xy %>% 
    #transform data to long format
    tidyr::gather(sendynpar, value, onsen_lin_SI_780_700:tsen_lin_SnsCnp, factor_key = TRUE) %>% 
    #reorder columns
    select(sendynpar, value, GPC) %>%
    #create nested tibble
    tidyr::nest(value, GPC) %>%
    tibble::as.tibble() %>%
    #fit a linear model GY~x where x is each individual sendynpar
    mutate(fit_lin = purrr::map(data, lm)) %>%
    #extract the R squared and p-values of the linear regressions
    mutate(glance = purrr::map(fit_lin, broom::glance),
           rsq = purrr::map_dbl(glance, "r.squared"),
           pval = purrr::map_dbl(glance, "p.value")) %>%
    arrange(desc(rsq))

#Prepare data for rfe 

          # # create df holding the predictor data
          #   preds <- xy[-length(xy)]
          # # create vector holding the response data
          #   gy <- xy["GY"] %>% as.matrix()
          # # split data into groups based on percentiles
          # # sample within these subgroups (10 here) (stratified samplings). 
          #   split <- createDataPartition(gy, p = .8, list = FALSE, groups = 10)
          #   
          #   #RECORD THE SPLIT!!!
          #   write.csv(split, file = "O:/Projects/KP0011/1/Analysis/RFE/split_GPC_2017.csv")
          #   
          # # integrate data into single df
          #   dat <- preds
          #   dat$gpc <- gpc
          # # create training and testing datasets using the created resampling indices from split()
          #   train0 <- dat[ split, ]
          #   validate0  <- dat[-split, ]
    
#OR:: NO SPLIT
    
    train0 <- xy
    names(train0)[length(train0)] <- "gpc"
  
  # create the cross-validation files as a list to use with different 
  # functions
    index <- createMultiFolds(train0$gpc, k = 5, times = 6)

  # The candidate set of the number of predictors to evaluate
    subsets <- c(90, 85, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 17, 14, 12:1)
    # subsets <- c(20, 10, 2:1)
    
  # function to calculate RMSE  
  rmse = function(actual, predicted) {
    sqrt(mean((actual - predicted) ^ 2))
  }

################################################################################################# - 
################################################################################################# -
  
## PERFORM RERUSIVE FEATURE ELIMINATION 
## EMBEDDED IN AN OUTER RESAMPLING LOOP
  
  start.time <- Sys.time()  
  
  # OUTER RESAMPLING
  # loop over all resamples
  
    out <- list()
  
    for(i in 1:length(index)){
      
      print(paste("resample ", i, "/", length(index), sep = ""))
      
      #use indices to create train and test data sets for the resample
      ind <- as.numeric(index[[i]])
      train <- train0[ind,]
      test <- train0[-ind, ]
      
      keep_vars <- list()
      test_rmse <- NULL
      train_rmse <- NULL
      varimp_all <- vector(length = length(subsets)+1, mode = "list")
      
### TUNE/TRAIN RANDOM FOREST USING ALL PREDICTORS; CALCULATE VIP AND TEST RMSE ###

      #adjust mtry parameter to decreasing predictor set
      mtry <- c(ceiling((length(train)-1)/10), ceiling((length(train)-1)/5), 
                ceiling((length(train)-1)/3), ceiling((length(train)-1)/1.5))
      
      #specify model tuning parameters
      tune_grid <- expand.grid(mtry = mtry,
                               splitrule = "variance",
                               min.node.size = 5)
      
      #define inner resampling procedure
      ctrl <- trainControl(method = "repeatedcv",
                           number = 10,
                           rep = 1,
                           verbose = FALSE,
                           allowParallel = TRUE)
      
      #random forest
      rf_ranger <- train(gpc~.,
                         data = train,
                         preProc = c("center", "scale"),
                         method = "ranger",
                         tuneGrid = tune_grid,
                         importance = "permutation",
                         num.trees = 10000,
                         num.threads = 10, 
                         trControl = ctrl)
      
      #get variable importance and order variables according to importance
      vip <- varImp(rf_ranger)$importance
      vip$vars <- rownames(vip)
      vip <- vip[order(vip$Overall, decreasing = TRUE), ]

      #assign output to vector
      varimp_all[[1]] <- vip
      keep_vars[[1]] <- vip[1:subsets[1], "vars"] 
      new_vars <- c(vip[1:subsets[1], "vars"], "gpc")
      test_rmse[1] <- rmse(test$gpc, predict(rf_ranger, test))
      train_rmse[1] <- getTrainPerf(rf_ranger)
      
### DONE ###
      
      #reduce train data: 
      #exclude predefined number of most irrelevant predictors 
      #according to sizes sequence
      newtrain <- select(train, new_vars)
      
### FEATURE ELIMINATION LOOP ###

      #perform backwards feature selection using rfRFE
      #based on tuned(mtry parameter, 10-fold CV) rf models
      #sequentially remove predefined number of predictors with the lowest vip score
      for (j in 1:length(subsets)){
        
        print(paste("--> step ", j, "/", length(subsets), sep = ""))
        
        #adjust mtry parameter to decreasing predictor set
        mtry <- c(ceiling((length(newtrain)-1)/10), ceiling((length(newtrain)-1)/5), 
                  ceiling((length(newtrain)-1)/3), ceiling((length(newtrain)-1)/1.5))
        
        #specify model tuning parameters
        tune_grid <- expand.grid(mtry = mtry,
                                 splitrule = "variance",
                                 min.node.size = 5)
        
        #define resampling procedure
        ctrl <- trainControl(method = "repeatedcv",
                             number = 10,
                             rep = 1,
                             verbose = FALSE)
        
        #random forest
        rf_ranger <- train(gpc~.,
                           data = newtrain,
                           preProc = c("center", "scale"),
                           method = "ranger",
                           tuneGrid = tune_grid,
                           importance = "permutation",
                           num.trees = 10000,
                           num.threads = 10,
                           trControl = ctrl)
        
        #get variable importance and order variables according to importance
        vip <- varImp(rf_ranger)$importance
        vip$vars <- rownames(vip)
        vip <- vip[order(vip$Overall, decreasing = TRUE), ]
        vip[is.na(vip)] <- 100
        
        #assign output to vector
        varimp_all[[j+1]] <- vip
        if(j != length(subsets)){
          new_vars<- c(vip[1:subsets[j+1], "vars"], "gpc")
          }else NULL
        test_rmse[j+1] <- rmse(test$gpc, predict(rf_ranger, test))  
        train_rmse[j+1] <- getTrainPerf(rf_ranger)

        #reduce train data: 
        #exclude predefined number of most irrelevant predictors 
        #according to sizes sequence
        newtrain <- select(newtrain, new_vars)

      }#END OF FEATURE ELIMINATION LOOP
      
### ASSIGN RANKS TO VARIABLES DEPENDING ON VARIMP ###
      
      #assign ranks to predictors based on the variable importance 
      varranks_all <- list()
      for(k in 1:length(varimp_all)){
        
        #exception for the last iteration...
        if(k !=length(varimp_all)){
          var <- varimp_all[[k]][(subsets[k]+1):(nrow(varimp_all[[k]])), "vars"]
          rank <- rep(length(subsets)-k+2, length(var))
        } else {
          var <- varimp_all[[k]][, "vars"]
          rank <- 1
        }
          
        varranks_all[[k]] <- as.data.frame(cbind(var, rank))
        
      }
      
### TIDY UP results ###
      
      #create tibble holding varriable rankings
      varranks <- tibble::as.tibble(do.call("rbind", varranks_all))
      varranks$rank <- as.numeric(as.character(varranks$rank))
      
      #create tibble hodling RMSE and subsetsize
      subset_size <- c(length(train0)-1, subsets)
      RMSEs <- tibble::as.tibble(cbind(subset_size, train_rmse, test_rmse))
      
      #store in list
      out[[i]] <- list(varranks, RMSEs)
      
    }
    #END OF OUTER RESAMPLING
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
  
  saveRDS(out,"O:/Projects/KP0011/1/Analysis/RFE_New/gpc_2017.rds")
  
################################################################################################# -
################################################################################################# -
  
#Tidy up
allranks <- lapply(out, "[[", 1)
allRMSE <- lapply(out, "[[", 2) %>% lapply(as.data.frame)

ranks <- allranks %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="var"), .)

RMSE <- allRMSE %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="subset_size"), .)

  saveRDS(ranks, "O:/Projects/KP0011/1/Analysis/RFE_New/ranks_GPC_2017.rds")
  saveRDS(RMSE, "O:/Projects/KP0011/1/Analysis/RFE_New/rmse_GPC_2017.rds")

#################################################################################################  
#################################################################################################  
  
## Independent validation
  
  #split
  split <- read.csv("O:/Projects/KP0011/1/Analysis/RFE/split_GY_2016.csv")[-1] %>%
    as.matrix
    
  #features
  feats0 <- read.csv("O:/Projects/KP0011/1/Analysis/RFE/all_ranks.csv") %>%
    tibble::as.tibble()
  
  #Prepare data for validation
  
  # create df holding the predictor data
  preds <- xy[-length(xy)]
  # create vector holding the response data
  gy <- xy["GY"] %>% as.matrix()

  # integrate data into single df
  dat <- preds
  dat$gy <- gy
  # create training and testing datasets using the created resampling indices from split()
  train0 <- dat[split, ]
  validate0  <- dat[-split, ]
  
  # The candidate set of the number of features to evaluate
  subsets <- c(1:12)
  
  test_rmse <- NULL
  
  for (i in 1:length(subsets)) {
    
    #list according featurs
    feats <- feats0 %>%
      filter(year == 2016, trait == "GY") %>%
      select(subset_size)
    #select top i features
    feats_sel <- as.character(unlist(feats[1:i,]))
    
    train <- train0 %>%
      select(one_of(feats_sel), gy) 
    
    train <- train[complete.cases(train),]
    
    ### TUNE/TRAIN RANDOM FOREST USING ALL PREDICTORS; CALCULATE VIP AND TEST RMSE ###
    
    #adjust mtry parameter to decreasing predictor set
    mtry <- c(ceiling((length(train)-1)/10), ceiling((length(train)-1)/5), 
              ceiling((length(train)-1)/3), ceiling((length(train)-1)/1.5))
    
    #specify model tuning parameters
    tune_grid <- expand.grid(mtry = mtry,
                             splitrule = "variance",
                             min.node.size = 5)
    
    #define resampling procedure
    ctrl <- trainControl(method = "repeatedcv",
                         number = 10,
                         rep = 1,
                         verbose = TRUE,
                         allowParallel = TRUE)
    
    #random forest
    rf_ranger <- train(gy~.,
                       data = train,
                       preProc = c("center", "scale"),
                       method = "ranger",
                       tuneGrid = tune_grid,
                       importance = "permutation",
                       num.trees = 1000,
                       num.threads = 16,
                       trControl = ctrl)
    
    test_rmse[i] <- rmse(validate0$gy, predict(rf_ranger, validate0))  
    
    ### DONE ###
    
  }
  
  result <- as.data.frame(test_rmse)
  result$new_feature <- "test"
  result$rep <- "test"
  result$subset_size <- c(1:12)
  result$type <- "validation"
  
  
#################################################################################################    

# perform 10 repetitions of random feature selections
# to compare to the result of rfe
  
out <- list()
  
for (j in 1:10){
  
  print(paste("repetition ", j, sep = ""))
  
  #list according featurs
  feats <- feats0 %>%
    filter(year == 2016, trait == "GY") %>%
    select(subset_size)
  
  #select 1st feature at random
  feats <- as.character(unlist(feats))
  feats_sel <- sample(feats, 1)
  
  #prepare vectors to take up results
  test_rmse <- NULL
  new_feature <- NULL
  
  for (i in 1:length(subsets)) {
    
    print(paste("==> subset size: ", i, sep = ""))
    
    if(i > 1){
      #randomly select the next feature, excluding the already selected features
      next_feat <- sample(feats[!feats %in% feats_sel], 1)
      #add selected feature to the subset
      feats_sel <- c(feats_sel, next_feat)
    } else NULL
    
    #define the training data
    train <- train0 %>%
      select(one_of(feats_sel), gy) 
    train <- train[complete.cases(train),]
    
    ### TUNE/TRAIN RANDOM FOREST USING ALL PREDICTORS; CALCULATE VIP AND TEST RMSE ###
    
    #adjust mtry parameter to decreasing predictor set
    mtry <- c(ceiling((length(train)-1)/10), ceiling((length(train)-1)/5), 
              ceiling((length(train)-1)/3), ceiling((length(train)-1)/1.5))
    
    #specify model tuning parameters
    tune_grid <- expand.grid(mtry = mtry,
                             splitrule = "variance",
                             min.node.size = 5)
    
    #define resampling procedure
    ctrl <- trainControl(method = "repeatedcv",
                         number = 10,
                         rep = 1,
                         verbose = FALSE,
                         allowParallel = TRUE)
    
    #random forest
    rf_ranger <- train(gy~.,
                       data = train,
                       preProc = c("center", "scale"),
                       method = "ranger",
                       tuneGrid = tune_grid,
                       importance = "permutation",
                       num.trees = 1000,
                       num.threads = 16,
                       trControl = ctrl)
    
    test_rmse[i] <- rmse(validate0$gy, predict(rf_ranger, validate0)) 
    new_feature[i] <- ifelse(i == 1, feats_sel, next_feat)
    
    ### DONE ###
    
  }
  
  out[[j]] <- cbind(test_rmse, new_feature)
  
}
  
saveRDS(out, "O:/Projects/KP0011/1/Analysis/RFE/validation/random_subsets.rds")
  
tidy <- as.data.frame(do.call(rbind, out))
tidy$subset_size <- rep(1:12, 10)
tidy$rep <- rep(1:10, each = 12)
tidy$test_rmse <- as.numeric(as.character(tidy$test_rmse))
tidy$type <- "random"

all_data <- rbind(tidy, result)

ggplot(all_data, aes(x = subset_size, y = test_rmse, group = rep, color = type)) +
  geom_line() +
  geom_line() +
  theme_bw()
  
  
  
  