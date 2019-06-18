
#Helper functions
########################################################################################

#calculate RMSE
get_rmse = function(actual, predicted) {
  sqrt(mean((actual - predicted) ^ 2))
}

#perform up- or downsampling
perform_sampling <- function(data, method = "upsampling"){
  data[, Trait] <- as.factor(data[, Trait])
  if(method == "upsampling"){
    data <- upSample(x = data,
                     y = data[, Trait],
                     yname = Trait)
  } else if (method == "downsample"){
    data <- downSample(x = data,
                       y = data[, Trait],
                       yname = Trait)
  } else print("specify sampling procedure for test data")
  data[, Trait] <- as.numeric(as.character(data[, Trait]))
  data <- data[, -ncol(data)]
  return(data)
}

#calculate RMSE on upsampled validation dataset
get_RMSE_upsample <- function(data){
  data <- perform_sampling(data, method = "upsampling")
  ID_sub <- data %>% 
    dplyr::select(-contains("rflt"), -SnsCnp, -SnsFl0)
  data_sub <- data %>% 
    select(contains("rflt"), SnsCnp, -SnsFl0)
  #predict new data using trained model
  pred <- caret::predict.train(fit, newdata = data_sub, 
                               type = "raw", na.action = na.omit)
  predobs <- cbind(ID_sub, data[, Trait], pred) %>% as_tibble()
  names(predobs)[length(predobs)-1] <- "obs"
  RMSE <- round(get_rmse(predobs$obs, predobs$pred), 2)
  return(RMSE)
}

#generate predictions and extract corresponding real observations
extract_predobs <- function(data){
  ID <- data %>% 
    select(-contains("rflt"), -SnsCnp, -SnsFl0)
  dat <- data %>% 
    select(contains("rflt"), SnsCnp, -SnsFl0)
  #create predictions
  pred <- caret::predict.train(fit, newdata = dat, 
                               type = "raw", na.action = na.omit)
  #create output tibble
  predobs <- cbind(ID, dat[, Trait], pred) %>% as_tibble()
  names(predobs)[length(predobs)-1] <- "obs"
  return(predobs)
}


#Main function
########################################################################################

#evaluate model performance within and across years
eval_multivariate_cross <- function(Trait, method = "pls",
                                    maxcomp = NULL, tuneGrid = NULL, 
                                    data, train, test, 
                                    trainsample = "fullsample", 
                                    testsample = "upsample",
                                    preProc = "center_scale", 
                                    runParallel = FALSE){
  
  #set up a cluster
  if(runParallel){
    cluster <- makeCluster(12, outfile = "")
    registerDoParallel(cluster)
    clusterEvalQ(cluster, library(MASS))
    clusterEvalQ(cluster, library(caret))
    clusterEvalQ(cluster, library(tidyverse))
  }
  
  #for each of the possible training datasets
  #i.e. for each year and all combinations of years
  out_train <- out_all <- list()
  for (i in 1:length(train)) {
    
    #Verbose
    print(paste("> train on", paste(train[[i]], collapse = "_")))
    
    #====================================================================================== -
    #TRAIN ALGORITHM, VALIDATE ON THE SAME YEAR#
    #====================================================================================== -
    
    #extract data for training year(s)
    dat_train0 <- data %>% filter(Exp %in% train[[i]]) %>% as.data.frame()
    
    #sample data if required
    if(trainsample == "downsample"){
      dat_train <- perform_sampling(dat_train0, method = "downsampling")
    } else if (trainsample == "upsample"){
      dat_train <- perform_sampling(dat_train0, method = "upsampling")
    } else if(trainsample == "fullsample") {
      dat_train <- dat_train0
    } else print("specify sampling procedure for train data")
    
    #count samples used for training
    nsamp <- nrow(dat_train)
    
    #extract IDs of samples used for training
    ID_train <- dat_train %>% 
      select(-contains("rflt"), -SnsCnp, -SnsFl0)
    
    #keep only predictors
    dat <- dat_train %>% 
      select(Trait, contains("rflt"))
    
    if(preProc == "none"){
      preProcess <- NULL
    } else if (preProc == "center_scale") {
      preProcess <- c("center", "scale")
    } else if (preProc == "scale"){
      preProcess <- "scale"
    } else if (preProc == "center") {
      preProcess <- "center"
    } else print("preProcessing needs to be specified!")
    
    train_index <- caret::createMultiFolds(y = dat[, Trait], k = 10, times = 1)
    ctrl <- caret::trainControl(method = "repeatedcv", 
                                index = train_index,
                                savePredictions = TRUE, 
                                selectionFunction = "oneSE",
                                allowParallel = TRUE,
                                verboseIter = TRUE)
    if(method == "pls"){
      fit <- caret::train(as.formula(paste(Trait, "~.", sep = "")),
                          data = dat, 
                          preProcess = preProcess,
                          method = method,
                          tuneLength = maxcomp, 
                          trControl = ctrl,
                          importance = TRUE)
    } else if(method == "cubist"){
      fit <- caret::train(as.formula(paste(Trait, "~.", sep = "")),
                          data = dat, 
                          preProcess = preProcess,
                          method = method,
                          tuneGrid = tuneGrid, 
                          trControl = ctrl,
                          importance = TRUE)
    }

    #====================================================================================== -

    #extract predobs and train RMSE
    predobs <- plyr::match_df(fit$pred, fit$bestTune, on = names(fit$bestTune)) %>% 
      group_by(rowIndex) %>%
      summarize(obs = mean(obs),
                mean_pred = mean(pred))
    ID_train$rowIndex <- 1:nrow(ID_train)
    predobs <- dplyr::full_join(ID_train, predobs, by = "rowIndex") %>% 
      dplyr::select(-rowIndex) %>% as_tibble()
    RMSE_train <- caret::getTrainPerf(fit) %>% pull(TrainRMSE)
    
    #calcualte RMSE for upsamples
    if(isTRUE(runParallel)){
      clusterExport(cluster, c("dat_train0", "get_RMSE_upsample", "perform_sampling", 
                               "get_rmse", "Trait", "fit"), 
                    envir=environment())
      RMSE_train_up <- parSapply(cluster, 1:10, function(i, ...) {x <- get_RMSE_upsample(dat_train0)}) %>% 
        mean() %>% round(., 2)
    } else {
      RMSE_train_up <- replicate(10, get_RMSE_upsample(dat_train0)) %>% mean()
    }
    print("==train RMSE calcualted==")

    #create predictions for all samples in the training dataset if required;
    #extract corresponding real observations
    if(trainsample == "downsample"){
      predobs_all <- extract_predobs(dat_train0)
    } else {
      predobs_all <- predobs
    }
    
    #create output tibble
    out_train[[i]] <-  tibble(Trait = Trait, 
                              train = paste(train[[i]], collapse = "_"),
                              trainsample = trainsample,
                              method = method, 
                              tunepar = list(fit$bestTune[1, 1]), 
                              preProc = preProc,
                              RMSE_train = RMSE_train,
                              RMSE_train_up = RMSE_train_up,
                              nsamp_train = nsamp,
                              predobs_train_all = list(predobs_all))
    
    #====================================================================================== -
    #VALIDATION ON DIFFERENT EXPERIMENT(S)#
    #====================================================================================== -
    
    out_test <- list()
    #for each validation dataset
    for(j in 1:length(test[[i]])){
      
      #get data
      dat_test0 <- data %>% 
        filter(Exp %in% test[[i]][[j]]) %>% as.data.frame()
      
      #count samples used for validation
      nsamp_all <- nrow(dat_test0)
      
      #create predictions,
      #extract corresponding real observations
      predobs <- extract_predobs(dat_test0)
      
      #calculate RMSE on upsamples
      #or on fullsamples
      if(isTRUE(runParallel)){
        clusterExport(cluster, c("dat_test0"), envir=environment())
        RMSE_test_up <- parSapply(cluster, 1:10, function(i, ...) {x <- get_RMSE_upsample(dat_test0)}) %>%
          mean() %>% round(., 2)
      } else {
        RMSE_test_up <- replicate(10, get_RMSE_upsample(dat_test0)) %>% 
        mean() %>% round(., 2)
      }
      predobs_all <- extract_predobs(dat_test0)
      RMSE_test <- round(get_rmse(predobs_all$obs, predobs_all$pred), 2)
      print(paste("==test RMSE calculated,", "(", j, ")", "==", sep = ""))
      
      #create output tibble
      out_test[[j]] <-  tibble(test = paste(test[[i]][[j]], collapse = "_"),
                               RMSE_test = RMSE_test,
                               RMSE_test_up = RMSE_test_up,
                               nsamp_test = nsamp,
                               predobs_test_all = list(predobs_all))
    }#end of validation loop
    out_test_done <- bind_rows(out_test)
    out_all[[i]] <- tibble(train = list(out_train[[i]]), test = list(out_test_done)) %>% 
      unnest(train) %>% unnest(test, .preserve = c(predobs_train_all, tunepar))
  }#end of training set loop
  out <- bind_rows(out_all)
  if(isTRUE(runParallel)){
    stopCluster(cluster)
    registerDoSEQ()
  }
  return(out)
}#end of function


#train:       A list of training datasets; 
#             Each list element consists of a vector with the Exp(s) used for training
#test:        A list of validation sets corresponding to ONE training set; 
#             Each list element consists of a list of vectors; 
#             Vectors contain all possible combinations of validation data sets. 
#trainsample: Procedure to assure class balance in the training set if required; 
#             Either "downsampling" or "fullsample".
#             It is strongly recommended to use "fullsample";
#             Downsampling negatively affects model performance; 
#             Upsampling results in a large increase in computation time. 
#preProcess:  Preprocessing of spectral data;
#             Either "none", "center", "scale", or "center_scale".

