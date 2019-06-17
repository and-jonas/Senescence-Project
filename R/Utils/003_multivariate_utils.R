
#calculate RMSE
get_rmse = function(actual, predicted) {
  sqrt(mean((actual - predicted) ^ 2))
}

#calculate RMSE on upsampled validation dataset
get_RMSE_upsample <- function(data){
  data[, Trait] <- as.factor(data[, Trait])
  data <- upSample(x = data,
                   y = data[, Trait],
                   yname = Trait)
  data[, Trait] <- as.numeric(as.character(data[, Trait]))
  data <- data[, -ncol(data)]
  #count samples used for validation
  nsamp <- nrow(data)
  ID_sub <- data %>% 
    dplyr::select(-contains("rflt"), -SnsCnp, -SnsFl0)
  #remove IDs from validation data
  data_sub <- data %>% 
    select(contains("rflt"), SnsCnp, -SnsFl0)
  #predict new data using trained model
  pred <- caret::predict.train(fit, newdata = data_sub, 
                               type = "raw", na.action = na.omit)
  predobs <- cbind(ID_sub, data[, Trait], pred) %>% as_tibble()
  names(predobs)[length(predobs)-1] <- "obs"
  RMSE <- round(rmse(predobs$obs, predobs$pred), 2)
  return(RMSE)
}

#evaluate model performance within and across years
eval_multivariate_cross <- function(Trait, method = "pls",
                                    maxcomp = NULL, tuneGrid = NULL, 
                                    data, train, test, 
                                    trainsample, 
                                    preProc, 
                                    bdir){
  
  #train:       A list of training datasets. Each list element consists of a vector with the year(s) used for training
  #test:        A list of validation sets corresponding to ONE training set; 
  #             Each list element consists of a list of vectors; 
  #             Vectors contain all possible combinations of validation data sets. 
  #data:        A list of dataframes, each list element contains all the data of a particular type. 
  #trainsample: Procedure to assure class balance in the training set. Either "downsampling" or "fullsample".
  #preProcess:  Preprocessing of spectral data: Either "none", "center", "scale", or "center_scale"
  #bdir:        Sink directory for output
  
  #for each of the possible training datasets
  #i.e. for each year and all combinations of years
  out_train <- out_all <- list()
  for (i in 1:length(train)) {
    
    #Verbose
    print(paste("> train on", paste(train[[i]], collapse = "_")))
    
    #====================================================================================== -
    
    # SELECT TRAINING DATA
    
    #extract data for training year(s)
    dat_train0 <- data %>% filter(Exp %in% train[[i]]) %>% as.data.frame()
    
    #select training data:
    ##sample data if required
    if(trainsample == "downsample"){
      dat_train <- dat_train0
      dat_train[, Trait] <- dat_train[, Trait] %>% as.factor()
      dat_train <- downSample(x = dat_train,
                              y = dat_train[, Trait],
                              yname = Trait)
      dat_train[, Trait] <- as.numeric(as.character(dat_train[, Trait]))
      dat_train <- dat_train[,-ncol(dat_train)]
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
    
    #====================================================================================== -
    
    # TRAIN ALGORITHM
    
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
      tuneGrid <- expand.grid(committees = c(1, 2, 5, 10, 20),
                              neighbors = 0)
      fit <- caret::train(as.formula(paste(Trait, "~.", sep = "")),
                          data = dat, 
                          preProcess = preProcess,
                          method = method,
                          tuneGrid = tuneGrid, 
                          trControl = ctrl,
                          importance = TRUE)
    }
    
    #====================================================================================== -
    
    # PREDOBS FOR TRAINING SAMPLES
    
    predobs <- plyr::match_df(fit$pred, fit$bestTune, on = names(fit$bestTune)) %>% 
      group_by(rowIndex) %>%
      summarize(obs = mean(obs),
                mean_pred = mean(pred))
    
    #add observation ID
    ID_train$rowIndex <- 1:nrow(ID_train)
    predobs <- dplyr::full_join(ID_train, predobs, by = "rowIndex") %>% 
      dplyr::select(-rowIndex) %>% as_tibble()
    
    RMSE_train <- caret::getTrainPerf(fit) %>% pull(TrainRMSE)
    RMSE_train_up <- replicate(10, get_RMSE_upsample(dat_train0)) %>% mean()
    
    #====================================================================================== -
    
    # PREDOBS FOR ALL SAMPLES OF TRAINING EXP
    
    #extract observation IDs
    ID_train0 <- dat_train0 %>% 
      select(-contains("rflt_"), -SnsCnp, -SnsFl0)
    
    #predict new data using trained model
    pred <- caret::predict.train(fit, newdata = dat_train0, 
                                 type = "raw", na.action = na.omit)
    
    predobs_all <- cbind(ID_train0, dat_train0[, Trait], pred) %>% as_tibble()
    names(predobs_all)[length(predobs_all)-1] <- "obs"
    
    #====================================================================================== -
    
    # CREATE OUTPUT OF TRAINING
    
    out_train[[i]] <-  tibble(Trait = Trait, 
                              train = paste(train[[i]], collapse = "_"),
                              datatype = k,
                              trainsample = trainsample,
                              method = method, 
                              tunepar = list(fit$bestTune[1, 1]), 
                              preProc = preProc,
                              RMSE_train = RMSE_train,
                              RMSE_train_up = RMSE_train_up,
                              nsamp_train = nsamp,
                              predobs_train_all = list(predobs_all)
                              )
    
    #====================================================================================== -
    #TRAINING AND VALIDATION ON THE SAME YEAR COMPLETE#
    #====================================================================================== -
    
    # Model validation on validation data: different Experiment(s)
    
    #create empty list to take up predobs data
    out_test <- list()
    
    #for each validation dataset
    for(j in 1:length(test[[i]])){
      
      #extractcorresponding data from the complete dataset
      dat_test0 <- data %>% 
        filter(Exp %in% test[[i]][[j]]) %>% as.data.frame()
      
      #count samples used for validation
      nsamp_all <- nrow(dat_test0)
      
      ID_test_sub <- dat_test0 %>% 
        select(-contains("rflt"), -SnsCnp, -SnsFl0)
      
      #remove IDs from validation data
      dat_test_sub <- dat_test0 %>% 
        select(contains("rflt"), SnsCnp, -SnsFl0)
      
      #predict new data using trained model
      pred <- caret::predict.train(fit, newdata = dat_test_sub, 
                                   type = "raw", na.action = na.omit)
      
      predobs <- cbind(ID_test_sub, dat_test0[, Trait], pred) %>% as_tibble()
      names(predobs)[length(predobs)-1] <- "obs"
      
      RMSE_test_up <- replicate(10, get_RMSE_upsample(dat_test0)) %>% mean()
      
      #######################################################################
      #extract predictions for all observations in validation year
      #######################################################################
      
      #extract observation IDs
      ID_test0 <- dat_test0 %>% 
        select(-contains("rflt"), -SnsCnp, -SnsFl0)
      
      #predict new data using trained model
      pred <- caret::predict.train(fit, newdata = dat_test0, 
                                   type = "raw", na.action = na.omit)
      
      predobs_all <- cbind(ID_test0, dat_test0[, Trait], pred) %>% as_tibble()
      names(predobs_all)[length(predobs_all)-1] <- "obs"
      
      RMSE_all <- round(rmse(predobs_all$obs, predobs_all$pred), 2)
  
      #create a list containing all the output for the training dataset and the internal validation
      
      out_test[[j]] <-  tibble(test = paste(test[[i]][[j]], collapse = "_"),
                               RMSE_test = RMSE_all,
                               RMSE_test_up = RMSE_test_up,
                               nsamp_test = nsamp,
                               predobs_test_all = list(predobs_all)
                               )
                               
    } #end of validation sets loop
    
    #combine the training and validation results
    out_test_done <- bind_rows(out_test)
    out_all[[i]] <- tibble(train = list(out_train[[i]]), test = list(out_test_done)) %>% 
      unnest(train) %>% unnest(test, .preserve = c(predobs_train_all, tunepar))
  
  }#end of training set loop
  
  out <- bind_rows(out_all)
    
} #end of function
