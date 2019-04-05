#====================================================================================== -

#HEADER ----

# AUTHOR: Jonas Anderegg

# Train cubist model to predict senescence scorings
# Perform recurisve feature elimination with cubist as base learner
# to identify the most important wavelengths

#====================================================================================== -

.libPaths("T:/R3UserLibs")

library(caret)

dir <- "O:/Projects/KP0011/1/"
setwd(dir)

#get data for each of the three years
train0 <- readRDS("Analysis/rfe_cubist/train0_fpww012.rds")
train0 <- readRDS("Analysis/rfe_cubist/train0_fpww018.rds")
train0 <- readRDS("Analysis/rfe_cubist/train0_fpww022.rds")

# create the cross-validation files as a list to use with different 
# functions
index <- caret::createMultiFolds(train0$SnsCnp, k = 5, times = 6)

# function to calculate RMSE  
rmse = function(actual, predicted) {
  sqrt(mean((actual - predicted) ^ 2))
}

#################################################################################################  
#################################################################################################

## PERFORM RERUSIVE FEATURE ELIMINATION 
## EMBEDDED IN AN OUTER RESAMPLING LOOP

# SET THRESHOLD FOR VARIMP

  th0 <- 10
  th <- 10

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
    train_rmse <- NULL
    test_rmse <- NULL
    varimp_all <- NULL
    
### TUNE/TRAIN CUBIST USING ALL PREDICTORS; CALCULATE VIP AND TEST RMSE ###
    
    #define tuning parameter grid
    train_grid <-   train_grid <- expand.grid(committees = c(1, 2, 5, 10, 20, 40),
                                              neighbors = 0)
    
    #define inner resampling procedure
    tr_ctrl <- trainControl(method = "repeatedcv", 
                            number = 10,
                            repeats = 1,
                            savePredictions = TRUE, 
                            selectionFunction = "oneSE",
                            verboseIter = FALSE,
                            allowParallel = TRUE)
  
    #set up a cluster
    library(parallel)
    library(doParallel)
    cluster <- makeCluster(11, outfile = "")
    registerDoParallel(cluster)
    clusterEvalQ(cluster, library(doParallel))
    clusterEvalQ(cluster, .libPaths("T:/R3UserLibs"))
    clusterEvalQ(cluster, library(doParallel))
    
    #tune/train a cubist regression model
    train_cubist <- train(
      SnsCnp ~ .,
      data = train,
      preProcess = c("center", "scale"),
      method = "cubist",
      tuneGrid = train_grid,
      trControl = tr_ctrl
    )
    
    #get variable importance
    vip <- varImp(train_cubist)$importance
    vip$vars <- rownames(vip)
    vip <- vip[order(vip$Overall, decreasing = TRUE), ]
    
    #assign output to vectors
    varimp_all[[1]] <- vip
    keep_vars[[1]] <- vip[vip$Overall > th0, "vars"] 
    train_rmse[1] <- getTrainPerf(train_cubist)[,"TrainRMSE"]
    test_rmse[1] <- rmse(test$SnsCnp, predict(train_cubist, test))  

    ### DONE ###
    
    #reduce train data: 
    #exclude the least important features
    #according to varImp output
    new_vars <- c(vip[vip$Overall > th0, "vars"], "SnsCnp")
    newtrain <- select(train, new_vars)
    
### FEATURE ELIMINATION LOOP ###
    
    #perform backwards feature selection using RFE with cubist as base-learner
    #number of committees are tuned using a constant grid,
    #instance-based correction is not done
    #sequentially remove predictors with a varImp value of <= 10
    #until only 1 predictors remain
    while(length(newtrain) > 50){
      
      print(paste("numb. remaining predictors = ", length(newtrain)-1, sep = ""))
      
      #train
      train_cubist <- train(
        SnsCnp ~ .,
        data = newtrain,
        preProcess = preProcess,
        method = "cubist",
        tuneGrid = train_grid,
        trControl = tr_ctrl
      )
      
      #get variable importance
      vip <- varImp(train_cubist)$importance
      vip$vars <- rownames(vip)
      vip <- vip[order(vip$Overall, decreasing = TRUE), ]

      #assign output to vector
      varimp_all <- c(varimp_all, list(vip))
      keep_vars <- c(keep_vars, list(vip[vip$Overall > th, "vars"]))
      
      if(nrow(vip) > 1){
        new_vars <- c(vip[vip$Overall > th, "vars"], "SnsCnp")
      } else {
        new_vars <- "SnsCnp"
      }
      
      train_rmse <- c(train_rmse, getTrainPerf(train_cubist)[,"TrainRMSE"])
      test_rmse <- c(test_rmse, rmse(test$SnsCnp, predict(train_cubist, test)))
  
      #reduce train data: 
      #exclude predefined number of most irrelevant predictors 
      #according to sizes sequence
      newtrain <- select(newtrain, new_vars)
  
    }#END OF FEATURE ELIMINATION LOOP #1
    
    while(length(newtrain) > 1){
      
      print(paste("numb. remaining predictors = ", length(newtrain)-1, sep = ""))
      
      #train
      train_cubist <- train(
        SnsCnp ~ .,
        data = newtrain,
        preProcess = preProcess,
        method = "cubist",
        tuneGrid = train_grid,
        trControl = tr_ctrl
      )
      
      #get variable importance
      vip <- varImp(train_cubist)$importance
      vip$vars <- rownames(vip)
      vip <- vip[order(vip$Overall, decreasing = TRUE), ]
      
      #assign output to vector
      varimp_all <- c(varimp_all, list(vip))
      keep_vars <- c(keep_vars, list(vip[vip$Overall > 1, "vars"]))
      
      if(nrow(vip) > 1){
        new_vars <- c(vip[vip$Overall > th, "vars"], "SnsCnp")
      } else {
        new_vars <- "SnsCnp"
      }
      
      train_rmse <- c(train_rmse, getTrainPerf(train_cubist)[,"TrainRMSE"])
      test_rmse <- c(test_rmse, rmse(test$SnsCnp, predict(train_cubist, test)))
      
      #reduce train data: 
      #exclude predefined number of most irrelevant predictors 
      #according to sizes sequence
      newtrain <- select(newtrain, new_vars)
      
    }#END OF FEATURE ELIMINATION LOOP #2
    
### ASSIGN RANKS TO VARIABLES DEPENDING ON THE STEP OF ELIMINATION ###
    
    #assign ranks to predictors based on the variable importance 
    varranks_all <- list()
    for(k in 1:length(varimp_all)){
      # exception for the last iteration...
      if(k !=length(varimp_all)){
        var <- varimp_all[[k]][varimp_all[[k]]$Overall <= th, "vars"]
        rank <- rep(length(varimp_all)-k+1, length(var))
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
    RMSEs <- tibble::as.tibble(cbind(sapply(varimp_all, nrow), train_rmse, test_rmse)) %>%
      rename(subset_size = V1)
  
    #store in list
    saveObj <- list(varranks, RMSEs)
    saveRDS(saveObj, paste("Analysis/rfe_cubist/Resamples/resample_", i, ".rds", sep = ""))
    out[[i]] <- list(varranks, RMSEs)
    
  }#END OF OUTER RESAMPLING
  
  saveRDS(out, "Analysis/rfe_cubist/Resamples/resample_all.rds")
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
  
  stopCluster(cluster)
  registerDoSEQ()
  