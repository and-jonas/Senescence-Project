
.libPaths("T:/R3UserLibs")

library(tidyverse)
library(readxl)
library(stringi)
library(stringr)
library(data.table)
library(caret)


#Prepare data
#################################################################

#define directories
adir <- "O:/Projects/KP0011/1/Analysis/multvarmod/SpectralData/"
bdir <- "O:/Projects/KP0011/1/Analysis/multvarmod/Results/PLSR/"
setwd(adir)

#list data type files
files <- as.list(list.files(adir))[2]
#read all data into list (all tpyes of spectral data)
dat1 <- lapply(files, readRDS)
#rename list elements
names(dat1) <- sapply(strsplit(unlist(files), "\\."), "[[", 1)
#reorder columns
dat1 <- lapply(dat1, function(x) x %>% 
                 dplyr::select(-contains("rflt"), everything()))

subset <- paste("rflt_", c(677, 686, 695, 710, 725, 755, 767), sep = "")
datsub <- lapply(dat1, function(x) x %>% dplyr::select(1:19, one_of(subset)))

# define train and validation sets
#################################################################

# Experiments
  Exp <- c("FPWW012", "FPWW018", "FPWW022")
    
# create train list of length 7
  train <- list(Exp[1], Exp[2], Exp[3], Exp[1:2], Exp[2:3], Exp[c(1, 3)], Exp[1:3])
    
# create validation list of length 7
  valid <- list(c(list(Exp[2], Exp[3], Exp[2:3])), 
                c(list(Exp[1], Exp[3], Exp[c(1, 3)])),
                c(list(Exp[1], Exp[2], Exp[1:2])),
                c(list(Exp[3])),
                c(list(Exp[1])),
                c(list(Exp[2])),
                c(list(Exp[1], Exp[2], Exp[3])))

###########################################################################################

#Example Run
  
library(parallel)
library(doParallel)
cluster <- makeCluster(12, outfile = "") # convention to leave 1 core for OS
registerDoParallel(cluster)
clusterEvalQ(cluster, library(doParallel))
clusterEvalQ(cluster, .libPaths("T:/R3UserLibs"))
clusterEvalQ(cluster, library(doParallel))

start.time <- Sys.time()
  
eval_pls_cross(train = train, 
               valid = valid, 
               data = dat1, 
               sampling_train = "fullsample",
               sampling_val = "downsample",
               preProc = "none",
               maxcomp = 12,
               bdir = bdir) 

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

stopCluster(cluster)
registerDoSEQ()

###########################################################################################
  
# create predobs plots: for one directory

  preProc <- sapply(strsplit(unlist(dirs), "\\/"), "[[", 2)

  dir <- "O:/FIP/2017/WW018/ASD/Renamed/RData/Results/PLSR/fullsample/center_scale"
  
  preProc <- "center_scale"
  
  files <- as.list(list.files(paste(dir), pattern = ".rds"))

  setwd(dir)

  data <- lapply(files, readRDS)
  names(data) <- sapply(strsplit(unlist(files), "\\."), "[[", 1)
  
  for(i in names(data)){
    
    d <- data[[i]]
    
    D <- unlist(d, recursive = FALSE)
    
    plot_predobs <- function(data) {
      
      data$data$obsf <- as.factor(data$data$obs)
      
      Plot <- ggplot() +
        geom_boxplot(aes(x = data$data$obsf, y = data$data$pred),
                     size = 0.5, alpha = 1) +
        stat_summary(mapping = aes(x = data$data$obs + 1, y = data$data$pred), fun.y = "median", colour = "red", size = 3, geom = "point", shape=18) +
        scale_y_continuous(breaks = seq(0,10,2), limits = c(-2, 12)) +
        geom_abline(intercept = -1, slope = 1, color = "red", size = 1, linetype="dashed") +
        xlab("observed") + ylab("predicted (5-times repeated 10-fold CV)") +
        annotate("text", x = 9, y = 0, label= paste("RMSE = ", data$RMSE, sep = ""), size = 2.5) +
        annotate("text", x = 9, y = -0.5, label = paste("ncomp = ", data$ncomp, sep = ""), size = 2.5) +
        annotate("text", x = 9, y = -1, label = paste("nsamp = ", data$nsamp, sep = ""), size = 2.5) +
        ggtitle(paste("training = ", paste(data$training, collapse = ", "),
                      "\nvalidation = ", paste(data$validation, collapse = ", "),
                      paste("\ndata = ", i, sep = ""), paste("\npreProc = ", preProc, sep = ""), sep = "")) +
        theme_bw() +
        theme(plot.title = element_text(size = 7, face = "bold"))
      
      return(Plot)
      
    }
    
    Plots <- lapply(D, plot_predobs)
    
    library(gridExtra)
    # arrange them in a 2x2 grid
    
    ml <- marrangeGrob(Plots, nrow=2, ncol=2)
    
    ggsave(paste("O:/FIP/2017/WW018/ASD/Renamed/RData/Results/PLSR/fullsample/", 
                 preProc, "/", i, ".pdf", sep = ""), 
           ml, 
           width = 297, 
           height = 210, 
           units = "mm")
    
  }
  
###########################################################################################
  
# Graph for publication
  
  data <- readRDS("O:/FIP/2017/WW018/ASD/Renamed/RData/Results/PLSR/downsample/center_scale/smth_avg_rflt.rds")
  example <- data[[6]]
  
  preProc <- "center_scale"
  
  Plots <- list()
  
  for (i in 1:length(example)){
    
    d <- example[[i]]
    
    Plots[[i]] <- plot_predobs(d)

  }
  
  ml <- marrangeGrob(Plots, nrow=1, ncol=2)
  
  
###########################################################################################
  
# create predobs plots: Loop through all directories
  
  dirs <- as.list(list.files(bdir, recursive = TRUE, pattern = ".rds"))
  dirs1 <- unique(sub("(/[^/]+)/.*", "\\1", dirs))
  
  for(j in dirs1){
    
    preProc <- unlist(strsplit(j, "\\/"))[2]
    sampling <- unlist(strsplit(j, "\\/"))[1]
    
    dir <- paste("O:/FIP/2017/WW018/ASD/Renamed/RData/Results/PLSR", sampling, preProc, sep = "/")
    
    files <- as.list(list.files(paste(dir), pattern = ".rds"))
    
    setwd(dir)
    
    data <- lapply(files, readRDS)
    names(data) <- sapply(strsplit(unlist(files), "\\."), "[[", 1)
    
    for(i in names(data)){
      
      d <- data[[i]]
      
      D <- unlist(d, recursive = FALSE)
      
      plot_predobs <- function(data) {
        
        data$data$obsf <- as.factor(data$data$obs)
        
        Plot <- ggplot() +
          geom_boxplot(aes(x = data$data$obsf, y = data$data$pred),
                       size = 0.5, alpha = 1) +
          stat_summary(mapping = aes(x = data$data$obs + 1, y = data$data$pred), fun.y = "median", colour = "red", size = 3, geom = "point", shape=18) +
          scale_y_continuous(breaks = seq(0,10,2), limits = c(-2, 12)) +
          geom_abline(intercept = -1, slope = 1, color = "red", size = 1, linetype="dashed") +
          xlab("observed") + ylab("predicted (5-times repeated 10-fold CV)") +
          annotate("text", x = 9, y = 0, label= paste("RMSE = ", data$RMSE, sep = ""), size = 2.5) +
          annotate("text", x = 9, y = -0.5, label = paste("ncomp = ", data$ncomp, sep = ""), size = 2.5) +
          annotate("text", x = 9, y = -1, label = paste("nsamp = ", data$nsamp, sep = ""), size = 2.5) +
          ggtitle(paste("training = ", paste(data$training, collapse = ", "),
                        "\nvalidation = ", paste(data$validation, collapse = ", "),
                        paste("\ndata = ", i, sep = ""), paste("\npreProc = ", preProc, sep = ""), sep = "")) +
          theme_bw() +
          theme(plot.title = element_text(size = 7, face = "bold"))
        
        return(Plot)
        
      }
      
      Plots <- lapply(D, plot_predobs)
      
      library(gridExtra)
      # arrange them in a 2x2 grid
      
      ml <- marrangeGrob(Plots, nrow=2, ncol=2)
      
      ggsave(paste("O:/FIP/2017/WW018/ASD/Renamed/RData/Results/PLSR/", sampling, "/", preProc, "/", i, ".pdf", sep = ""), ml)
      
    }
    
  }
  
###########################################################################################

# further analysis steps
  ## for all types of data and all pre-treatments,
  ## extract model information and obtained predictions

  # load data
    basedir <- "O:/Projects/KP0011/1/Analysis/multvarmod/Results/PLSR/SnsCnp"  
    fulldirs <- as.list(list.files(basedir, pattern = ".rds", recursive = TRUE, full.names = TRUE))
    data_all <- lapply(fulldirs, readRDS)
  
  # function to extract model information
  
    ex_inf <- function(data) {
      
      RMSEs <- sapply(data, function(x) sapply(x, '[[', 5)) %>% Reduce(c, .)
      trainexp0 <- lapply(data, function(x) lapply(x, '[[', 1))
      trainexp <-  as.vector(sapply(trainexp0, function(x) sapply(x, paste, collapse = "_"))) %>% Reduce(c, .)
      validexp0 <- lapply(data, function(x) lapply(x, '[[', 2))
      validexp <- as.vector(sapply(validexp0, function(x) sapply(x, paste, collapse = "_"))) %>% Reduce(c, .)
      preProc <- as.vector(sapply(data, function(x) sapply(x, '[[', 8))) %>% Reduce(c, .)
      datatype <- as.vector(sapply(data, function(x) sapply(x, '[[', 3))) %>% Reduce(c, .)
      sampling <- as.vector(sapply(data, function(x) sapply(x, '[[', 6))) %>% Reduce(c, .)
      inf <- cbind(sampling, datatype, preProc, trainexp, validexp, RMSEs)
    
      return(inf)
      
    }
    
  # extract model information
  
    all_inf <- lapply(data_all, ex_inf)
  
  # extract model predictions

    preds <- lapply(data_all, function(y) lapply(y, function(x) lapply(x, '[[', 10)) %>% Reduce(c, .))

###########################################################################################

# wrapper to evaluate scaled predictions obtained from all models 
    # Generate correlation information, mean error

###########################################################################################

# Load scoring data
  sen_d <- readRDS("T:/PhD/DataAnalysis/FPWW012_FPWW018/Data/sen_d_clean.rds")

# define parameters to be evaluated
  ind <- "pred"

# wrapper function
  
  evaluate_plspreds <- function(data){
    
    #this function uses several costum functions!!
    
    #remove obs values from dataframe
    eval_data <- data %>% 
      select(-obs) %>%
      #scale predictions to range [0:10]
      f_scale_pred() %>%
      #select required variables
      select(Plot_ID, spc_ID, meas_date, pred) %>%
      #join to scaled scorings
      f_match_join(sen_d) %>%
      #reorder columns as required by main evaluation function
      select(heading_DAS, heading_GDDAS, grading_DAS, grading_DAH, grading_GDDAS, 1:12, SnsFl0, SnsCnp, everything()) %>%
      #evaluate predictions using custom function
      f_calc_err(method = "lin", ind = ind, traits = c("SnsCnp"))
      #extract correlation coefficients for dynamics parameters and total error
      means2 <- extract_inf(eval_data, trait = "SnsCnp", out = "means_all")
      stat_data1 <- extract_inf(eval_data, trait = "SnsCnp", out = "stat_data_all")
      
      out <- list(means2, stat_data1)
      
    return(out)
    
  }

# use the wrapper to extract all the information for all data types
  
  OUT <- lapply(preds, function(x) lapply(x,  evaluate_plspreds))
  
    saveRDS(OUT, "O:/Projects/KP0011/1/Analysis/multvarmod/Results/PLSR/Final_Result/SnsCnp/result_pars.rds")
    
  
###########################################################################################
    
  OUT <- readRDS("O:/Projects/KP0011/1/Analysis/multvarmod/Results/PLSR/Final_Result/SnsCnp/result_pars.rds")

# Rearrange output    
    
  tidyresults <- list()
  
  for (i in 1:length(OUT)){
    
    #get result each data type
    result <- OUT[[i]]
    #get the data type information (extracted above)
    tidyresults[[i]] <- lapply(all_inf, data.frame)[[i]]
    #append result to data type information
    tidyresults[[i]]$result <- result
    #transform to tibble
    tidyresults[[i]] <- as_tibble(tidyresults[[i]])
  
  }
  
  #collapse to data frame
  tidyresults <- do.call("rbind", tidyresults)
  #correct the sampling information
  ##the training dataset matters, not the evaluation dataset, which is always downsampled!
  tidyresults$sampling[221:440] <- "fullsample"
  
  #seperate error from stats
  tidyresults <- tidyresults %>% 
    unnest()
  errs <- tidyresults[seq(1, nrow(tidyresults), 2), ]
  stats <- tidyresults[seq(0, nrow(tidyresults), 2), ]
  
  stats <- stats %>%
    mutate(cor_onsen = map(result, "cor_onsen"),
           cor_midsen = map(result, "cor_midsen"),
           cor_endsen = map(result, "cor_endsen"),
           cor_tsen = map(result, "cor_tsen")) %>%
    select(-result) %>%
    unnest()
  
  #function to extract errors
  get_inf <- function(data, Trait){
    
    d <- data[data$Trait == Trait, ]
    out <- d$mean_error
    
  }
  
  errs <- errs %>%
    mutate(err = map(result, get_inf, "SnsCnp")) %>%
    select(-result) %>%
    unnest()
  
  #combine all performance data
  perf_data <- full_join(stats, errs) %>% 
    #to apply filter() to each row
    rowwise() %>% 
    #keep only out of season evaluations
    filter(!grepl(validexp, trainexp)) %>%
    # #keep only in-season evaluations
    # filter(validexp == trainexp) %>%
    #to cancel rowwise
    ungroup() %>%
    arrange(desc(cor_onsen))
  
  write.csv(perf_data, "O:/Projects/KP0011/1/Analysis/multvarmod/leaveout.csv", row.names = FALSE)
  
  unique(perf_data$datatype)
  
  #group-wise mean
  perf_data %>%
    group_by(datatype) %>%
    summarise(mean_cor = mean(cor_onsen))
  
  perf_data %>%
    group_by(validexp) %>%
    summarise(mean_cor = mean(cor_onsen))
  
  perf_data %>%
    group_by(trainexp) %>%
    summarise(mean_cor = mean(cor_onsen))
  
  perf_data %>%
    group_by(sampling) %>%
    summarise(mean_cor = mean(cor_onsen))

### DONE ###   
###########################################################################################
###########################################################################################

###########################################################################################
#  CUBIST
###########################################################################################

#define tuning parameter grid
train_grid <-   train_grid <- expand.grid(committees = c(1, 5, 10, 20, 50, 100),
                                          neighbors = 0)

bdir <- "O:/Projects/KP0011/1/Analysis/multvarmod/Results/cubist/"

eval_cubist_cross <- function(train, valid, data, 
                              sampling_train, 
                              sampling_val, preProc, 
                              train_grid, bdir){???
  
  #@train: is a list of training sets. Each list element consists of a vector with the year(s) used for training
  #@valid: is a list of validation sets corresponding to ONE training set. Each list element consists of a list of vectors. These vectors contain all possible combinations of validation data sets. 
  #@data: is a list of dataframes, each list element contains all the data of a particular type. 
  #@sampling_train: defines the procedure to assure class balance in the training set. Either "upsampling", "downsampling", or "fullsample"
  #@sampling_val: defines the procedure to assure class balance in the validation set. Either "upsampling", "downsampling", or "fullsample"
  #@preProcess: specifies the preprocessing of spectral data: Either "none", "center", "scale", or "center_scale"
  #@maxcomp: specifies the tunelength
  #@bdir: defines the directory where output data is written to
  
  #for all types of spectral input data
  for(k in names(data)){
    
    #iter
    print(paste("initializing training for following data type:", k, sep = " "))
    
    #select the respective data from the list
    full_df <- data[[k]]
    
    #for each of the possible training datasets
    #i.e. for each year and all combinations of years
    
    out_train <- list()
    out_all <- list()
    
    for (i in 1:length(train)) {
      
      #######################################################################
      #Parameter tuning on training data
      #######################################################################
      
      #iter
      print(paste("==> start train on", train[[i]], sep = " "))
      
      #extract data for training year(s)
      dat_train0 <- full_df %>% filter(Exp %in% train[[i]])
      
      #select training data:
      ##sample data if required
      if(sampling_train == "downsample"){
        dat_train <- dat_train0
        dat_train$SnsCnp <- as.factor(dat_train$SnsCnp)
        dat_train <- downSample(x = dat_train[, -ncol(dat_train)],
                                y = dat_train$SnsCnp,
                                yname = "SnsCnp")
        dat_train$SnsCnp <- as.numeric(as.character(dat_train$SnsCnp))
        dat_train <- dat_train[,-ncol(dat_train)]
      } else if(sampling_train == "upsample"){
        dat_train <- dat_train0
        dat_train$SnsCnp <- as.factor(dat_train$SnsCnp)
        dat_train <- upSample(x = dat_train[, -ncol(dat_train)],
                              y = dat_train$SnsCnp,
                              yname = "SnsCnp")
        dat_train$SnsCnp <- as.numeric(as.character(dat_train$SnsCnp))
        dat_train <- dat_train[,-ncol(dat_train)]
      } else if(sampling_train == "fullsample") {
        dat_train <- dat_train0
      } else print("specify sampling procedure for train data")
      
      #count samples used for training
      nsamp <- nrow(dat_train)
      
      #extract IDs of samples used for training
      ID_train <- dat_train %>% 
        select(-contains("rflt"), -SnsCnp, -SnsFl0)
      
      #keep only predictors
      dat <- dat_train %>% 
        select(contains("rflt"), SnsCnp)
      
      ##### Training data ready #####
      
      #Train model
    
      # Define cross validation procedure
      ##  10-fold cross validation is used to identify the model with the 
      ##  highest prediction accuracy on the hold-out samples. 
      ##  The best model is reduced by removing committees until the increase in the RMSEP exceeds one standard error. 
      
      tr_ctrl <- trainControl(method = "repeatedcv", 
                              number = 10,
                              repeats = 1,
                              savePredictions = TRUE, 
                              selectionFunction = "oneSE",
                              verboseIter = TRUE,
                              allowParallel = TRUE)
      
      # Train PLSR model, limiting the number of components to a maximum 30
      # Difference spectra are centered and scaled prior to model fitting. 
      # oscorepls is used to fit the model (differences among algorithms neglible with large training data sets)
      
      #fix preProcessing
      if(preProc == "none"){
        preProcess <- NULL
      } else if (preProc == "center_scale") {
        preProcess <- c("center", "scale")
      } else if (preProc == "scale"){
        preProcess <- "scale"
      } else if (preProc == "center") {
        preProcess <- "center"
      } else print("preProcessing needs to be specified!")
      
      #train
      train_cubist <- train(
        SnsCnp ~ .,
        data = dat,
        preProcess = preProcess,
        method = "cubist",
        tuneGrid = train_grid,
        trControl = tr_ctrl
      )
      
      # #create final directory name for each function run
      # full_bdir <- paste(bdir, "SnsCnp", sampling_train,  preProc, "", sep = "/")
      # 
      # #create directories if necessary
      # if(!file.exists(full_bdir)) {
      #   dir.create(full_bdir, recursive = TRUE)}
      # 
      # saveRDS(train_cubist, file = paste(full_bdir, k, "_train.rds", sep = ""))
      
      #######################################################################
      #Model evaluation
      #######################################################################
      
      #######################################################################
      #extract predictions for observations that were used for training
      #######################################################################
      
      predobs_cv <- plyr::ldply(list(train_cubist),
                                function(x) plyr::match_df(x$pred, x$bestTune),
                                .id = "model")
      
      #Average predictions of the hold-out samples and the corresponding observations of the k folds
      obs <- aggregate(obs ~ rowIndex, 
                       data = predobs_cv, mean)
      
      pred <- aggregate(pred ~ rowIndex,
                        data = predobs_cv, mean)
      
      predobs <- full_join(pred, obs, 
                           by = "rowIndex")
      
      predobs_cv_avg <- semi_join(predobs, predobs_cv, 
                                  by = "rowIndex")
      
      #add observation ID
      ID_train$rowIndex <- 1:nrow(ID_train) #add rowIndex to spectral data
      
      predobs <- dplyr::full_join(ID_train, predobs_cv_avg, 
                                  by = "rowIndex")
      
      #Replace extreme values by NA
      predobs[predobs$pred < -5 | predobs$pred > 15, "pred"] <- NA
      
      min_max <- c(min(predobs$pred, na.rm = TRUE), max(predobs$pred, na.rm = TRUE))
      
      #calculate RMSE
      RMSE <- round(plyr::match_df(train_cubist$result, train_cubist$bestTune)[,"RMSE"], 2)
      
      #######################################################################
      #extract predictions for all observations in training year
      #######################################################################
      
      #extract observation IDs
      ID_train0 <- dat_train0 %>% 
        select(-contains("rflt"), -SnsCnp, -SnsFl0)
      
      #predict new data using trained model
      pred <- caret::predict.train(train_cubist, newdata = dat_train0, 
                                   type = "raw", na.action = na.omit)
      
      predobs_all <- cbind(ID_train0,pred, dat_train0["SnsCnp"])
      names(predobs_all) <- gsub("SnsCnp", "obs", names(predobs_all))
      
      #Replace extreme values by NA
      predobs_all[predobs_all$pred < -5 | predobs_all$pred > 15, "pred"] <- NA
      
      min_max <- c(min(predobs_all$pred, na.rm = TRUE), max(predobs_all$pred, na.rm = TRUE))
      
      #######################################################################
      
      #create a list containing all the output for the training dataset and the internal validation
      out_train[[i]] <- list(train[[i]], train[[i]], 
                             k,
                             train_cubist$bestTune, 
                             RMSE, sampling_train, nsamp, preProc,
                             predobs,
                             predobs_all) 
      
      #rename list elements
      names(out_train[[i]]) <- c("training", "validation", 
                                 "datatype",
                                 "bestTune", 
                                 "RMSE", "sampling","nsamp", "preProc",
                                 "data",
                                 "data_all")
      
      #############Training and validation on the same year DONE#############
      
      #######################################################################
      # Model validation on validation data (different Experiment) 
      #######################################################################
      
      #create empty list to take up predobs data
      out_validate <- list()
      
      #for each validation dataset
      for(j in 1:length(valid[[i]])){
        
        #extractcorresponding data from the complete dataset
        dat_val0 <- full_df %>% 
          filter(Exp %in% valid[[i]][[j]])
        
        #sampling procedure to assure class balance
        #downsampling
        if(sampling_val == "downsample"){
          dat_val <- dat_val0
          dat_val$SnsCnp <- as.factor(dat_val$SnsCnp)
          dat_val <- downSample(x = dat_val, ############## CORRECTED!!?!!!!?
                                y = dat_val$SnsCnp,
                                yname = "SnsCnp")
          dat_val$SnsCnp <- as.numeric(as.character(dat_val$SnsCnp))
          dat_val <- dat_val[, -ncol(dat_val)]
        } else if (sampling_val == "upsample"){
          dat_val <- dat_val0
          dat_val$SnsCnp <- as.factor(dat_val$SnsCnp)
          dat_val <- upSample(x = dat_val[, -ncol(dat_val)],
                              y = dat_val$SnsCnp,
                              yname = "SnsCnp")
          dat_val$SnsCnp <- as.numeric(as.character(dat_val$SnsCnp))
          dat_val <- dat_val[, -ncol(dat_val)]
        } else if(sampling_val == "fullsample") {
          dat_val <- dat_val0
        } else print("specify sampling procedure for validation data")
        
        #count samples used for validation
        nsamp <- nrow(dat_val)
        
        ID_val_sub <- dat_val %>% 
          select(-contains("rflt"), -SnsCnp, -SnsFl0)
        
        #remove IDs from validation data
        dat_val_sub <- dat_val %>% 
          select(contains("rflt"), SnsCnp, -SnsFl0)
        
        #predict new data using trained model
        pred <- caret::predict.train(train_cubist, newdata = dat_val_sub, 
                                     type = "raw", na.action = na.omit)
        
        predobs <- cbind(ID_val_sub, pred, dat_val["SnsCnp"])
        names(predobs) <- gsub("SnsCnp", "obs", names(predobs))
        
        #Replace extreme values by NA
        predobs[predobs$pred < -5 | predobs$pred > 15, "pred"] <- NA
        
        min_max <- c(min(predobs$pred, na.rm = TRUE), max(predobs$pred, na.rm = TRUE))
        
        RMSE <- round(hydroGOF::rmse(predobs$pred, predobs$obs),2)
        
        #######################################################################
        #extract predictions for all observations in validation year
        #######################################################################
        
        #extract observation IDs
        ID_val0 <- dat_val0 %>% 
          select(-contains("rflt"), -SnsCnp, -SnsFl0)
        
        #predict new data using trained model
        pred <- caret::predict.train(train_cubist, newdata = dat_val0, 
                                     type = "raw", na.action = na.omit)
        
        predobs_all <- cbind(ID_val0, pred, dat_val0["SnsCnp"])
        names(predobs_all) <- gsub("SnsCnp", "obs", names(predobs))
        
        min_max <- c(min(predobs_all$pred, na.rm = TRUE), max(predobs_all$pred, na.rm = TRUE))
        
        #Replace extreme values by NA
        predobs_all[predobs_all$pred < -5 | predobs_all$pred > 15, "pred"] <- NA
        
        #create a list containing all the output for the training dataset and the internal validation
        out_validate[[j]] <-  list(train[[i]], valid[[i]][[j]],
                                   k,
                                   train_cubist$bestTune, 
                                   RMSE,
                                   sampling_val, nsamp, preProc,
                                   predobs,
                                   predobs_all)
        
        #rename list elements
        names(out_validate[[j]]) <- c("training", "validation", 
                                      "datatype",
                                      "bestTune", 
                                      "RMSE",
                                      "sampling", "nsamp", "preProc",
                                      "data",
                                      "data_all")
        
      } #end of validation sets loop
      
      #combine the training and validation results
      out_all[[i]] <- c(list(out_train[[i]]), out_validate)
      
    }#end of training set loop
    
    #save output
    
      #create final directory name for each function run
      full_bdir <- paste(bdir, "SnsCnp", sampling_train,  preProc, "", sep = "/")

      #create directories if necessary
      if(!file.exists(full_bdir)) {
        dir.create(full_bdir, recursive = TRUE)}

      saveRDS(out_all, file = paste(full_bdir, k, ".rds", sep = ""))

  } #end of data type loop
  
} #end of function

###########################################################################################

#Example run

library(parallel)
library(doParallel)
cluster <- makeCluster(10, outfile = "") # convention to leave 1 core for OS
registerDoParallel(cluster)
clusterEvalQ(cluster, library(doParallel))
clusterEvalQ(cluster, .libPaths("T:/R3UserLibs"))
clusterEvalQ(cluster, library(doParallel))

start.time <- Sys.time()

eval_cubist_cross(train = train,
                  valid = valid, 
                  data = dat1,
                  sampling_train = "fullsample",
                  sampling_val = "downsample",
                  preProc = "center_scale",
                  train_grid = train_grid,
                  bdir = bdir)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

stopCluster(cluster)
registerDoSEQ()

###########################################################################################

# create predobs plots: for one directory

preProc <- sapply(strsplit(unlist(dirs), "\\/"), "[[", 2)

dir <- "O:/Projects/KP0011/1/Analysis/multvarmod/Results/cubist/SnsCnp/fullsample/center_scale"

preProc <- "center_scale"

files <- as.list(list.files(paste(dir), pattern = ".rds"))

setwd(dir)

data <- lapply(files, readRDS)
names(data) <- sapply(strsplit(unlist(files), "\\."), "[[", 1)

for(i in names(data)){
  
  d <- data[[i]]
  
  D <- unlist(d, recursive = FALSE)
  
  plot_predobs <- function(data) {
    
    data$data$obsf <- as.factor(data$data$obs)
    
    Plot <- ggplot() +
      geom_boxplot(aes(x = data$data$obsf, y = data$data$pred),
                   size = 0.5, alpha = 1) +
      stat_summary(mapping = aes(x = data$data$obs + 1, y = data$data$pred), fun.y = "median", colour = "red", size = 3, geom = "point", shape=18) +
      scale_y_continuous(breaks = seq(0,10,2), limits = c(-2, 12)) +
      geom_abline(intercept = -1, slope = 1, color = "red", size = 1, linetype="dashed") +
      xlab("observed") + ylab("predicted (10-fold CV)") +
      annotate("text", x = 9, y = 0, label= paste("RMSE = ", data$RMSE, sep = ""), size = 2.5) +
      annotate("text", x = 9, y = -0.5, label = paste("Coms = ", data$bestTune$committees, " ", 
                                                      "Nbrs = ", data$bestTune$neighbors, sep = ""), size = 2.5) +
      annotate("text", x = 9, y = -1, label = paste("nsamp = ", data$nsamp, sep = ""), size = 2.5) +
      ggtitle(paste("training = ", paste(data$training, collapse = ", "),
                    "\nvalidation = ", paste(data$validation, collapse = ", "),
                    paste("\ndata = ", i, sep = ""), paste("\npreProc = ", preProc, sep = ""), sep = "")) +
      theme_bw() +
      theme(plot.title = element_text(size = 7, face = "bold"))
    
    return(Plot)
    
  }
  
  Plots <- lapply(D, plot_predobs)
  
  library(gridExtra)
  # arrange them in a 2x2 grid
  
  ml <- marrangeGrob(Plots, nrow=2, ncol=2)
  
  ggsave(paste("O:/FIP/2017/WW018/ASD/Renamed/RData/Results/cubist/fullsample/", 
               preProc, "/", i, ".pdf", sep = ""), 
         ml, 
         width = 297, 
         height = 210, 
         units = "mm")
  
}

###########################################################################################

# further analysis steps
## for all types of data and all pre-treatments,
## extract model information and obtained predictions

# load data
basedir <- "O:/Projects/KP0011/1/Analysis/multvarmod/Results_selvars/cubist"
basedir <- "O:/Projects/KP0011/1/Analysis/multvarmod/Results/cubist/SnsCnp/fullsample/center_scale"
fulldirs <- as.list(list.files(basedir, pattern = ".rds", recursive = TRUE, full.names = TRUE))[1] #make sure the "smth_avg... is
data_all <- lapply(fulldirs, readRDS)

# function to extract model information

ex_inf <- function(data) {
  
  RMSEs <- sapply(data, function(x) sapply(x, '[[', 5)) %>% Reduce(c, .)
  trainexp0 <- lapply(data, function(x) lapply(x, '[[', 1))
  trainexp <-  as.vector(sapply(trainexp0, function(x) sapply(x, paste, collapse = "_"))) %>% Reduce(c, .)
  validexp0 <- lapply(data, function(x) lapply(x, '[[', 2))
  validexp <- as.vector(sapply(validexp0, function(x) sapply(x, paste, collapse = "_"))) %>% Reduce(c, .)
  preProc <- as.vector(sapply(data, function(x) sapply(x, '[[', 8))) %>% Reduce(c, .)
  datatype <- as.vector(sapply(data, function(x) sapply(x, '[[', 3))) %>% Reduce(c, .)
  sampling <- as.vector(sapply(data, function(x) sapply(x, '[[', 6))) %>% Reduce(c, .)
  inf <- cbind(sampling, 
               datatype,
               preProc, trainexp, validexp, RMSEs)
  
  return(inf)
  
} #check whether datatype is available. Else adjust <===== !!!!!!!!!!!!!!!!!!!!!

# extract model information

all_inf <- lapply(data_all, ex_inf)

# extract model predictions

preds <- lapply(data_all, function(y) lapply(y, function(x) lapply(x, '[[', 10)) %>% Reduce(c, .)) #check wheter datatype is available. Else adjust <===== !!!!!!!!!!!!!!!!!!!!!

# select the validation data: only experiment used to train
# and the combination of the experiments not used to train
subset <- lapply(preds, "[", c(1, 4, 5, 8, 9, 12, 13:19))
subset_inf <- lapply(all_inf, "[", c(1, 4, 5, 8, 9, 12, 13:19), )

###########################################################################################

# wrapper to evaluate scaled predictions obtained from all models 
# Generate correlation information, mean error

########################################################################################### -

# Load scoring data
sen_d <- readRDS("T:/PhD/DataAnalysis/FPWW012_FPWW018/Data/sen_d_clean.rds")

# define parameters to be evaluated
ind <- "pred"

# wrapper function

evaluate_plspreds <- function(data){
  
  #this function uses several costum functions!!
  
  #remove obs values from dataframe
  eval_data <- data %>% 
    select(-obs) %>%
    #scale predictions to range [0:10]
    f_scale_pred() %>%
    #select required variables
    select(Plot_ID, spc_ID, meas_date, pred) %>%
    #join to scaled scorings
    f_match_join(sen_d) %>%
    #reorder columns as required by main evaluation function
    dplyr::select(heading_DAS, heading_GDDAS, grading_DAS, grading_DAH, grading_GDDAS, 1:12, SnsFl0, SnsCnp, everything()) %>%
    #evaluate predictions using custom function
    f_calc_err(method = "lin", ind = ind, traits = c("SnsFl0", "SnsCnp"))
  #extract correlation coefficients for dynamics parameters and total error
  means2 <- extract_inf(eval_data, trait = "SnsCnp", out = "means_all")
  stat_data1 <- extract_inf(eval_data, trait = "SnsCnp", out = "stat_data_all")
  
  out <- list(means2, stat_data1)
  
  return(out)
  
}

# use the wrapper to extract all the information for all data types

OUT <- lapply(preds, function(x) lapply(x,  evaluate_plspreds))

saveRDS(OUT, "O:/Projects/KP0011/1/Analysis/multvarmod/Results/cubist/SnsCnp/fullsample/center_scale/corr/corr.rds")


###########################################################################################

OUT <- readRDS("O:/Projects/KP0011/1/Analysis/multvarmod/Results/cubist/SnsCnp/fullsample/center_scale/corr/corr_restr.rds")

# Rearrange output    

tidyresults <- list()

for (i in 1:length(OUT)){
  
  #get result each data type
  result <- OUT[[i]]
   #get the data type information (extracted above)
  tidyresults[[i]] <- lapply(all_inf, data.frame)[[i]]
  #append result to data type information
  tidyresults[[i]]$result <- result
  #transform to tibble
  tidyresults[[i]] <- as_tibble(tidyresults[[i]])
  
}

#collapse to data frame
tidyresults <- do.call("rbind", tidyresults)
#correct the sampling information
##the training dataset matters, not the evaluation dataset, which is always downsampled!
tidyresults$sampling[177:198] <- "fullsample"

#seperate error from stats
tidyresults <- tidyresults %>% 
  unnest()
errs <- tidyresults[seq(1, nrow(tidyresults), 2), ]
stats <- tidyresults[seq(0, nrow(tidyresults), 2), ]

stats <- stats %>%
  mutate(cor_onsen = map(result, "cor_onsen"),
         cor_midsen = map(result, "cor_midsen"),
         cor_endsen = map(result, "cor_endsen"),
         cor_tsen = map(result, "cor_tsen")) %>%
  select(-result) %>%
  unnest()

#function to extract errors
get_inf <- function(data, Trait){
  
  d <- data[data$Trait == Trait, ]
  out <- d$mean_error
  
}

errs <- errs %>%
  mutate(err_Cnp = map(result, get_inf, "SnsCnp")) %>%
  select(-result) %>%
  unnest()

#combine all performance data
perf_data <- full_join(stats, errs) %>% 
  #to apply filter() to each row
  rowwise() %>% 
  #keep only out of season evaluations
  filter(!grepl(validexp, trainexp)) %>%
  # #keep only in-season evaluations
  # filter(validexp == trainexp) %>%
  #to cancel rowwise
  ungroup() %>%
  arrange(desc(cor_onsen))

write.csv(perf_data, "O:/Projects/KP0011/1/Analysis/multvarmod/final_table/across_rfe.csv", row.names = FALSE)

unique(perf_data$datatype)

#group-wise mean
group_by(datatype) %>%
  summarise(mean_cor = mean(cor_onsen))

### DONE ###   


