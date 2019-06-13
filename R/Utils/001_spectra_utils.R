# MAIN FUNCTIONS -----

# read senescence rating data
f_sen_read <- function(dir, file_names) {
  
  #list files
  setwd(dir)
  files <- as.list(file_names)
  
  #read files
  d <- lapply(files, function(i){
    x <- read.csv(i)
    x
  })
  
  #add Experiment ID
  #build final tibble
  data <- do.call("rbind", d) %>% mutate(Exp = str_sub(Plot_ID, 1, 7)) %>% 
    dplyr::select(Exp, everything()) %>% tibble::as_tibble() %>% 
    #transform to dates
    mutate(grading_date = grading_date %>% as.Date(),
           heading_date = heading_date %>% as.Date())
  
  return(data)
  
}

# invert ratings
f_invert_sen <- function(data) {
  
  data <- data %>% mutate(SnsFl0 = SnsFl0 %>% revert(),
                          SnsCnp = SnsCnp %>% revert()) %>% 
    arrange(Plot_ID, grading_date)

  return(data)
  
}

# scale ratings
f_scale_sen <- function(data) {
  
  #drop plots which contain missing values in scoring time series
  data <- data %>% group_by(Plot_ID) %>% filter(!any(is.na(SnsCnp))) %>% ungroup() %>% droplevels()
  #select and group
  df_sc <- data %>% arrange(Plot_ID) %>% dplyr::select(Plot_ID, SnsFl0, SnsCnp) %>% 
    nest(-Plot_ID) %>% 
    #scale scorings on a plot level
    mutate(data = purrr::map(data, col_mapping)) %>% 
    unnest() %>% 
    #complete tibble
    bind_cols(data %>% dplyr::select(-SnsCnp, -SnsFl0, -Plot_ID), .) %>% 
    #reorder columns
    dplyr::select(Plot_ID, everything())

  return(df_sc)
  
}

# read *.asd files
f_spc_read <- function(dir, Exp) {
  
  #files have to be named measDate_"ASD"_PlotID_"Can"_rep.asd
  
  #load data
  setwd(dir)
  fnames <- list.files(dir, pattern = ".asd")
  ASD <- prospectr::readASD(fnames, 'binary', 'list')
  rflt <- sapply(ASD, "[[", "reflectance")
  wvlt <- sapply(ASD, "[[", "wavelength")
  names <- sapply(fnames, basename)
  
  #transpose to wide format
  df <- cbind(wvlt[,1], rflt)
  df <- t(df)[-1,] %>% as.data.frame()

  #rename columns
  colnames(df) <- unique(paste("rflt", wvlt, sep = "_"))
  
  #extract info from spectra name
  pn0 <- names %>% strsplit(., "_", fixed = TRUE) %>% unlist()
  Plot_ID <- grep(Exp, pn0, value = TRUE)
  replicate <- stringr::str_sub(names,-5,-5)
  meas_date <- stringr::str_sub(strptime(grep("201", names, value = TRUE), "%Y%m%d"), 1, 10) %>% as.Date()

  #create data frame
  df <- data.frame(Plot_ID, meas_date, replicate, df)
  
  #drop reference from data frame
  df <- df[!grepl("Ref", df$Plot_ID),]
  df$Plot_ID <- df$Plot_ID %>% droplevels()
  df <- tibble::as_tibble(df)
  
  return(df)
  
}

# perform continuum removal
f_cont_rem <- function(data, interpol = "linear", method = "substraction") {
  
  #Prepare required format
  #convert to df and remove info
  data <- data %>% as.data.frame() #to enable rownames
  dat <- data %>% dplyr::select(contains("rflt_"))
  #strip rflt from names
  names(dat) <- as.character(gsub("_", "", str_sub(names(dat), -4, -1)))
  #add rownames
  rownames(dat) <- paste(1:nrow(dat))
  dat <- as.matrix(dat)
  # Define wavelengths
  wvlt <- as.numeric(colnames(dat))
  
  # continuum removal
  spc_cr <- prospectr::continuumRemoval(X = dat, wvlt, type = "R", 
                                        interpol = interpol, method = method) %>% tibble::as_tibble()
  
  #reassemble tibble
  #correct names
  names(spc_cr) <- paste("rflt_", names(spc_cr), sep = "")
  #add info and convert to tibble
  dat <- data %>% dplyr::select(-contains("rflt_")) %>% 
    dplyr::bind_cols(spc_cr) %>% tibble::as_tibble()
  
  return(dat)

}

# define matching dates and join data sets
f_match_join <- function(spc, sen, gddah, matches) {
  
  #matching dates lookup table
  match_dates <- matches %>% transmute(scor_date = scor %>% as.Date(),
                                       meas_date = meas %>%  as.Date(),
                                       ref_date = reference_date %>% as.Date())
  
  #assign scoring dates their measurement date, 
  #if possible (max 1 day difference)
  sen_data <- sen %>% full_join(., match_dates, by = c("grading_date" = "scor_date"))

  #convert gddah data to tibble
  gddah_data <- gddah %>% tibble::as_tibble() %>% 
    mutate(meas_date = meas_date %>% as.Date(),
           heading_date = heading_date %>% as.Date())
  
  #add gddah data to spectral data
  spc_data <- right_join(gddah_data, spc, by = c("Plot_ID", "meas_date")) %>% 
    #Plots for which heading date is missing must be excluded  
    filter(!is.na(heading_date)) 
  
  #join spectral and scoring data
  data <- full_join(spc_data, sen_data) %>% 
    #reorder columns
    dplyr::select(-contains("rflt_"), everything()) %>% 
    dplyr::select(Exp, Plot_ID, everything()) %>% 
    arrange(Plot_ID, ref_date)

}

# extract errors and senescence dynamics parameters
get_errors_and_dynpars <- function(data, method){
  
  if (method != "lin"){
    
    if (method == "log"){
      
      #fit logistic model to scoring and SVI values
      m1 <- nls(formula = as.formula("scoring ~ A + C/(1+exp(-b*(grading_GDDAH-M)))"), data = data, 
                start = list(A = 10, C = -10, b = 0.01, M = 675), na.action = na.exclude)
      m2 <- nls(formula = as.formula("value ~ A + C/(1+exp(-b*(meas_GDDAH-M)))"), data = data, 
                start = list(A = 10, C = -10, b = 0.01, M = 675), na.action = na.exclude)
      
    } else if (method == "cgom"){
      
      #fit constrained Gompertz model to scoring and SVI values
      m1 <- nls(formula = as.formula("scoring ~10*exp(-exp(-b*(grading_GDDAH-M)))"), data = data, 
                start = list(b = 0.01, M = 675), na.action = na.exclude)
      m2 <- nls(formula = as.formula("value ~ 10*exp(-exp(-b*(meas_GDDAH-M)))"), data = data, 
                start = list(b = 0.01, M = 675), na.action = na.exclude)
      
    } else if (method == "gom"){
      
      #fit flexible Gompertz model to scoring and SVI values
      m1 <- nls(formula = as.formula("scoring ~ A+C*exp(-exp(-b*(grading_GDDAH-M)))"), data = data, 
                start = list(A = 0, C = 11, b = 0.01, M = 675), na.action = na.exclude)
      m2 <- nls(formula = as.formula("value ~ A+C*exp(-exp(-b*(meas_GDDAH-M)))"), data = data, 
                start = list(A = 0, C = 11, b = 0.01, M = 675), na.action = na.exclude)
    }
    
    #get min and max and create sequence of values to predict
    r1 <- range(data$grading_GDDAH, na.rm = TRUE)
    xNew1 <- seq(r1[1],r1[2],length.out = 1000)
    
    #create predictions
    y <- predict(m1, list(grading_GDDAH = xNew1))
    y2 <- predict(m2, list(meas_GDDAH = xNew1))
    
    #create output identical to linear interpolation (see below)
    l <- list(xNew1, y)
    ll <- list(xNew1, y2)
    names(l) <- names(ll) <- c("x", "y")
    
  } else if (method == "lin"){
    
    x <- data$meas_GDDAH
    y <- data$value
    
    #linear interpolation of SVI values
    l <- approx(x = x, y = y, 
                xout = seq(round(min(x, na.rm = TRUE),0), 
                           round(max(x, na.rm = TRUE), 0)))
    
    x <- data$grading_GDDAH
    y <- data$scoring
    
    #linear interpolation of Scorings
    ll <- approx(x = x, y = y, 
                 xout = seq(round(min(x, na.rm = TRUE), 0), 
                            round(max(x, na.rm = TRUE), 0), 1))
  }
  
  #Calculate error
  d1 <- as.data.frame(do.call("cbind", l))
  d2 <- as.data.frame(do.call("cbind", ll))
  d3 <- merge(d1, d2, by = "x")
  d3$y <- abs(d3$y.x - d3$y.y)
  d3$y.x[is.na(d3$y.x)] <- 0
  d3$y.y[is.na(d3$y.y)] <- 0
  f1 <- approxfun(d3$x, d3$y.x-d3$y.y)     # piecewise linear function
  f2 <- function(x) abs(f1(x))
  Error <- integrate(f2, min(d3$x), max(d3$x), subdivisions = 2000)$value
  
  #Extract senescence dynamics parameters
  onsen_SI <- l[[1]][which(l[[2]] < 8)[1]]
  onsen_Trait <- ll[[1]][which(ll[[2]] < 8)[1]]
  midsen_SI <- l[[1]][which(l[[2]] < 5)[1]]
  midsen_Trait <- ll[[1]][which(ll[[2]] < 5)[1]]
  endsen_SI <- l[[1]][which(l[[2]] < 2)[1]]
  endsen_Trait <- ll[[1]][which(ll[[2]] < 2)[1]]
  tsen_SI <- endsen_SI - onsen_SI
  tsen_Trait <- endsen_Trait - onsen_Trait
  
  #create data frame
  func_out <- do.call(rbind, Map(tibble, "onsen_SI" = onsen_SI, "midsen_SI" = midsen_SI,
                                 "endsen_SI" = endsen_SI, "tsen_SI" = tsen_SI,
                                 "onsen_Trait" = onsen_Trait, "midsen_Trait" = midsen_Trait,
                                 "endsen_Trait" = endsen_Trait, "tsen_Trait" = tsen_Trait)) %>% 
    #add bias
    mutate(d_onsen = onsen_SI - onsen_Trait,
           d_midsen = midsen_SI - midsen_Trait,
           d_endsen = endsen_SI - endsen_Trait,
           d_tsen = tsen_SI - tsen_Trait) %>% 
    #add error
    mutate(Error = Error)
  
  return(func_out)
  
}

#calculate correlation and p-values
do_cor_test <- function(data, x, y, use = "pairwise.complete.obs", 
                        method = "pearson", return = "estimate") {
  cor <- cor.test(data %>% pull(x), data %>% pull(y), 
                  use = use, method = method) %>% 
    broom::tidy(.) %>% pull(return)
  return(cor)
}

# Extract Index dynamics parameters from linear interpolations and parametric models
# DEPRECATED
f_ind_dyn <- function(data, methods, ind) {
  
  method <- list()
  
  for (i in methods) {
    
    print(paste("starting with interpolation method:", i))
    
    trait <- par <- list()
    
    list_traits <- ind
    
    design <- read.csv("T:/PhD/DataAnalysis/design.csv", sep = ";")
    
    for(k in list_traits){
      
      print(paste("starting with trait:", k ))
      
      onsen <- midsen <- endsen <- tsen <- sen_par <- list()
      
      for (j in unique(data$Plot_ID)){
        
        tryCatch({
          
          # interpolate
          
          if (i != "linear"){
            
            if (i == "logistic"){
              
              m1 <- nls(formula = as.formula(paste(k, "~ A + C/(1+exp(-b*(meas_GDDAH-M)))")), data = data[data$Plot_ID == j,], 
                        start = list(A = 10, C = -10, b = 0.01, M = 675), na.action = na.exclude)
              
            } else if (i == "constrained gompertz"){
              
              m1 <- nls(formula = as.formula(paste(k, "~10*exp(-exp(-b*(meas_GDDAH-M)))")), data = data[data$Plot_ID == j,], 
                        start = list(b = 0.01, M = 675), na.action = na.exclude)
              
            } else if (i == "gompertz"){
              
              m1 <- nls(formula = as.formula(paste(k, "~A+C*exp(-exp(-b*(meas_GDDAH-M)))")), data = data[data$Plot_ID == j,], 
                        start = list(A = 0, C = 11, b = 0.01, M = 675), na.action = na.exclude)
            }
            
            r1 <- range(data[data$Plot_ID == j,]$grading_GDDAH, na.rm = TRUE)
            xNew1 <- seq(r1[1],r1[2],length.out = 1000)
            
            y <- predict(m1, list(grading_GDDAH = xNew1))
            
            l <- list(xNew1, y)
            
            names(l) <- c("x", "y")
            
          } else if (i == "linear"){
            
            x <- data[data$Plot_ID == j,]$meas_GDDAH
            y <- data[data$Plot_ID == j,] %>% pull(k)
            
            l <- approx(x = x, y = y, 
                        xout = seq(round(min(x, na.rm = TRUE),0), round(max(x, na.rm = TRUE), 0)))
          }
          
          # Extract senescence dynamics parameters
          
          onsen[[j]] <- l[[1]][which(l[[2]] < 8)[1]]
          midsen[[j]] <- l[[1]][which(l[[2]] < 5)[1]]
          endsen[[j]] <- l[[1]][which(l[[2]] < 2)[1]]
          tsen[[j]] <- endsen[[j]] - onsen[[j]]
          
        }, error=function(e){cat(paste("ERROR_", j, ":", sep = ""),conditionMessage(e), "\n")})
        
      }
      
      
      # assemble data frames  
      
      sen_par <- do.call(rbind, Map(data.frame, "onsen" = onsen, "midsen" = midsen,
                                    "endsen" = endsen, "tsen" = tsen))
      
      par[[k]] <- sen_par
      
    } 
    
    method[[i]] <- par
    
  }
  
  fun <- function(x) {
    
    all <- do.call("rbind", x)
    id <- rownames(all)
    dings <- strsplit(id, ".", fixed=TRUE)
    trait <- sapply(dings, `[[`, 1)
    pn <- sapply(dings, `[[`, 2)
    all$Plot_ID <- pn
    all$Trait <- trait
    all <- all %>% dplyr::select("Trait", "Plot_ID", everything())
    
    return(all)
    
  }  
  
  all <- lapply(method, fun)
  
  all <- all %>%
    Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by=c("Plot_ID", "Trait")), .)
  
  suf <- str_sub(methods, 1, 3)
  trait <- c("onsen", "midsen", "endsen", "tsen")
  
  colnames <- paste(rep(trait, 4), rep(suf, each = 4), sep = "_")
  colnames(all)[-c(1:2)] <- colnames
  
  return(all)
  
}

# SI perf analysis
# DEPRECATED
f_calc_err <- function(data, method, ind, traits) {
  
  trait <- par <- list()
  
  design <- read.csv("O:/Projects/KP0011/1/Senescence-Project/Data/Helper_files/design.csv", sep = ",")
  
  for(k in traits){
    
    error <- par_0 <- list()
    
    for (i in ind){
      
      print(paste("start with Index", i))
      
      Error <- onsen_SI <- onsen <- onsen_Trait <- midsen_SI <- midsen_Trait <- 
        endsen_SI <-endsen_Trait <- tsen_SI <- tsen_Trait <- sen_par <- list()
      
      dat <- data %>% dplyr::select(-contains("SI_"), i)
      
      for (j in unique(dat$Plot_ID)){
        
        tryCatch({
          
          if (method != "lin"){
            
            if (method == "log"){
              
              m1 <- nls(formula = as.formula(paste(k, "~ A + C/(1+exp(-b*(grading_GDDAH-M)))")), data = data[data$Plot_ID == j,], 
                        start = list(A = 10, C = -10, b = 0.01, M = 675), na.action = na.exclude)
              m2 <- nls(formula = as.formula(paste(i, "~ A + C/(1+exp(-b*(meas_GDDAH-M)))")), data = dat[dat$Plot_ID == j,], 
                        start = list(A = 10, C = -10, b = 0.01, M = 675), na.action = na.exclude)
              
            } else if (method == "cgom"){
              
              m1 <- nls(formula = as.formula(paste(k, "~10*exp(-exp(-b*(grading_GDDAH-M)))")), data = dat[dat$Plot_ID == j,], 
                        start = list(b = 0.01, M = 675), na.action = na.exclude)
              m2 <- nls(formula = as.formula(paste(i, "~10*exp(-exp(-b*(meas_GDDAH-M)))")), data = dat[dat$Plot_ID == j,], 
                        start = list(b = 0.01, M = 675), na.action = na.exclude)
              
            } else if (method == "gom"){
              
              m1 <- nls(formula = as.formula(paste(k, "~A+C*exp(-exp(-b*(grading_GDDAH-M)))")), data = dat[dat$Plot_ID == j,], 
                        start = list(A = 0, C = 11, b = 0.01, M = 675), na.action = na.exclude)
              m2 <- nls(formula = as.formula(paste(i, "~A+C*exp(-exp(-b*(meas_GDDAH-M)))")), data = dat[dat$Plot_ID == j,], 
                        start = list(A = 0, C = 11, b = 0.01, M = 675), na.action = na.exclude)
            }
            
            r1 <- range(dat[dat$Plot_ID == j,]$grading_GDDAH, na.rm = TRUE)
            xNew1 <- seq(r1[1],r1[2],length.out = 1000)
            
            y <- predict(m1, list(grading_GDDAH = xNew1))
            y2 <- predict(m2, list(meas_GDDAH = xNew1))
            
            l <- list(xNew1, y)
            ll <- list(xNew1, y2)
            
            names(l) <- names(ll) <- c("x", "y")
            
          } else if (method == "lin"){
            
            x <- dat[dat$Plot_ID == j,]$meas_GDDAH
            y <- dat[dat$Plot_ID == j,] %>% pull(i)
            
            l <- approx(x = x, y = y, 
                        xout = seq(round(min(x, na.rm = TRUE),0), round(max(x, na.rm = TRUE), 0)))
            
            x <- dat[dat$Plot_ID == j,]$grading_GDDAH
            y <- dat[dat$Plot_ID == j,] %>% pull(k)
            
            ll <- approx(x = x, y = y, 
                         xout = seq(round(min(x, na.rm = TRUE), 0), round(max(x, na.rm = TRUE), 0), 1))
          }
          
          #Calculate errors
          d1 <- as.data.frame(do.call("cbind", l))
          d2 <- as.data.frame(do.call("cbind", ll))
          d3 <- merge(d1, d2, by = "x")
          d3$y <- abs(d3$y.x - d3$y.y)
          d3$y.x[is.na(d3$y.x)] <- 0
          d3$y.y[is.na(d3$y.y)] <- 0
          f1 <- approxfun(d3$x, d3$y.x-d3$y.y)     # piecewise linear function
          f2 <- function(x) abs(f1(x))
          Error[[j]] <- integrate(f2, min(d3$x), max(d3$x), subdivisions = 2000)$value
          
          #Extract senescence dynamics parameters
          onsen_SI[[j]] <- l[[1]][which(l[[2]] < 8)[1]]
          onsen_Trait[[j]] <- ll[[1]][which(ll[[2]] < 8)[1]]
          midsen_SI[[j]] <- l[[1]][which(l[[2]] < 5)[1]]
          midsen_Trait[[j]] <- ll[[1]][which(ll[[2]] < 5)[1]]
          endsen_SI[[j]] <- l[[1]][which(l[[2]] < 2)[1]]
          endsen_Trait[[j]] <- ll[[1]][which(ll[[2]] < 2)[1]]
          tsen_SI[[j]] <- endsen_SI[[j]] - onsen_SI[[j]]
          tsen_Trait[[j]] <- endsen_Trait[[j]] - onsen_Trait[[j]]

        }, error=function(e){cat}) #("ERROR :",conditionMessage(e), "\n")
        
      }
      
      #Errors
      Error <- as.data.frame(do.call("rbind", Error))
      rownames <- rownames(Error)
      Error <- cbind(rownames, Error)
      colnames(Error)[1] <- "Plot_ID"
      colnames(Error)[2] <- i
      error[[i]] <- Error
      
      #Parameters
      sen_par <- do.call(rbind, Map(data.frame, "onsen_SI" = onsen_SI, "midsen_SI" = midsen_SI,
                                    "endsen_SI" = endsen_SI, "tsen_SI" = tsen_SI,
                                    "onsen_Trait" = onsen_Trait, "midsen_Trait" = midsen_Trait,
                                    "endsen_Trait" = endsen_Trait, "tsen_Trait" = tsen_Trait))
      
      design$Plot_ID <- as.character(design$Plot_ID)
      
      setDT(sen_par, keep.rownames = TRUE)[]
      setnames(sen_par, 1, "Plot_ID")
      
      dat <- inner_join(design, sen_par, by = "Plot_ID") %>% 
        mutate(d_onsen = onsen_SI - onsen_Trait,
               d_midsen = midsen_SI - midsen_Trait,
               d_endsen = endsen_SI - endsen_Trait,
               d_tsen = tsen_SI - tsen_Trait) %>% 
        arrange(Gen_Name)

      par_0[[i]] <- dat
      
    } 
    
    all_errors <- error %>%
      Reduce(function(dtf1, dtf2) left_join(dtf1, dtf2, by = "Plot_ID"), .)
    
    #Reassemble final dataframe
    trait[[k]] <- all_errors
    par[[k]] <- par_0
    
  }
  
  final <- list(trait, par)
  
  return(final)
  
}

# Tidy up function output
# DEPRECATED
extract_inf <- function(data, trait, out) {
  
  #data: list output from the f_calc_err function
  #trait: either SnsFl0 for flag leaf or SnsCnp for canopy 
  #out: object to be returned: either means1, means2, stat_data0 or stat_data_all
  
  ####  READ THE DESIGN FILE
  
  #read the design file
  design <- read.csv("T:/PhD/DataAnalysis/design.csv", sep = ",")
  
  #### EXTRACT THE ERROR DATA
  
  #extract the error data from the list output
  errs <- data[[1]]
  
  #add the trait to each dataframe in the list based on list element name
  err_data <- mapply(cbind, errs, "Trait" = names(errs), SIMPLIFY=F)
  
  #Collapse list to data frame
  err_data <- do.call("rbind", err_data) %>%
    dplyr::select(Plot_ID, Trait, everything()) %>%
    tidyr::gather(Index, "Error", -c(1:2)) %>%
    #join to design
    inner_join(design) %>%
    dplyr::select(Plot_ID, Exp, Lot, RangeLot, RowLot, Gen_Name, everything())
  
  #exclude Lot 6 from FPWW018, which is uncomplete
  err_data <- err_data[err_data$Lot != 6,] # exclude incomplete data
  
  #get experiment-wise means of errors
  means_exp <- err_data %>%
    dplyr::group_by(Exp, Index, Trait) %>%
    dplyr::summarise(mean_error = mean(Error)) %>%
    #sort according to error
    dplyr::arrange(Exp, Trait, mean_error) %>%
    #extract only the data for the requested trait
    dplyr::filter(Trait == trait)
  
  #get overall means of errors
  means_all <- err_data %>%
    dplyr::group_by(Index, Trait) %>%
    dplyr::summarise(mean_error = mean(Error)) %>%
    #sort according to error
    dplyr::arrange(Trait, mean_error) %>%
    #extract only the data for the requested trait
    dplyr::filter(Trait == trait)
  
  #### GET THE CORRELATIONS
  
  #extract the parameter data
  pars <- data[[2]][[trait]]
  Indices <- names(pars)
  
  #add the index to each dataframe in the list based on list element name
  stat_data <- mapply(cbind, pars, "Index" = Indices, SIMPLIFY = FALSE)
  
  #collapse list into single dataframe
  stat_data <- do.call("rbind", stat_data)
  
  #add trait to dataframe
  stat_data$Trait <- trait
  
  #reorder columns
  stat_data <- stat_data %>% dplyr::select(Plot_ID, Index, Trait, everything())
  
  # exclude incomplete data
  stat_data <- stat_data[stat_data$Lot != 6,] 
  
  # Dataset per experiment
  stat_data_exp <- stat_data %>% 
    group_by(Exp, Trait, Index) %>%
    tidyr::nest() 
  
  # Dataset over all experiments
  stat_data_all <- stat_data %>% 
    group_by(Trait, Index) %>%
    tidyr::nest()
  
  # add correlation coefficients and corresponding p-values to the datasets
  
  stat_data_exp$cor_onsen <- unlist(sapply(stat_data_exp$data, function(x) 
    unname(unlist(broom::tidy(cor.test(x$onsen_SI, x$onsen_Trait, method = "pearson"))["estimate"]))))
  stat_data_exp$p.value_onsen <- unlist(sapply(stat_data_exp$data, function(x) 
    unname(unlist(broom::tidy(cor.test(x$onsen_SI, x$onsen_Trait, method = "pearson"))["p.value"]))))
  stat_data_exp$cor_midsen<- unlist(sapply(stat_data_exp$data, function(x) 
    unname(unlist(broom::tidy(cor.test(x$midsen_SI, x$midsen_Trait, method = "pearson"))["estimate"]))))
  stat_data_exp$p.value_midsen <- unlist(sapply(stat_data_exp$data, function(x) 
    unname(unlist(broom::tidy(cor.test(x$midsen_SI, x$midsen_Trait, method = "pearson"))["p.value"]))))
  stat_data_exp$cor_endsen <- unlist(sapply(stat_data_exp$data, function(x) 
    unname(unlist(broom::tidy(cor.test(x$endsen_SI, x$endsen_Trait, method = "pearson"))["estimate"]))))
  stat_data_exp$p.value_endsen <- unlist(sapply(stat_data_exp$data, function(x) 
    unname(unlist(broom::tidy(cor.test(x$endsen_SI, x$endsen_Trait, method = "pearson"))["p.value"]))))
  stat_data_exp$cor_tsen <- unlist(sapply(stat_data_exp$data, function(x) 
    unname(unlist(broom::tidy(cor.test(x$tsen_SI, x$tsen_Trait, method = "pearson"))["estimate"]))))
  stat_data_exp$p.value_tsen <- unlist(sapply(stat_data_exp$data, function(x) 
    unname(unlist(broom::tidy(cor.test(x$tsen_SI, x$tsen_Trait, method = "pearson"))["p.value"]))))
  
  stat_data_exp$mean_d_onsen <- unlist(sapply(stat_data_exp$data, function(x) mean(x$d_onsen, na.rm = TRUE)))
  stat_data_exp$sd_d_onsen <- unlist(sapply(stat_data_exp$data, function(x) sd(x$d_onsen, na.rm = TRUE)))
  stat_data_exp$mean_d_midsen <- unlist(sapply(stat_data_exp$data, function(x) mean(x$d_midsen, na.rm = TRUE)))
  stat_data_exp$sd_d_midsen <- unlist(sapply(stat_data_exp$data, function(x) sd(x$d_midsen, na.rm = TRUE)))
  stat_data_exp$mean_d_endsen <- unlist(sapply(stat_data_exp$data, function(x) mean(x$d_endsen, na.rm = TRUE)))
  stat_data_exp$sd_d_endsen <- unlist(sapply(stat_data_exp$data, function(x) sd(x$d_endsen, na.rm = TRUE)))
  stat_data_exp$mean_d_tsen <- unlist(sapply(stat_data_exp$data, function(x) mean(x$d_tsen, na.rm = TRUE)))
  stat_data_exp$sd_d_tsen <- unlist(sapply(stat_data_exp$data, function(x) sd(x$d_tsen, na.rm = TRUE)))
  
  stat_data_exp <- stat_data_exp[-4]
  
  # add correlation coefficients and corresponding p-values to the datasets
  
  stat_data_all$cor_onsen <- unlist(sapply(stat_data_all$data, function(x) 
    unname(unlist(broom::tidy(cor.test(x$onsen_SI, x$onsen_Trait, method = "pearson"))["estimate"]))))
  stat_data_all$p.value_onsen <- unlist(sapply(stat_data_all$data, function(x) 
    unname(unlist(broom::tidy(cor.test(x$onsen_SI, x$onsen_Trait, method = "pearson"))["p.value"]))))
  stat_data_all$cor_midsen<- unlist(sapply(stat_data_all$data, function(x) 
    unname(unlist(broom::tidy(cor.test(x$midsen_SI, x$midsen_Trait, method = "pearson"))["estimate"]))))
  stat_data_all$p.value_midsen <- unlist(sapply(stat_data_all$data, function(x) 
    unname(unlist(broom::tidy(cor.test(x$midsen_SI, x$midsen_Trait, method = "pearson"))["p.value"]))))
  stat_data_all$cor_endsen <- unlist(sapply(stat_data_all$data, function(x) 
    unname(unlist(broom::tidy(cor.test(x$endsen_SI, x$endsen_Trait, method = "pearson"))["estimate"]))))
  stat_data_all$p.value_endsen <- unlist(sapply(stat_data_all$data, function(x) 
    unname(unlist(broom::tidy(cor.test(x$endsen_SI, x$endsen_Trait, method = "pearson"))["p.value"]))))
  stat_data_all$cor_tsen <- unlist(sapply(stat_data_all$data, function(x) 
    unname(unlist(broom::tidy(cor.test(x$tsen_SI, x$tsen_Trait, method = "pearson"))["estimate"]))))
  stat_data_all$p.value_tsen <- unlist(sapply(stat_data_all$data, function(x) 
    unname(unlist(broom::tidy(cor.test(x$tsen_SI, x$tsen_Trait, method = "pearson"))["p.value"]))))
  
  stat_data_all$mean_d_onsen <- unlist(sapply(stat_data_all$data, function(x) mean(x$d_onsen, na.rm = TRUE)))
  stat_data_all$sd_d_onsen <- unlist(sapply(stat_data_all$data, function(x) sd(x$d_onsen, na.rm = TRUE)))
  stat_data_all$mean_d_midsen <- unlist(sapply(stat_data_all$data, function(x) mean(x$d_midsen, na.rm = TRUE)))
  stat_data_all$sd_d_midsen <- unlist(sapply(stat_data_all$data, function(x) sd(x$d_midsen, na.rm = TRUE)))
  stat_data_all$mean_d_endsen <- unlist(sapply(stat_data_all$data, function(x) mean(x$d_endsen, na.rm = TRUE)))
  stat_data_all$sd_d_endsen <- unlist(sapply(stat_data_all$data, function(x) sd(x$d_endsen, na.rm = TRUE)))
  stat_data_all$mean_d_tsen <- unlist(sapply(stat_data_all$data, function(x) mean(x$d_tsen, na.rm = TRUE)))
  stat_data_all$sd_d_tsen <- unlist(sapply(stat_data_all$data, function(x) sd(x$d_tsen, na.rm = TRUE)))
  
  stat_data_all <- stat_data_all[-3]
  
  final <- list(means_exp, means_all, stat_data_exp, stat_data_all)
  names(final) <- c("means_exp", "means_all", "stat_data_exp", "stat_data_all")
  
  out_function <- final[[out]]
  
  return(out_function)
  
}

#============================================================================================ -

# HELPER FUNCTIONS ----

# Contrained Gompertz equation
Gompertz_constrained <- function(b, M, tempsum) {
  grenness_decay <- 10*exp(-exp(-b*(tempsum-M)))
  return(grenness_decay)
}

# Logistic equation
logistic <- function(A, C, b, M, tempsum) {
  grenness_decay <- A + C/(1+exp(-b*(tempsum-M)))
  return(grenness_decay)
}

# Flexible Gompertz equation
Gompertz_flex <- function(A, C, b, M, tempsum) {
  grenness_decay <- A + C*exp(-exp(-b*(grading_GDDAH-M)))
  return(grenness_decay)
}

# Linear interpolation
lin_approx <- function(data){
  data <- as.data.frame(data)
  out <- approx(data[,"grading_GDDAH"], data[,"Score"], xout = seq(round(min(data[,"grading_GDDAH"], na.rm = TRUE), 0),
                                                                   round(max(data[,"grading_GDDAH"], na.rm = TRUE), 0), 1))
  names(out) <- c("grading_GDDAH", ".fitted")
  out <- list(out)
  return(out)
}

# Extraction of parameters from fits
extract_pars <- function(data){
  onsen <- data[which(data[2] < 8)[1], "grading_GDDAH"]
  midsen <- data[which(data[2] < 5)[1], "grading_GDDAH"]
  endsen <- data[which(data[2] < 2)[1], "grading_GDDAH"]
  tsen <- endsen - onsen
  pars <- cbind(onsen, midsen, endsen, tsen)
  names(pars) <- c("onsen", "midsen", "endsen", "tsen")
  return(pars)
}
