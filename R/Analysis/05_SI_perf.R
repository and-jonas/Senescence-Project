#====================================================================================== -

#HEADER ----

# AUTHOR: Jonas Anderegg

# Extract error and senescence dynamics parameters from SI and scorings
# Extract correlations, mean errors and SI bias

#====================================================================================== -

.libPaths("T:/R3UserLibs")
.libPaths("T:/R3UserLibs_rds")

library(tidyverse)
library(nls.multstart)
library(furrr)

source("O:/Projects/KP0011/1/Senescence-Project/R/Utils/001_Spectra_utils.R")

dir <- "O:/Projects/KP0011/1/Senescence-Project/"
setwd(dir)

dirfrom <- "Data/Matched/"

#====================================================================================== -

data <- readRDS(paste(dirfrom, "SI_scr_compl.rds", sep = ""))

#remove incomplete time series
data <- data[!data$Lot == 6, ] %>% droplevels()
data <- data[!(data$Lot == 3 & data$Exp == "FPWW022"), ]

#Define subset of SI for which to carry out the analysis
ind <- readRDS("T:/PhD/Data_Analysis_3/Trait_pred/Data/SI_decorr.rds")
#add narrow-band NDVI to the dataset as benchmark
ind <- c(ind, "SI_NDVI_nb_ASD")
  
#====================================================================================== -

# SI performance ----
# > Error, bias ----

# extract senescence dynamics parameters from spectral indices and scorings;
# from linear interpolations or parametric models;
# calculate area between the interpolation curves;
# calculate bias as difference;
perf <- data %>% gather(SVI, value, contains("SI_")) %>% 
  gather(Trait, scoring, contains("Sns")) %>%
  filter(SVI %in% ind) %>%
  group_by(Exp, Trait, SVI, Plot_ID) %>% 
  nest() %>% 
  mutate(eval = furrr::future_map(.x = data, method = "lin",
                                  .f = possibly(get_errors_and_dynpars, otherwise = NA_real_)))
  # mutate(eval = purrr::map(.x = data, method = "lin",
  #                          .f = possibly(get_errors_and_dynpars, otherwise = NA_real_)))

#get index performance metrics per experiment and trait
metrics_exp <- perf %>% dplyr::select(-data) %>% 
  unnest() %>% dplyr::select(1:4, contains("d_"), Error) %>%  
  gather(., metric, value, d_onsen:Error) %>% 
  group_by(Exp, Trait, SVI, metric) %>% 
  summarize(mean = mean(value),
            sd = sd(value))

#get errors per experiment
err_exp <- metrics_exp %>% filter(metric == "Error") %>% 
  arrange(Exp, Trait, mean)

#get index performance metrics per trait (across Experiments)
metrics_overall <- perf %>% dplyr::select(-data) %>% 
  unnest() %>% dplyr::select(1:4, contains("d_"), Error) %>%  
  gather(., metric, value, d_onsen:Error) %>% 
  group_by(Trait, SVI, metric) %>% 
  summarize(mean = mean(value),
            sd = sd(value))

#get errors overall
err_overall <- metrics_overall %>% filter(metric == "Error") %>% 
  arrange(Trait, mean)

#====================================================================================== -

# > Correlations ----

#get correlations per experiment and trait
corr_exp <- perf %>% dplyr::select(-data) %>% 
  unnest() %>% dplyr::select(-contains("d_"), -Error) %>%  
  gather(., Dynpar, value, onsen_SI:tsen_Trait) %>% 
  mutate(level = lapply(strsplit(Dynpar, "_"), "[[", 2) %>% unlist(),
         Dynpar = lapply(strsplit(Dynpar, "_"), "[[", 1) %>% unlist()) %>% 
  mutate(level = ifelse(level == "Trait", "scoring", level)) %>% 
  spread(., level, value) %>% 
  arrange(Dynpar) %>% 
  group_by(Exp, Trait, SVI, Dynpar) %>% 
  nest() %>% 
  mutate(pearson_r = purrr::map_dbl(.x = data, x = "scoring", y = "SI", 
                                    .f = do_cor_test),
         p_val = purrr::map_dbl(.x = data, x = "scoring", y = "SI", return = "p.value",
                                .f = do_cor_test)) %>% 
  dplyr::select(-data) %>% 
  arrange(Trait, Exp, Dynpar, desc(abs(pearson_r)))

#get correlations overall
corr_overall <- perf %>% dplyr::select(-data) %>% 
  unnest() %>% dplyr::select(-contains("d_"), -Error) %>%  
  gather(., Dynpar, value, onsen_SI:tsen_Trait) %>% 
  mutate(level = lapply(strsplit(Dynpar, "_"), "[[", 2) %>% unlist(),
         Dynpar = lapply(strsplit(Dynpar, "_"), "[[", 1) %>% unlist()) %>% 
  mutate(level = ifelse(level == "Trait", "scoring", level)) %>% 
  spread(., level, value) %>% 
  arrange(Dynpar) %>% 
  group_by(Trait, SVI, Dynpar) %>% 
  nest() %>% 
  mutate(pearson_r = purrr::map_dbl(.x = data, x = "scoring", y = "SI", 
                                    .f = do_cor_test),
         p_val = purrr::map_dbl(.x = data, x = "scoring", y = "SI", return = "p.value",
                                .f = do_cor_test)) %>% 
  dplyr::select(-data) %>% 
  arrange(Trait, Dynpar, desc(abs(pearson_r)))

#====================================================================================== -
