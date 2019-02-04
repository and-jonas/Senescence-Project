#====================================================================================== -

#HEADER ----

# AUTHOR: Jonas Anderegg

# Prepare spectral Index dataset
# Extract error and senescence dynamics parameters from SI and scorings
# Extract correlations, mean errors and SI bias

#====================================================================================== -

.libPaths("T:/R3UserLibs")

library(dplyr)
library(prospectr)
library(tidyr)
library(stringr)
library(readxl)
library(nls.multstart)
library(purrr)
library(broom)
library(data.table)

source("O:/Projects/KP0011/3/Project/Septoria-Project/R/Utils/001_spectra_Utils.R")
source("O:/Projects/KP0011/1/Senescence-Project/R/Utils/001_A_spectra_utils.R")

dir <- "O:/Projects/KP0011/1/Senescence-Project/Data/"
setwd(dir)

#====================================================================================== -

#load smoothed, averaged, trimmed spectra 
spc <- readRDS("Spectra/smth_avg_rflt_untrimmed.rds")

#exclude 2018-06-10
spc <- spc[!spc$meas_date == "2018-06-10",]

#calculate and scale SI
SI <- f_calc_si(spc) %>%
  f_scale_si(sub = "post_heading")

saveRDS(SI, file = "SI/SI_sc_posthead.rds")

#====================================================================================== -

#load gddah_data and match_dates
gddah <- read.csv("Helper_files/gddah_data.csv") %>% 
  mutate(meas_date = meas_date %>% as.Date()) %>% 
  dplyr::select(Plot_ID, meas_date, meas_GDDAH)
match_dates <- read_excel("Helper_files/match_dates.xlsx")

#add design and measurement information to data 
spc <- readRDS("SI/SI_sc_posthead.rds")
scr <- readRDS("Scorings/scorings_scaled.rds")

data <- f_match_join(spc, scr, gddah_data, match_dates) %>% 
  dplyr::select(-contains("SI_"), everything()) %>% 
  filter(complete.cases(.))

#Define subset of SI for which to carry out the analysis
ind <- readRDS("T:/PhD/Data_Analysis_3/Trait_pred/Data/SI_decorr.rds")
#add narrow-band NDVI to the dataset as benchmark
ind <- c(ind, "SI_NDVI_nb_ASD")

# extract senescence dynamics parameters from spectral indices and scorings;
# from linear interpolations or parametric models;
# calculate area between the interpolation curves;
# calculate bias as difference;
# do for each SI, plot, trait
result <- f_calc_err(data = data, 
                     method = "lin", 
                     ind = ind, 
                     traits = c("SnsFl0", "SnsCnp")) 


# extract information 
means1 <- extract_inf(result, trait = "SnsCnp", out = "means_exp")
# saveRDS(means1, "O:/Projects/KP0011/1/Analysis/Senescence dynamics/Data/means_exp.rds")
means2 <- extract_inf(result, trait = "SnsCnp", out = "means_all")
# saveRDS(means2, "O:/Projects/KP0011/1/Analysis/Senescence dynamics/Data/means_all.rds")
means_Fl0 <- extract_inf(result, trait = "SnsFl0", out = "means_exp")
# saveRDS(means_Fl0, "O:/Projects/KP0011/1/Analysis/Senescence dynamics/Data/means_all_Fl0.rds")
stat_data0 <- extract_inf(result, trait = "SnsCnp", out = "stat_data_exp")
# saveRDS(stat_data0, "O:/Projects/KP0011/1/Analysis/Senescence dynamics/Data/stat_data_exp.rds")
stat_data1 <- extract_inf(result, trait = "SnsCnp", out = "stat_data_all")
# saveRDS(stat_data1, "O:/Projects/KP0011/1/Analysis/Senescence dynamics/Data/stat_data_all.rds")
stat_data_Fl0 <- extract_inf(result, trait = "SnsFl0", out = "stat_data_all")
# saveRDS(stat_data_Fl0, "O:/Projects/KP0011/1/Analysis/Senescence dynamics/Data/stat_data_all_Fl0.rds")



