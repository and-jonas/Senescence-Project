#====================================================================================== -

#HEADER ----

# AUTHOR: Jonas Anderegg

# Prepare spectral datasets for subsequent analyses

#====================================================================================== -

.libPaths("T:/R3UserLibs")

library(dplyr)
library(prospectr)
library(tidyr)
library(stringr)
library(readxl)

source("O:/Projects/KP0011/3/Project/Septoria-Project/R/Utils/001_spectra_Utils.R")
source("O:/Projects/KP0011/1/Senescence-Project/R/Utils/001_A_spectra_utils.R")

#====================================================================================== -

#Read spectral data ----

#read 2016 and 2017 data
spc1 <- f_spc_read(dir = "O:/FIP/2017/WW018/ASD/Renamed/all", 
                       Exp = "FP")

#read 2018 data (pre-processed in python due to disfunctionalities of the prospectr::readASD() function)
spc2 <- f_spc_read_fromPy(dir = "O:/FIP/2018/WW022/ASD/out_Py", 
                              Exp = "FP")
#combine all data
spc <- rbind(spc1, spc2)

saveRDS(spc, file = "O:/Projects/KP0011/1/Senescence-Project/Data/Spectra/spectra_raw.rds")

#====================================================================================== -

#Read Senescence data ----

scr <- f_sen_read(dir = "O:/Projects/KP0011/1/Data preparation/Data/out",
                    file_names = c("Senescence_FPWW012_DAS.csv", 
                                   "Senescence_FPWW018_DAS.csv",
                                   "Senescence_FPWW022_DAS.csv")) %>%
  f_invert_sen() %>% 
  f_scale_sen()

#Exclude lodging plots
lodg <- read_excel("O:/Projects/KP0011/1/Data preparation/Data/raw/Additional_PhenInfo_ww012.xls")  
lodg_plots <- pull(lodg[!is.na(lodg$Lodging), 1])
scr <- scr %>% filter(!Plot_ID %in% lodg_plots) %>% droplevels()

saveRDS(scr, file = "O:/Projects/KP0011/1/Senescence-Project/Data/Scorings/scorings_scaled.rds")

#====================================================================================== -

#Prep for SI and MultVarMod ----
#> load data ----

dir <- "O:/Projects/KP0011/1/Senescence-Project/Data/"
setwd(dir)

#load raw spectral data (50650 x 2154, 1790 Plots)
spc <- readRDS("Spectra/spectra_raw.rds")

#load unscaled scoring data (2104 Plots)
scr <- readRDS("Scorings/scorings_scaled.rds")

#load gddah_data and match_dates
gddah_data <- read.csv("Helper_files/gddah_data.csv")
match_dates <- read_excel("Helper_files/match_dates.xlsx")

#====================================================================================== -

#> pre-process spectra ----

# 0) mean of smoothed reflectance spectra, untrimmed for SI (14532 x 2143, 1790 Plots)
smth_avg_rflt <- spc %>%
  f_spc_smooth(3, 11, 0) %>%
  f_spc_avg()
saveRDS(smth_avg_rflt, file = "Spectra/smth_avg_rflt_untrimmed.rds")

# 1) mean reflectance spectra
avg_rflt <- spc %>% 
  f_spc_avg() %>%
  f_match_join(., scr, gddah_data, match_dates) %>%
  filter(complete.cases(.))
saveRDS(avg_rflt, file = "Spectra/avg_rflt.rds")

# 2) mean reflectance spectra trimmed
avg_rflt_trim <- spc %>%
  f_spc_avg() %>%
  f_spc_trim() %>%
  f_match_join(., scr, gddah_data, match_dates) %>%
  filter(complete.cases(.))
saveRDS(avg_rflt_trim, file = "Spectra/avg_rflt_trim.rds")

# 2b) mean reflectance spectra restricted to 500nm-700nm
avg_rflt_trim <- spc %>%
  f_spc_avg() %>%
  f_spc_trim() %>%
  f_match_join(., scr, gddah_data, match_dates) %>%
  filter(complete.cases(.)) %>%
  dplyr::select(-contains("rflt_"), rflt_500:rflt_700)
saveRDS(avg_rflt_trim, file = "Spectra/avg_rflt_restr.rds")

# 3) mean of smoothed reflectance spectra
smth_avg_rflt <- spc %>%
  f_spc_smooth(3, 11, 0) %>%
  f_spc_avg() %>%
  f_spc_trim() %>%
  f_match_join(., scr, gddah_data, match_dates) %>%
  filter(complete.cases(.))
saveRDS(smth_avg_rflt, file = "Spectra/smth_avg_rflt.rds")

# 3b) mean of smoothed reflectance spectra, bin
smth_avg_rflt <- spc %>%
  f_spc_smooth(3, 11, 0) %>%
  f_spc_avg() %>%
  f_spc_bin(bin_size = 3) %>%
  f_spc_trim() %>%
  f_match_join(., scr, gddah_data, match_dates) %>%
  filter(complete.cases(.))
saveRDS(smth_avg_rflt, file = "Spectra/smth_avg_rflt_bin3.rds")     

# 4) mean of smoothed 1st derivative spectra
smth_der1_avg <- spc %>%
  f_spc_smooth(3, 11, 1) %>%
  f_spc_avg() %>%
  f_spc_trim() %>%
  f_match_join(., scr, gddah_data, match_dates) %>%
  filter(complete.cases(.))
saveRDS(smth_der1_avg, file = "Spectra/smth_der1_avg.rds")

# 5) remove continuum from smoothed and trimmed reflectance spectra
avg_rflt_cr_smth_der1 <- spc %>%
  f_spc_avg() %>%
  f_spc_trim() %>%
  f_cont_rem() %>%
  f_spc_smooth(3, 11, 0) %>%
  f_match_join(., scr, gddah_data, match_dates) %>%
  filter(complete.cases(.))
saveRDS(avg_rflt_cr_smth_der1, file = "Spectra/avg_rflt_cr_smth.rds")

# 6) remove continuum from smoothed and trimmed reflectance spectra
avg_rflt_cr_smth_der1 <- spc %>%
  f_spc_avg() %>%
  f_spc_trim() %>%
  f_cont_rem() %>%
  f_spc_smooth(3, 11, 1) %>%
  f_match_join(., scr, gddah_data, match_dates) %>%
  filter(complete.cases(.))
saveRDS(avg_rflt_cr_smth_der1, file = "Spectra/avg_rflt_cr_smth_der1.rds")

#====================================================================================== -
