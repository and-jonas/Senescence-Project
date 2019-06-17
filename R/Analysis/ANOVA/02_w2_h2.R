#====================================================================================== -

#HEADER ----

# AUTHOR: Jonas Anderegg

# Calculate within-year repeatbility (w2)
# Calculate across-year heritability (h2)

#====================================================================================== -

# .libPaths("T:/R3UserLibs")
.libPaths("T:/R3UserLibs_rds")
# update.packages(lib.loc = "T:/R3UserLibs_rds")

library(tidyverse)
library(SpATS)
library(asreml)
library(furrr)
library(tictoc)

#set working directory
wd <- "O:/Projects/KP0011/1/Senescence-Project/"
setwd(wd)

#load required function
source("R/Utils/002_h2_BLUEs_utils.R")

#define source and sink directories
dirfrom <- "Data/AllTraits_ANOVA/"
dirto <- "Analysis/ANOVA/"

#load data
Data <- readRDS(paste(dirfrom, "AllTraits_ANOVA.rds", sep = ""))

#====================================================================================== -

#w2, BLUEs ----

#calculate within-year repeatability and extract BLUEs for all senesence dynamics traits,
#for GY, GPC and heading date
plan("multiprocess")
Exp_wise <- Data %>% 
  dplyr::select(-heading_date) %>% 
  gather(Trait, value, 14:length(.)) %>%
  filter(Trait %in% c("heading_GDDAS", "GY", "GPC") | grepl('sen', Trait)) %>%
  # filter(Trait %in% c("GY", "GPC", "onsen_gom_Cnp")) %>% 
  group_by(Exp, Trait) %>% 
  nest() %>% 
  #calculate within-year repeatablity
  mutate(w2 = furrr::future_map(.x = data,
                                .f = possibly(f_spats, otherwise = NA_real_)) %>%
           furrr::future_map_dbl(.x = .,
                                 .f = possibly(getHeritability, otherwise = NA_real_))) %>%
  #extract BLUEs
  mutate(BLUE = furrr::future_map(.x = data, genotype.as.random = FALSE,
                                  .f = possibly(f_spats, otherwise = NA_real_)) %>%
           furrr::future_map(.x = .,
                             .f = possibly(extract_BLUE, otherwise = NA_real_))) %>% 
  #get spatial component
  mutate(sp = furrr::future_map(.x = data,
                                .f = possibly(f_spats, otherwise = NA_real_)) %>% 
           furrr::future_map(.x = .,
                             .f = possibly(get_spatial, otherwise = NA_real_))) %>% 
  #plot the spatial trend 
  mutate(plot = furrr::future_map(.x = sp,  
                                  form = formula(spatial ~ RangeLot + RowLot | Lot),
                                  .f = possibly(desplot::desplot, otherwise = NA_real_)))
plan("sequential")

saveRDS("Analysis/ANOVA/Expwise_furrr.rds")

# saveRDS(Exp_wise, paste(dirto, "Expwise.rds", sep = ""))
w2_table <- Exp_wise %>% dplyr::select(Exp, Trait, w2) %>% unnest()

#====================================================================================== -

#h2 ----

Exp_wise <- readRDS("Analysis/ANOVA/Expwise.rds")

#calculate across-year heritability for all senescence dynamics traits,
#for GY, GPC, and heading date
BLUEs <- Exp_wise %>% 
  dplyr::select(Exp, Trait, BLUE) %>% 
  filter(!is.na(BLUE)) %>% unnest()

saveRDS(BLUEs, paste(dirto, "BLUEs.rds", sep = ""))

#get h2 and BLUE
#Requires asreml                <== !!!!!!!!!!
Trait_wise <- BLUEs %>% 
  group_by(Trait) %>% nest() %>% 
  #genotype as random to estimate variance component and calculate h2
  mutate(h2 = purrr::map_dbl(.x = data, fixed = "Exp", random = "Gen_Name",
                             .f = possibly(get_h2_years, otherwise = NA_real_)),
         #genotype as fixed to estimate the genotypic effect and get BLUE
         BLUE = purrr::map_dbl(.x = data, fixed = "Exp + Gen_Name", random = "NULL",
                               .f = possibly(get_BLUE, otherwise = NA_real_)))

saveRDS(Trait_wise, paste(dirto, "Traitwise.rds", sep = ""))

Trait_wise <- readRDS()

h2_table <- Trait_wise %>% dplyr::select(Trait, h2) %>% unnest()
BLUE <- Trait_wise %>% dplyr::select(Trait, BLUE) %>% filter(!is.na(BLUE)) %>% unnest()

#====================================================================================== -


