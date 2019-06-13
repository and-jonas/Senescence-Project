#====================================================================================== -

#HEADER ----

# AUTHOR: Jonas Anderegg

# Calculate within-year repeatbility (w2)
# Calculate across-year heritability (h2)

#====================================================================================== -

.libPaths("T:/R3UserLibs")
# .libPaths("T:/R3UserLibs_rds")
# update.packages(lib.loc = "T:/R3UserLibs_rds")

library(tidyverse)
library(SpATS)
library(asreml)

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
Exp_wise <- Data %>% 
  dplyr::select(-heading_date) %>% 
  gather(Trait, value, 14:length(.)) %>%
  filter(Trait %in% c("heading_GDDAS", "GY", "GPC") | grepl('sen', Trait)) %>% 
  group_by(Exp, Trait) %>% 
  nest() %>% 
  #calculate within-year repeatablity
  mutate(w2 = purrr::map(.x = data, 
                         .f = possibly(f_spats, otherwise = NA_real_)) %>% 
           purrr::map(.x = ., 
                      .f = possibly(getHeritability, otherwise = NA_real_))) %>% 
  #extract BLUEs
  mutate(BLUE = purrr::map(.x = data, genotype.as.random = FALSE, 
                           .f = possibly(f_spats, otherwise = NA_real_)) %>% 
           purrr::map(.x = ., 
                      .f = possibly(extract_BLUE, otherwise = NA_real_)))

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

Trait_wise <- BLUEs %>% 
  group_by(Trait) %>% nest() %>% 
  #genotype as random to estimate variance component and calculate h2
  mutate(h2 = purrr::map(.x = data, fixed = "Exp", random = "Gen_Name",
                         .f = possibly(get_h2_years, otherwise = NA_real_)),
         #genotype as fixed to estimate the genotypic effect and get BLUE
         BLUE = purrr::map(.x = data, fixed = "Exp + Gen_Name", random = "NULL",
                           .f = possibly(get_BLUE, otherwise = NA_real_)))

saveRDS(Trait_wise, paste(dirto, "Traitwise.rds", sep = ""))
h2_table <- Trait_wise %>% dplyr::select(Trait, h2) %>% unnest()
BLUE <- Trait_wise %>% dplyr::select(Trait, BLUE) %>% filter(!is.na(BLUE)) %>% unnest()

#====================================================================================== -


