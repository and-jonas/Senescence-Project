#====================================================================================== -

#HEADER ----

# AUTHOR: Jonas Anderegg

# Preparation agronomic data and senescence dynamics data 
# for calculation of heritability and BLUEs

#====================================================================================== -

.libPaths("T:/R3UserLibs")
library(tidyverse)

source("O:/Projects/KP0011/1/Senescence-Project/R/Utils/001_A_spectra_utils.R")
source("O:/Projects/KP0011/3/Project/Septoria-Project/R/Utils/001_spectra_Utils.R")

dirfrom <- "T:/PhD/DataAnalysis/FPWW012_FPWW018/"
setwd(dirfrom)
#====================================================================================== -

#Load data ----

#> heading ----

head16 <- read.csv("ANOVA_HD/Data/Heading_FPWW012_DAS.csv")
head17 <- read.csv("ANOVA_HD/Data/Heading_FPWW018_DAS.csv")
head18 <- read.csv("ANOVA_HD/Data/Heading_FPWW022_DAS.csv", sep = ";")

head <- rbind(head16, head17, head18) %>% 
  dplyr::select(Plot_ID, heading_date, heading_GDDAS) %>% 
  mutate(heading_date = as.Date(heading_date)) %>% 
  as_tibble()
#====================================================================================== -

#> GY ----

grnyld16 <- read.csv("ANOVA_yield/Data/FPWW012_GrnYld_2016_08_31.csv")
names(grnyld16)[c(1, 3)] <- c("Plot_ID","GY")
grnyld17 <- read.csv("ANOVA_yield/Data/FPWW018_GrnYld_2017_08_18.csv")
names(grnyld17)[c(1, 3)] <- c("Plot_ID","GY")
grnyld18 <- read.csv("ANOVA_yield/Data/FPWW022_GrnYld_2018_09_12.csv")
names(grnyld18)[c(1, 3)] <- c("Plot_ID","GY")

grnyld <- rbind(grnyld16, grnyld17, grnyld18)
grnyld$Plot_ID <- as.character(grnyld$Plot_ID)

grnyld <- grnyld %>% 
  dplyr::select(Plot_ID, GY) %>% 
  as_tibble()
#====================================================================================== -

#> GPC ----

grnprt16 <- read.csv("ANOVA_yield/Data/FPWW012_NIRS_2017_01_30.csv") %>%
  filter(trait == "GPC")
grnprt17 <- read.csv("ANOVA_yield/Data/FPWW018_NIRS_2018_02_28.csv") %>%
  filter(trait == "GPC")

grnprt <- rbind(grnprt16, grnprt17)
grnprt$plot_id <- as.character(grnprt$plot_id)

names(grnprt)[c(1, 3)] <- c("Plot_ID", "GPC")

grnprt <- grnprt %>% 
  dplyr::select("Plot_ID", "GPC") %>% 
  as_tibble()
#====================================================================================== -

#> Scr ---- 

scr <- f_sen_read(dir = "O:/Projects/KP0011/1/Data preparation/Data/out", 
                    file_names = c("Senescence_FPWW012_DAS.csv", 
                                   "Senescence_FPWW018_DAS.csv", 
                                   "Senescence_FPWW022_DAS.csv")) %>% 
  #transform to wide
  #in order to analyze scoring at each date seperately
  #i.e. as a separate trait
  ##first to FULL long format
  dplyr::select(Plot_ID, grading_date, SnsFl0:SnsCnp) %>% 
  gather(., level, scoring, SnsFl0:SnsCnp, factor_key = TRUE) %>% 
  mutate(Trait = paste(level, grading_date, sep = "_")) %>% 
  dplyr::select(-level, -grading_date) %>% 
  ##then to wide
  spread(., Trait, scoring)
#====================================================================================== -

#> SVIs ----

SVI <- readRDS("O:/Projects/KP0011/1/Senescence-Project/Data/Spectra/smth_avg_rflt_untrimmed.rds") %>% 
  filter(meas_date != "2018-06-10") %>% 
  f_calc_si() %>% 
  #transform to wide
  #in order to analyze scoring at each date seperately
  #i.e. as a separate trait
  ##first to FULL long format
  gather(., SI, value, contains("SI_"), factor_key = TRUE) %>% 
  mutate(Trait = paste(SI, meas_date, sep = "_")) %>% 
  dplyr::select(-SI, -meas_date) %>% 
  ##then to wide
  spread(., Trait, value)
#====================================================================================== -

#> SenDynPars_Scr ---- 

dynpars_scr <- readRDS("O:/Projects/KP0011/1/Senescence-Project/Analysis/Scorings/senpars_scr.R") %>% 
  dplyr::select(-convInfo) %>% 
  #transform to wide
  #in order to analyze scoring at each date seperately
  #i.e. as a separate trait
  ##first to FULL long format
  gather(., parameter, value, onsen_gom:tsen_lin, factor_key = TRUE) %>% 
  mutate(level = gsub("Sns", "", Trait)) %>% 
  mutate(Trait = paste(parameter, level, sep = "_")) %>% 
  dplyr::select(-level, -parameter) %>% 
  ##then to wide
  spread(., Trait, value)
#====================================================================================== -

#> SenDynPars_SVI ----

dynpars_SVI <- readRDS("O:/Projects/KP0011/1/Senescence-Project/Analysis/SI/senpars_SI.rds") %>% 
  #transform to wide
  #in order to analyze scoring at each date seperately
  #i.e. as a separate trait
  ##first to FULL long format
  gather(., parameter, value, onsen_lin:tsen_lin, factor_key = TRUE) %>% 
  mutate(Trait = paste(parameter, Trait, sep = "_")) %>% 
  dplyr::select(-parameter) %>% 
  ##then to wide
  spread(., Trait, value) %>% 
  as_tibble()
#====================================================================================== -

# Combine all ----

Data <- list(head, grnyld, grnprt, scr, SVI, dynpars_scr, dynpars_SVI) %>% 
  Reduce(function(dtf1,dtf2) full_join(dtf1, dtf2, by = "Plot_ID"), .)

#====================================================================================== -

# Add design ----

# Read customized design file for the FPWW012 and FPWW018 experiments
# Must contain a Row and Range variable (Experiment, not Lot!)
design <-  read.csv("T:/PhD/DataAnalysis/design.csv", sep = ",")

# Asreml cannot handle holes in the design;
# therefore add dummy plots
dummy_Plots <- read.csv("T:/PhD/DataAnalysis/FPWW018/dummy_plots.csv", sep = ",") %>% 
  dplyr::select(-Onsen, -Midsen, -Endsen, -tsen)

full_design <- rbind(design, dummy_Plots) %>% 
  as_tibble()

Y2016 <- dplyr::filter(full_design, Exp == "FPWW012") %>% 
  mutate(Rep = ifelse(Lot == 3, 1, 2) %>% as.factor(),
         RowBL = ifelse(Lot == 3, max(RowLot) + 5 + RowLot, RowLot),
         RangeBL = ifelse(Lot == 3, max(RangeLot) + 5 + RangeLot, RangeLot))

Y2017 <- dplyr::filter(full_design, Exp == "FPWW018") %>% 
  mutate(Rep = ifelse(Lot == 2, 1, 2) %>% as.factor(),
         RowBL = ifelse(Lot == 2, max(RowLot) + 5 + RowLot, RowLot),
         RangeBL = ifelse(Lot == 2, max(RangeLot) + 5 + RangeLot, RangeLot))

Y2018 <- dplyr::filter(full_design, Exp == "FPWW022") %>% 
  mutate(Rep = ifelse(Lot == 1, 1, 2) %>% as.factor(),
         RowBL = ifelse(Lot == 1, max(RowLot) + 5 + RowLot, RowLot),
         RangeBL = ifelse(Lot == 1, max(RangeLot) + 5 + RangeLot, RangeLot))

full_design <- rbind(Y2016, Y2017, Y2018) %>% 
  mutate(Xf = RangeBL %>% as.factor(),
         Yf = RowBL %>% as.factor())

Data <- full_join(full_design, Data, by = "Plot_ID") %>%
  as_tibble() %>% 
  arrange(Exp, Lot, RangeLot, RowLot)

#====================================================================================== -

# Adjust data type ----

Data[c(1:7)] <- lapply(Data[c(1:7)], as.factor)

saveRDS(Data, "O:/Projects/KP0011/1/Senescence-Project/Data/AllTraits_ANOVA/AllTraits_ANOVA.rds")

#====================================================================================== -
#====================================================================================== -
#====================================================================================== -
