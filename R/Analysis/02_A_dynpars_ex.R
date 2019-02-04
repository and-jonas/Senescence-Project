#====================================================================================== -

#HEADER ----

# AUTHOR: Jonas Anderegg

# Get senescence dynamics parameters from scorings
# using linear interpolation and gompertz models

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

source("O:/Projects/KP0011/1/Senescence-Project/R/Utils/001_A_spectra_utils.R")

#====================================================================================== -

# extract dynamics parameters ----

data <- readRDS("O:/Projects/KP0011/1/Senescence-Project/Data/Scorings/scorings_scaled.rds") %>% 
  tidyr::gather(Trait, Score, SnsFl0:SnsCnp, factor_key = TRUE) %>% filter(complete.cases(.)) %>% 
  data.frame()

# do linear interpolation of scorings
# and fit nls using multstart
# Performs 250 times repeated NLS fitting (Levenberg-Marquardt algorithm)
# with random-search start parameter sets randomly sampled from a uniform
# distribution between upper and lower starting parameter bounds
data_fits <- data %>%
  nest(-c(Plot_ID, Trait)) %>%
  group_by(Plot_ID, Trait) %>%
  mutate(fit_lin = purrr::map(data, lin_approx)) %>% 
  mutate(fit_cgom = purrr::map(data, ~ nls_multstart(Score ~ Gompertz_constrained(b, M, tempsum = grading_GDDAH),
                                                     data = .x,
                                                     iter = 250,
                                                     start_lower = c(b = -0.1, M = 550),
                                                     start_upper = c(b = 0, M = 750),
                                                     convergence_count = 100)))

# new data frame of predictions
new_preds <- data %>%
  do(., data.frame(grading_GDDAH = seq(min(.$grading_GDDAH), max(.$grading_GDDAH), length.out = 1000), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(data, Plot_ID) %>%
  summarise(., min_gGDDAH = min(grading_GDDAH), max_gGDDAH = max(grading_GDDAH)) %>%
  ungroup()

# create new predictions
preds2 <- data_fits %>%
  tidyr::unnest(fit_cgom %>% purrr::map(broom::augment, newdata = new_preds)) %>%
  merge(., max_min, by = 'Plot_ID') %>%
  group_by(., Plot_ID) %>%
  filter(., grading_GDDAH > unique(min_gGDDAH) & grading_GDDAH < unique(max_gGDDAH)) %>%
  arrange(., Plot_ID, Trait, grading_GDDAH) %>%
  rename(., Score = .fitted) %>%
  ungroup()

# check whether model converged
convInfo <- data_fits %>%
  transmute(convInfo = purrr::map(purrr::map(fit_cgom, "convInfo"), "isConv"))

# # plot
# ggplot() +
#   geom_point(aes(grading_GDDAH, Score, color = Trait), size = 2, data) +
#   geom_line(aes(grading_GDDAH, Score, group = Trait, color = Trait), alpha = 0.5, preds2) +
#   facet_wrap(~ Plot_ID, labeller = labeller(.multi_line = FALSE)) +
#   scale_colour_manual(values = c("green4", "black")) +
#   scale_y_continuous(breaks = seq(0,10,2), limits = c(0, 10)) +
#   theme_bw(base_size = 12, base_family = "Helvetica") +
#   ylab("Visual senescence Score") +
#   xlab("°C days after heading")

#extract parameters from nls fits
preds_gom <- preds2 %>% 
  group_by(Plot_ID, Trait) %>%
  nest() %>%
  mutate(pars_gom = purrr::map(data, extract_pars)) %>%
  select(-data)

#extract parameters from linear interpolation
preds_lin <- data_fits %>%
  unnest(fit_lin) %>% 
  mutate(data_lin = purrr::map(fit_lin, bind_cols)) %>%
  transmute(pars_lin = purrr::map(data_lin, extract_pars))

#combine predictions and convInfo
pred <- list(preds_gom, preds_lin, convInfo) %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by=c("Plot_ID", "Trait")), .)

#combine predictions and convInfo
df <- pred %>% unnest(pars_gom, pars_lin, convInfo)
names(df)[4:length(df)] <- c("onsen_gom", "midsen_gom", "endsen_gom", "tsen_gom",
                             "onsen_lin", "midsen_lin", "endsen_lin", "tsen_lin")

#replace parameter value with NA if model did not converge
df[df$convInfo == FALSE, c(4:7)] <- NA

saveRDS(df, file = "O:/Projects/KP0011/1/Senescence-Project/Analysis/Scorings/senpars_scr.R")

#====================================================================================== -



#Define subset of SI for which to carry out the analysis
ind <- readRDS("T:/PhD/Data_Analysis_3/Trait_pred/Data/SI_decorr.rds")
#add narrow-band NDVI to the dataset as benchmark
ind <- c(ind, "SI_NDVI_nb_ASD")

#calculate Error and pearson correlation coefficients for dynamics parameters
result <- f_calc_err(data = data, method = "lin", ind = ind, traits = c("SnsFl0", "SnsCnp")) 
# saveRDS(result, "O:/Projects/KP0011/1/Analysis/Senescence dynamics/Data/Function_output/Result.rds")

result <- readRDS("O:/Projects/KP0011/1/Analysis/Senescence dynamics/Data/Function_output/Result.rds")  
