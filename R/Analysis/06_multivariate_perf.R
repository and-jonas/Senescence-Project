#====================================================================================== -

#HEADER ----

# AUTHOR: Jonas Anderegg

# Fit multivariate models to predict senescence scorings from reflectance spectra
# Evaluate across-year applicability of models

#====================================================================================== -

.libPaths("T:/R3UserLibs_rds")
library(tidyverse)
library(purrrlyr)
library(caret)
library(doParallel)
library(Cubist)
library(tictoc)

rm(list = ls())

dir <- "O:/Projects/KP0011/1/Senescence-Project/"
setwd(dir)

source("R/Utils/003_multivariate_Utils.R")

dirfrom <- "Data/Spectra/MM/"
dirto <- NULL

#====================================================================================== -

#Prepare func args ----
data_rflt <- readRDS(paste(dirfrom, "smth_avg_rflt_bin6_scrunscaled.rds", sep = ""))
data_der1 <- readRDS(paste(dirfrom, "avg_bin_smth_der1.rds", sep = ""))
data_crem <- readRDS(paste(dirfrom, "avg_rflt_cr_smth.rds", sep = ""))
data <- list("rflt" = data_rflt, "der1" = data_der1, "crem" <- data_crem)

# Experiments
Exp <- c("FPWW012", "FPWW018", "FPWW022")

# create train list of length 7
train <- list(Exp[1]
              , Exp[2]
              # , Exp[3], Exp[1:2], Exp[2:3], Exp[c(1, 3)], Exp[1:3]
              )

# create validation list of length 7
test <- list(c(list(Exp[2], Exp[3], Exp[2:3]))
             , c(list(Exp[1], Exp[3], Exp[c(1, 3)]))
             # , 
             # c(list(Exp[1], Exp[2], Exp[1:2])),
             # c(list(Exp[3])),
             # c(list(Exp[1])),
             # c(list(Exp[2])),
             # c(list(Exp[1], Exp[2], Exp[3]))
             )

Grid <- expand.grid(committees = c(1, 2, 5, 10, 20),
                    neighbors = c(0, 1, 2, 4))

#====================================================================================== -

result <- eval_multivariate_cross(Trait = "SnsCnp", method = "pls", 
                                  data_type = "rflt",
                                  testsample = "upsample")
result <- eval_multivariate_cross(Trait = "SnsCnp", method = "cubist", 
                                  data_type = "rflt",
                                  testsample = "upsample")

#====================================================================================== -

#Run all in parallel ----
##Create a grid with all possible argument combinations
args <- expand.grid(Trait = "SnsCnp", method = "pls", 
                    data_type = c("rflt", "der1", "crem"),
                    testsample = "upsample",
                    stringsAsFactors = FALSE)

# Run the model on the different types of data
tic()
output <- args %>% purrrlyr::invoke_rows(.f = eval_multivariate_cross)
toc()

#====================================================================================== -
