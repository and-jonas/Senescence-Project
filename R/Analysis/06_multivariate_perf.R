#====================================================================================== -

#HEADER ----

# AUTHOR: Jonas Anderegg

# Fit multivariate models to predict senescence scorings from reflectance spectra
# Evaluate across-year applicability of models

#====================================================================================== -

.libPaths("T:/R3UserLibs_rds")
library(tidyverse)
library(caret)
library(doParallel)
library(Cubist)
library(tictoc)

dir <- "O:/Projects/KP0011/1/Senescence-Project/"
setwd(dir)

source("R/Utils/003_multivariate_Utils.R")

dirfrom <- "Data/Spectra/"
dirto <- NULL

#====================================================================================== -

#Prepare data ----

data <- readRDS(paste(dirfrom, "smth_avg_rflt_bin6_scrunscaled.rds", sep = ""))

# Experiments
Exp <- c("FPWW012", "FPWW018", "FPWW022")

# create train list of length 7
train <- list(Exp[1], Exp[2]
              # , Exp[3], Exp[1:2], Exp[2:3], Exp[c(1, 3)], Exp[1:3]
              )

# create validation list of length 7
test <- list(c(list(Exp[2], Exp[3], Exp[2:3])),
             c(list(Exp[1], Exp[3], Exp[c(1, 3)]))
             # , 
             # c(list(Exp[1], Exp[2], Exp[1:2])),
             # c(list(Exp[3])),
             # c(list(Exp[1])),
             # c(list(Exp[2])),
             # c(list(Exp[1], Exp[2], Exp[3]))
             )

tuneGrid <- expand.grid(committees = c(1, 2, 5, 10, 20),
                        neighbors = 0)

#====================================================================================== -

#Run sequentially ----

tic()
result <- eval_multivariate_cross(Trait = "SnsCnp", method = "pls", 
                                  tuneGrid = NULL, maxcomp = 15,
                                  train = train, test = test, data = data,
                                  trainsample = "fullsample", 
                                  testsample = "upsample",
                                  preProc = "center_scale", runParallel = FALSE)
toc()
#====================================================================================== -

#Run in parallel ----

tic()
result <- eval_multivariate_cross(Trait = "SnsCnp", method = "cubist", 
                                  tuneGrid = tuneGrid, maxcomp = NULL,
                                  train = train, test = test, data = data, 
                                  trainsample = "fullsample", 
                                  testsample = "upsample",
                                  preProc = "center_scale", 
                                  runParallel = TRUE)
toc()

#====================================================================================== -

#Run all in parallel ----
##Create a grid with all possible argument combinations
args1 <- expand.grid(Trait = "SnsCnp",
                     method = c("pls", "cubist"), 
                     preProcess = "center_scale",
                     maxcomp = 25, 
                     dir = "T:/PhD/DataAnalysis/FPWW012_FPWW018/Correlation Analyses/Correlation_PLSR/Data/",
                     subset = c("FPWW012", "FPWW018"))
