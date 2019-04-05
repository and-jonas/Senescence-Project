#====================================================================================== -

#HEADER ----

# AUTHOR: Jonas Anderegg

# Summarise results of feature selection

#====================================================================================== -


.libPaths("T:/R3UserLibs")

library(caret)
library(dplyr)
library(data.table)

##get all output data
dir <- "O:/Projects/KP0011/1/Analysis/rfe_cubist/Resamples"
setwd(dir)

paths <- as.list(list.files(dir, recursive = TRUE, full.names = TRUE, pattern = "all"))
data <- lapply(paths, readRDS)

out <- readRDS("O:/Projects/KP0011/1/Analysis/rfe_cubist/Resamples/FPWW012/resample_all.rds")  

colnames <- c("subset_size", paste("Resample", 1:30, sep = ""))

#Function to tidy up the results of recursive feature elimination iterations
#collapses results to two dataframes per year
tidyres <- function(d) {
  
  allranks <- lapply(d, "[[", 1)
  allRMSE <- lapply(d, "[[", 2)
  allRMSETrain <- lapply(allRMSE, "[", c(1:2))
  allRMSETest <- lapply(allRMSE, "[", c(1, 3))
  
  colnames <- c("subset_size", paste("Resample", 1:30, sep = ""))
  
  #create ranks table
  ranks <- allranks %>%
    Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="var"), .)
  names(ranks) <- colnames

  #create test RMSE table
  TestRMSE <- allRMSETest %>%
    Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="subset_size"), .)
  names(TestRMSE) <- colnames

  #create train RMSE table
  TrainRMSE <- allRMSETrain %>%
    Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="subset_size"), .)
  names(TrainRMSE) <- colnames
  
  #average over resamples, get sd of means
  Testperf <- TestRMSE %>% 
    tibble::as.tibble() %>%
    tidyr::gather(resample, RMSE, 2:31) %>%
    group_by(subset_size) %>%
    arrange(subset_size) %>%
    summarise_at(vars(RMSE), funs(mean, sd), na.rm = TRUE) %>%
    mutate(set = "Test") %>% 
    as.data.frame()
  Trainperf <- TrainRMSE %>% 
    tibble::as.tibble() %>%
    tidyr::gather(resample, RMSE, 2:31) %>%
    group_by(subset_size) %>%
    arrange(subset_size) %>%
    summarise_at(vars(RMSE), funs(mean, sd), na.rm = TRUE) %>%
    mutate(set = "Train") %>% 
    as.data.frame()
  
  perf <- rbind(Testperf, Trainperf)
  
  #mean ranks and sd
  robranks <- ranks %>% 
    tibble::as.tibble() %>%
    tidyr::gather(resample, rank, 2:31) %>%
    group_by(subset_size) %>%
    arrange(subset_size) %>%
    summarise_at(vars(rank), funs(mean, sd), na.rm = TRUE) %>%
    mutate(wvlt = as.numeric(gsub("rflt_", "", subset_size))) %>% 
    as.data.frame()
  
  return <- list(perf, robranks)
  return(return)

}
tidy_out <- lapply(data, tidyres)

#extract performance  
perf <- lapply(tidy_out, "[[", 1) %>% rbindlist()
perf$Exp <- rep(c("2016", "2017", "2018"), each = 84) 

#plot performance profiles
profile <- ggplot(perf, aes(x = subset_size, y = mean, group = set, colour = set)) +
  ggtitle("A") +
  geom_point() + geom_line() + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=2.5) + 
  xlab("subset size") + ylab("RMSE") +
  scale_x_continuous(limits = c(-2.5, 32.5)) +
  facet_wrap(~Exp, nrow = 1, ncol = 3) +
  theme_bw() +
  theme(legend.title=element_blank(),
        plot.title = element_text(size=15, face="bold"),
        legend.position  = c(0.85, 0.85))
  
#results per resample,  for each Experiement
tidyres2 <- function(d) {
  
  vars <- lapply(d,"[[", 1) %>%
    lapply(tail, 12) %>%
    dplyr::bind_rows() %>% 
    droplevels()
  
  #top features for each Resample
  glance <- lapply(d,"[[", 1) %>%
    lapply(tail, 12) %>%
    dplyr::bind_cols() %>% 
    dplyr::select(contains("var")) %>% 
    droplevels()
  names(glance) <- paste("Resample", 1:30, sep = "")
  
  #frequency of the unique features across Resamples
  table <- table(vars$var) %>% data.frame %>% 
    filter(Freq > 2) %>% 
    mutate(wvlt = gsub("rflt_", "", Var1) %>% as.numeric()) %>% 
    arrange(wvlt)

  return <- list(glance, table)
  return(return)
  
}
  
vars <- lapply(data, tidyres2)

resamples <- lapply(vars, "[[", 1) %>% rbindlist()
resamples$Exp <- rep(c("FPWW012", "FPWW018", "FPWW022"), each = 12)
resamples <- resamples %>% 
  dplyr::select(Exp, everything()) %>% 
  tibble::as.tibble()

  write.csv(resamples, "O:/Projects/KP0011/1/Analysis/rfe_cubist/Final_results/topfeat_resamples.csv",
            row.names = FALSE)

Freq <- lapply(vars, "[[", 2) %>% rbindlist()
Freq$Exp <- c(rep("2016", 36), rep("2017", 32), rep("2018", 33))
Freq <- Freq %>% 
  dplyr::select(Exp, everything()) %>% 
  tibble::as.tibble()

  write.csv(Freq, "O:/Projects/KP0011/1/Analysis/rfe_cubist/Final_results/topfeat_Freq.csv",
            row.names = FALSE)
  
Freq$Freq <- as.numeric((Freq$Freq/30))

#re-load mean reflectance for a senescence score of 8
dat <-  readRDS("O:/Projects/KP0011/1/Analysis/multvarmod/SpectralData/smth_avg_rflt.rds")

data_long <- dat %>%
  filter(Exp == "FPWW022") %>%
  select(SnsCnp, rflt_355:rflt_2399) %>%
  tidyr::gather(wvlt, reflectance, contains("rflt_"), factor_key = TRUE) %>%
  group_by(SnsCnp, wvlt) %>%
  summarise(mean_rflt = mean(reflectance)) %>% 
  filter(SnsCnp %in% c(8))

data_long$wvlt <- as.numeric(gsub("rflt_", "", data_long$wvlt)) 

data_long <- as.data.frame(data_long)
data_long$SnsCnp <- as.factor(data_long$SnsCnp)
data_long$range <- ifelse(data_long$wvlt < 1350, "visnir", ifelse(data_long$wvlt < 1781, "swir1", "swir2"))

data <- right_join(Freq, data_long, by = "wvlt")
data$cat <- ifelse(data$SnsCnp == 0, 1, ifelse(data$SnsCnp == 5, 2, 3))

freq_plot <- ggplot() +   
  ggtitle("B") +
  scale_shape_identity() +
  geom_point(data = data[!is.na(data$Exp),], aes(x = wvlt, y = Freq, group = Exp, color = Exp), size = 1) +
  geom_line(data = data, aes(x = wvlt, y = mean_rflt, group = cat), col = "grey", size = 1) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0.2, 0.4, 0.6, 0.8, 1.0),
                     sec.axis = sec_axis(~.*1, name = "Reflectance")) +
  scale_x_continuous(limits = c(500, 900)) +
  xlab("wavelength (nm)") + ylab("Frequency") +
  theme_bw() +
  geom_abline(slope = 0, intercept = 0.1, lty = 2) +
  theme(legend.title=element_blank(), 
        legend.position = c(0.9, 0.8),
        plot.title = element_text(size=15, face="bold"))

tiff("O:/Projects/KP0011/1/Figures/rfe_snscnp1.tiff", width = 10, height = 4, units = 'in', res = 300)
grid.arrange(freq_plot, profile,
             widths = c(1, 1),
             layout_matrix = rbind(c(2, 1)))
dev.off()

freq_plot <- ggplot(Freq) +   
  ggtitle("A") +
  scale_shape_identity()+
  geom_point(aes(x = wvlt, y = Freq, group = Exp, color = Exp), size = 1) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0.2, 0.4, 0.6, 0.8, 1.0)) +
  xlab("wavelength (nm)") + ylab("Frequency") +
  theme_bw() +
  geom_abline(slope = 0, intercept = 0.1, lty = 2) +
  theme(legend.title=element_blank(), 
        plot.title = element_text( size=20, face="bold"))
  
  



