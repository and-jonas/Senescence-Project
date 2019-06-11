#====================================================================================== -

#HEADER ----

# AUTHOR: Jonas Anderegg

# Prepare spectral Index dataset
# Extract error and senescence dynamics parameters from SI and scorings
# Extract correlations, mean errors and SI bias

#====================================================================================== -

.libPaths("T:/R3UserLibs")

library(tidyverse)
library(gridExtra)

source("O:/Projects/KP0011/3/Project/Septoria-Project/R/Utils/001_spectra_Utils.R")
source("O:/Projects/KP0011/1/Senescence-Project/R/Utils/001_A_spectra_utils.R")

#====================================================================================== -

###################################################################################### -
#COMPARISON PSRI vs. NDVI vs. SnsCnp ----
###################################################################################### -

#get data
data <- readRDS("O:/Projects/KP0011/1/Analysis/SI/ind_sc_d_posthead_complete.rds") %>% 
  dplyr::select(1:17, SI_NDVI_nb_ASD, SI_PSRI_r, SnsCnp) %>% 
  filter(Exp == "FPWW022") %>% 
  rowwise() %>% 
  mutate(we_GDDAH = ifelse(is.na(grading_GDDAH), meas_GDDAH, grading_GDDAH))
  

#create list of plots
list <- split(data, data$Plot_ID)
plist <- lapply(list, function(data) {
  d <- tidyr::gather(data, trait, value, 18:20)
  ggplot() + 
    geom_point(data = d, aes(x = we_GDDAH, y = value, shape = trait, col = trait), size = 0.3) +
    geom_line(data = d[!is.na(d$value),], aes(x = we_GDDAH, y = value, shape = trait, col = trait), size = 0.3) +
    ggtitle(paste(unique(data$Plot_ID))) +
    xlab("Thermal time after heading (°C days)") +
    ylab("Visual senescence score") +   
    geom_abline(intercept = 8, slope = 0, lty = 2) +
    geom_abline(intercept = 5, slope = 0, lty = 2) +
    geom_abline(intercept = 2, slope = 0, lty = 2) +
    theme_classic(base_size = 5)
  })

# plot to pdf
pdf("O:/Projects/KP0011/1/Figures/try.pdf", 7, 5)
for (i in seq(1, length(plist), 6)) {
  grid.arrange(grobs=plist[i:(i+5)], 
               ncol=3)
}
dev.off()

##done

###################################################################################### -
# > Create graph for contrasting plots ----
###################################################################################### -

d <- data %>% filter(Plot_ID %in% c("FPWW0220055")) %>% 
  tidyr::gather(Trait, value, 18:20)

tiff("O:/Projects/KP0011/1/Figures/time_series/New time series/FPWW0220055.tiff", width = 6, height = 3, units = 'in', res = 300)

ggplot() + 
  geom_line(data = d[!is.na(d$value),], aes(x = we_GDDAH, y = value, shape = Trait), size = 0.5) +
  theme(legend.position = "none") +
  geom_point(data = d, aes(x = we_GDDAH, y = value, shape = Trait), size = 2) +
  xlab("Thermal time after heading (°C days)") +
  ylab("Visual scoring \n Scaled SI") +   
  scale_shape_discrete(solid = FALSE,
                       name="Trait",
                       breaks=c("pred", "SI_NDVI_nb_ASD", "SI_PSRI_r", "SnsCnp"),
                       labels=c("PLS prediction", "NDVI", "PSRI", "Visual Canopy Score")) +
  scale_y_continuous(breaks = seq(0,10,2), limits = c(-0.5, 10)) +
  geom_abline(intercept = 8, slope = 0, lty = 1, col = "grey") +
  geom_abline(intercept = 5, slope = 0, lty = 1, col = "grey") +
  geom_abline(intercept = 2, slope = 0, lty = 1, col = "grey") +
  annotate("text", label = "A", x = 185, y = -0.5, size = 5, hjust = 0.5) +
  annotate("text", label = "B", x = 390, y = -0.5, size = 5, hjust = 0.5) +
  annotate("text", label = "C", x = 470, y = -0.5, size = 5, hjust = 0.5) +
  annotate("text", label = "D", x = 655, y = -0.5, size = 5, hjust = 0.5) +
  annotate("text", label = "E", x = 775, y = -0.5, size = 5, hjust = 0.5) +
  annotate("text", label = "F", x = 910, y = -0.5, size = 5, hjust = 0.5) +
  scale_color_brewer(name = "", palette = "Paired") +
  theme_classic(base_size = 12)

dev.off()

###################################################################################### -

#for agritechday

pdf("O:/Projects/KP0011/1/Figures/time_series/New time series/FPWW0220055_agritech.pdf", width = 6, height = 4)

ggplot() + 
  geom_line(data = d[!is.na(d$value),], aes(x = we_GDDAH, y = value, shape = Trait, color = Trait), size = 0.5) +
  geom_point(data = d, aes(x = we_GDDAH, y = value, shape = Trait), size = 2) +
  xlab("Temperatursumme") +
  ylab("Visuelle Bonitur \n Spektralindex") +  
  scale_color_discrete(name="Trait",
                       breaks=c("SI_NDVI_nb_ASD", "SI_PSRI_r", "SnsCnp"),
                       labels=c("NDVI", "PSRI", "Visuelle Bonitur")) +
  scale_shape_discrete(solid = FALSE,
                       name="Trait",
                       breaks=c("SI_NDVI_nb_ASD", "SI_PSRI_r", "SnsCnp"),
                       labels=c("NDVI", "PSRI", "Visuelle Bonitur")) +
  scale_y_continuous(breaks = seq(0,10,2), limits = c(-0.5, 10)) +
  geom_abline(intercept = 8, slope = 0, lty = 1, col = "grey") +
  geom_abline(intercept = 5, slope = 0, lty = 1, col = "grey") +
  geom_abline(intercept = 2, slope = 0, lty = 1, col = "grey") +
  annotate("text", label = "A", x = 185, y = -0.5, size = 5, hjust = 0.5) +
  annotate("text", label = "B", x = 390, y = -0.5, size = 5, hjust = 0.5) +
  annotate("text", label = "C", x = 470, y = -0.5, size = 5, hjust = 0.5) +
  annotate("text", label = "D", x = 655, y = -0.5, size = 5, hjust = 0.5) +
  annotate("text", label = "E", x = 775, y = -0.5, size = 5, hjust = 0.5) +
  annotate("text", label = "F", x = 910, y = -0.5, size = 5, hjust = 0.5) +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

dev.off()

###################################################################################### -
# MEAN RFLT FOR GIVEN SNS SCORING ----
###################################################################################### -

dat <-  readRDS("O:/Projects/KP0011/1/Analysis/multvarmod/SpectralData/smth_avg_rflt.rds")

data_long <- dat %>%
  filter(Exp == "FPWW022") %>%
  select(SnsCnp, rflt_355:rflt_2399) %>%
  tidyr::gather(wvlt, reflectance, contains("rflt_"), factor_key = TRUE) %>%
  group_by(SnsCnp, wvlt) %>%
  summarise(mean_rflt = mean(reflectance)) %>% 
  filter(SnsCnp %in% c(0, 2, 4, 6, 8, 10)) %>% 
  rename("Scoring" = "SnsCnp")

data_long$wvlt <- as.numeric(gsub("rflt_", "", data_long$wvlt)) 

data_long <- as.data.frame(data_long)
data_long$Scoring <- as.numeric(data_long$Scoring)
data_long$range <- ifelse(data_long$wvlt < 1350, "visnir", ifelse(data_long$wvlt < 1781, "swir1", "swir2"))

p0 <- ggplot(data_long) +
  geom_line(aes(x = wvlt, y = mean_rflt, group = interaction(Scoring, range) , color = Scoring), size = 0.75) +
  xlab("wavelength (nm)") + ylab("mean reflectance") +
  scale_x_continuous(breaks = c(500, 677, 750, 800, seq(1000, 2250, by = 250)), expand = c(0.01, 0.01)) +
  geom_vline(xintercept = 500) + geom_vline(xintercept = 677) +
  geom_vline(xintercept = 750) + geom_vline(xintercept = 800) +
  scale_color_gradient(low = "tan4", high = "darkgreen") +
  ggtitle("A") +
  theme(axis.line = element_line(colour = "black"),
        legend.position=c(0.9,0.7),
        panel.border = element_rect(fill = NA),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 15, angle = 90),
        axis.title.x = element_text(size = 15, angle = 0),
        axis.text.x = element_text(size = 12.5, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 12.5),
        legend.text = element_text(size = 12.5),
        legend.title = element_text(size = 15),
        plot.title = element_text( size=20, face="bold"))

###################################################################################### -
# COR WVLT ~ SNSCNP ----
###################################################################################### -

# load custom functions
library(stringr)
source("T:/PhD/Scripts/Functions/my_functions.R")
source("T:/PhD/Scripts/Grafiken/myPairs.R")

d_spctr <- readRDS("O:/Projects/KP0011/1/Analysis/multvarmod/SpectralData/smth_avg_rflt.rds")

# Drop measurement days for which no corresponding ratings are available
# Drop Plots which were measured only once during the senescence phase

dates <- c("2016-05-26", "2016-06-22", "2017-06-13", "2017-05-29")
plots <- paste("FPWW0180", c(491:594), sep = "")

spc <- d_spctr %>% 
  filter(!meas_date %in% as.Date(dates)) %>% 
  filter(!Plot_ID %in% plots) %>%
  droplevels()

#reorder columns
dat <- spc %>% select(heading_DAS:grading_GDDAS, 1:12, SnsFl0:SnsCnp, everything())
#remove rows containing missing data
dat <- dat[complete.cases(dat), ]

###################################################################################### -
# > across years ----

# subset years
dat$Exp <- stringr::str_sub(dat$Plot_ID, 1, 7)
d16 <- dat[dat$Exp == "FPWW012",]
d17 <- dat[dat$Exp == "FPWW018",]
d18 <- dat[dat$Exp == "FPWW022",]

#get correlations
cor_all <- -(as.vector(cor(as.matrix(dat[,19]), as.matrix(dat[,-c(1:19)]), method = "spearman")))
cor_16 <- -(as.vector(cor(as.matrix(d16[,19]), as.matrix(d16[,-c(1:19)]), method = "spearman")))
cor_17 <- -(as.vector(cor(as.matrix(d17[,19]), as.matrix(d17[,-c(1:19)]), method = "spearman")))
cor_18 <- -(as.vector(cor(as.matrix(d18[,19]), as.matrix(d18[,-c(1:19)]), method = "spearman")))

wvlt <- as.numeric(gsub("_", "", stringr::str_sub(names(dat)[-c(1:19)], -4, -1)))
COR <- as.data.frame(cbind(wvlt, cor_all, cor_16, cor_17, cor_18))

data_long <- tidyr::gather(COR, Exp, cor, cor_16:cor_18, factor_key = TRUE) %>% 
  mutate(range = ifelse(wvlt <= 1349, 1, ifelse(wvlt <= 1780, 2, 3)) %>% as.factor())

#create first part of the plot
p1 <- ggplot(data_long) + 
  geom_line(aes(x = wvlt, y = cor, group = interaction(Exp, range), color = Exp), size = 0.6) + 
  geom_vline(xintercept = 500) + geom_vline(xintercept = 677) +
  geom_vline(xintercept = 750) + geom_vline(xintercept = 800) +
  xlab("wavelength (nm)") + ylab("Correlation coefficient") +
  scale_x_continuous(breaks = seq(0,2500,500), limits = c(350, 2550), expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = seq(-1,1,0.2), limits = c(-1, 1)) +
  scale_color_discrete(name = "Year",
                       labels = c("2016", "2017", "2018")) +
  geom_abline(slope = 0, intercept = 0) +
  ggtitle("B") +
  theme(axis.line = element_line(colour = "black"),
        legend.position=c(0.8,0.2),
        panel.border = element_rect(fill = NA),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 15, angle = 90),
        axis.title.x = element_text(size = 15, angle = 0),
        axis.text.x = element_text(size = 12.5, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 12.5),
        legend.text = element_text(size = 12.5),
        legend.title = element_text(size = 15),
        plot.title = element_text( size=20, face="bold"))

##################################################################################### -

### > across classes ----

# subset classes
dat_Class1 <- dat[dat$SnsCnp < 4,]
dat_Class2 <- dat[dat$SnsCnp > 3 & dat$SnsCnp < 7,]
dat_Class3 <- dat[dat$SnsCnp > 7,]

#get correlations
late <- -c(as.vector(cor(as.matrix(dat_Class1[,19]), as.matrix(dat_Class1[,-c(1:19)]), method = "spearman")))
mid <- -c(as.vector(cor(as.matrix(dat_Class2[,19]), as.matrix(dat_Class2[,-c(1:19)]), method = "spearman")))
early <- -c(as.vector(cor(as.matrix(dat_Class3[,19]), as.matrix(dat_Class3[,-c(1:19)]), method = "spearman")))

max_early <- max(early, na.rm = TRUE)
max_mid <- max(mid, na.rm = TRUE)
max_late <- max(late, na.rm = TRUE)

wvlt <- as.numeric(gsub("_", "", stringr::str_sub(names(dat)[-c(1:19)], -4, -1)))
COR <- as.data.frame(cbind(wvlt, early, mid, late))

data_long <- tidyr::gather(COR, sub, cor, early:late, factor_key = TRUE) %>% 
  mutate(range = ifelse(wvlt <= 1349, 1, ifelse(wvlt <= 1780, 2, 3)) %>% as.factor())

#create second part of the plot
p2 <- ggplot(data_long) + 
  geom_line(aes(x = wvlt, y = cor, group = interaction(sub, range), col = sub), size = 0.6) + 
  geom_vline(xintercept = 500) + geom_vline(xintercept = 677) +
  geom_vline(xintercept = 750) + geom_vline(xintercept = 800) +
  scale_color_manual(values=c("green4", "darkolivegreen3", "tan4"), 
                     name = "Phase", 
                     labels = c("early", "mid", "late")) +
  xlab("wavelength (nm)") + ylab("Correlation coefficient") +
  scale_x_continuous(breaks = seq(0,2500,500), limits = c(350, 2550), expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = seq(-1,1,0.2), limits = c(-1, 1)) +
  geom_abline(slope = 0, intercept = 0) +
  ggtitle("C") +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill = NA),
        legend.position=c(0.8,0.2),
        panel.background = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 15, angle = 0),
        axis.text.x = element_text(size = 12.5, angle = 45, vjust = 0.5),
        axis.text.y = element_blank(),
        legend.text = element_text(size = 12.5),
        legend.title = element_text(size = 15),
        plot.title = element_text( size=20, face="bold"))

#combine all parts

pdf("O:/Projects/KP0011/1/Figures/spec_corrs.pdf", width = 10, height = 10)
gridExtra::grid.arrange(p0, p1, p2, widths = c(1.15,1),
                        layout_matrix = rbind(c(1),
                                              c(2, 3)))
dev.off()

###################################################################################### -
# WEATHER CONDITIONS ----
###################################################################################### -


library(dplyr)
library(ggplot2)

#Load Temperature data

load(file = "T:/PhD/DataAnalysis/Climate/MeteoData15_17.RDA")
WS15 <- c("2014-10-20 12:00:00", "2015-08-03 12:00:00")
WS16 <- c("2015-10-13 12:00:00", "2016-07-27 12:00:00")
WS17 <- c("2016-10-31 12:00:00", "2017-07-18 12:00:00")

meteo_WS15 <- filter(MeteoData15_17,Date >= WS15[1] & Date <= WS15[2]) %>% 
  mutate(.,GDD = cumsum(Tsum)) %>%
  select(Date, Tmean, Tmax)
meteo_WS16 <- filter(MeteoData15_17,Date >= WS16[1] & Date <= WS16[2]) %>% 
  mutate(.,GDD = cumsum(Tsum))%>%
  select(Date, Tmean, Tmax)
meteo_WS17 <- filter(MeteoData15_17,Date >= WS17[1] & Date <= WS17[2]) %>% 
  mutate(.,GDD = cumsum(Tsum))%>%
  select(Date, Tmean, Tmax)

meteo_WS18 <- readxl::read_xlsx("O:/Projects/KP0011/1/Data preparation/Data/raw/GDDAH_2018.xlsx") %>%
  select(day, mean, max)

names(meteo_WS18) <- c("Date", "Tmean", "Tmax")

# Only data from main growing period

meteo_WS16 <- meteo_WS16[meteo_WS16$Date > "2016-02-29 12:00:00",]
meteo_WS17 <- meteo_WS17[meteo_WS17$Date > "2017-02-28 12:00:00",]
meteo_WS18 <- meteo_WS18[meteo_WS18$Date > "2018-02-28 01:00:00",]

WS <- rbind(meteo_WS16, meteo_WS17, meteo_WS18)

# Load rainfall data

d16 <- read.csv("O:/Projects/KP0011/1/Data preparation/Data/raw/rainfall_2016.csv", sep = ";", header = FALSE)[-c(1:3),]
d17 <- read.csv("O:/Projects/KP0011/1/Data preparation/Data/raw/rainfall_2017.csv", sep = ";", header = FALSE)[-c(1:3),]
d18 <- read.csv("O:/Projects/KP0011/1/Data preparation/Data/raw/rainfall_2018.csv", sep = ";", header = FALSE)[-c(1:3),]

# Only data from main growing period

d16$V1 <- as.Date(d16$V1, format = "%d.%m.%Y")
d17$V1 <- as.Date(d17$V1, format = "%d.%m.%Y")
d18$V1 <- as.Date(d18$V1, format = "%d.%m.%Y")

d16 <- d16[d16$V1 > "2016-02-29",]
d17 <- d17[d17$V1 > "2017-02-28",]
d18 <- d18[d18$V1 > "2018-02-28",]


d <- rbind(d16, d17, d18)

d$V2 <- as.numeric(as.character(d$V2))

names(d) <- c("Date2", "Rainfall")
d <- d %>% filter(Date2 < "2018-07-21") %>%
  filter(Date2 != "2018-03-25") #somehow, this day is missing from temperature data...

WS <- cbind(WS, d)

WS$Year <- format(WS$Date, "%Y")

#load heading data

head16 <- read.csv("O:/Projects/KP0011/1/Data preparation/Data/out/Heading_FPWW012_DAS.csv")
head17 <- read.csv("O:/Projects/KP0011/1/Data preparation/Data/out/Heading_FPWW018_DAS.csv")
head18 <- read.csv("O:/Projects/KP0011/1/Data preparation/Data/out/Heading_FPWW022_DAS.csv", sep = ";")

head16$heading_date <- as.Date(head16$heading_date)
head17$heading_date <- as.Date(head17$heading_date)
head18$heading_date <- as.Date(head18$heading_date)

head <- rbind(head16, head17, head18)

head2 <- head %>% 
  mutate(Exp = str_sub(Plot_ID, 1, 7)) %>% 
  group_by(Exp) %>% 
  summarise(mean = mean(heading_date, na.rm = TRUE))

table(head$HeadingDAS)

count <- as.data.frame(table(head$heading_date), cut.names = FALSE)
colnames(count) <- c("Date", "head_freq")
count$Date <- as.Date(count$Date)

WS <- full_join(WS, count, by = c("Date2" = "Date"))

WS[is.na(WS$head_freq), "head_freq"] <- 0

WS$head_freq <- -WS$head_freq


#   p <- ggplot(WS, aes(Date2, head_freq))
#   p + geom_bar(stat = "identity") +
#     scale_x_date(date_breaks = "2 months") +
#     labs(x = "Date", y = "Frequency") +  
#     facet_wrap(~Year, scales = "free_x") + 
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           panel.background = element_blank(), axis.line = element_line(colour = "black"))

plot <- ggplot(WS) + 
  geom_line(aes(Date2, Tmean)) +
  geom_line(aes(Date2, Tmax), col = "red") +
  scale_x_date(date_breaks = "1 month") +
  geom_bar(aes(Date2, Rainfall/6), stat = "identity", col = "blue", fill = "blue") +    
  # geom_bar(stat = "identity", position = "identity", aes(Date2, head_freq)) +
  labs(x = "Date", y = "Temperature (°C)") +  
  scale_y_continuous(sec.axis = sec_axis(~.*6, name = "Rainfall (mm/m2)"), limits = c(-2,35)) +
  facet_wrap(~Year, scales = "free_x") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(size = 12.5, angle = 45, vjust = 1, hjust = 1))

pdf("O:/Projects/KP0011/weather_conditions.pdf", width = 10, height = 5)

print(plot)

dev.off()

###################################################################################### -
# MEAN INDEX DYNAMICS ----
###################################################################################### -

#Load data
spc_d <- readRDS("T:/PhD/DataAnalysis/FPWW012_FPWW018/Data/ind_d_posthead_complete.rds")[c(2, 10, 11, 13:146)]
sen_d <- readRDS("T:/PhD/DataAnalysis/FPWW012_FPWW018/Data/sen_d_clean.rds")

#Join spectral and scoring data
#and select only FPWW012 data
data <- f_match_join(spc_d, sen_d) %>% 
  dplyr::filter(Exp == "FPWW012")

#define SIs for which to create the plot
SI <- c("SI_VARIgreen", "SI_mND705", "SI_NDVI_nb_ASD", 
        "SI_780_740", "SI_PSRI_r", "SI_NDWI1",
        "SnsCnp") %>% as.list()

#> linintpol SI per plot ----
interp_si <- function(index, data){
  
  int <- list()
  print(index)
  
  #interpolate for each Plot -
  for (i in unique(data$Plot_ID)){
    tryCatch({
      int[[i]] <- approx(data[data$Plot_ID == i,]$meas_GDDAH, data[data$Plot_ID == i,][,index],
                         xout = seq(round(min(data[data$Plot_ID == i,]$meas_GDDAH, na.rm = TRUE), 0), 
                                    round(max(data[data$Plot_ID == i,]$meas_GDDAH, na.rm = TRUE), 0), 1))
    }, error=function(e){cat}) #("ERROR :",conditionMessage(e), "\n")
  }
  return(int)
}
d <- lapply(SI, interp_si, data = data)
names(d) <- unlist(SI)

#common range of GDDAH ----
complete_list <- function(data){
  diff_start <- data$x[1]-0
  diff_end <- 1400 - data$x[length(data$x)]
  data$x <- c(0:1400)
  data$y <- c(rep(NA, diff_start), data$y, rep(NA, diff_end))
  return(data)
}

do_it <- function(data){
  
  d_new <- data
  d_compl <- lapply(d_new, complete_list)
  
  #extract GDDAH
  GDDAH <- d_compl[[1]]$x
  
  #extract preditions
  d_preds <- lapply(d_compl, "[[", "y")
  
  #calcualte mean index value at each GDDAH across all Plots
  means <- plyr::aaply(plyr::laply(d_preds, as.matrix), 2, mean, na.rm = TRUE)
  #calculate standard deviation at each GDDAH across all plots
  SD <- plyr::aaply(plyr::laply(d_preds, as.matrix), 2, stats::sd, na.rm = TRUE)
  
  #build dataframe
  means <- as.data.frame(means)
  SD <- as.data.frame(SD)
  all <- cbind(GDDAH, means, SD)
  all$upper <- all$means + all$SD
  all$lower <- all$means - all$SD
  all <- all[c(1:2, 4:5)]
  return(all)
  
}

ALL <- lapply(d, do_it)

#adjust SI names
names(ALL) <- c("VARIgreen", "mND705", "NDVI", "R780/R740", "PSRI", "NDWI", "SnsCnp")

## PLOT MEAN INDEX DYNAMICS -
#################################################################### -

data_long <- do.call("rbind", ALL)

data_long$Index <- sapply(strsplit(rownames(data_long), ".", fixed = TRUE), "[", 1)
d_cnp <- data_long[data_long$Index == "SnsCnp",]

create_gg <- function(data, Index){
  
  d1 <- data[data$Index == Index,] 
  d_Cnp <-data[data$Index == "SnsCnp",]
  
  g1 <- ggplot(d1) + 
    stat_smooth(data = d1, aes(GDDAH, means), col = "black") +
    geom_smooth(data = d1, aes(GDDAH, upper), col = "black") +
    geom_smooth(data = d1, aes(GDDAH, lower), col = "black") +
    geom_smooth(data = d_cnp, aes(GDDAH, means), col = "darkgreen") +
    geom_smooth(data = d_cnp, aes(GDDAH, upper), col = "darkgreen") +
    geom_smooth(data = d_cnp, aes(GDDAH, lower), col = "darkgreen") +
    geom_abline(intercept = 8, slope = 0, lty = 2) +
    geom_abline(intercept = 5, slope = 0, lty = 2) +
    geom_abline(intercept = 2, slope = 0, lty = 2) +
    xlab("Thermal time after heading (°C days)") +
    ylab("Visual scoring \n Scaled SI") +   
    scale_y_continuous(breaks = seq(0,10,2), limits = c(-1, 11)) +
    theme(axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA),
          legend.position=c(0.8,0.2),
          panel.background = element_blank(),
          axis.title.x = element_text(size = 15, angle = 0),
          axis.text.x = element_text(size = 12.5, angle = 0, vjust = 0.5),
          legend.text = element_text(size = 12.5),
          legend.title = element_text(size = 15),
          plot.title = element_text(size=20, face="bold"))
    
  # build plot object for rendering 
  gg1 <- ggplot_build(g1)
  
  # extract data for the loess lines from the 'data' slot
  df2 <- data.frame(x = gg1$data[[1]]$x,
                    ymean <- gg1$data[[1]]$y,
                    ymin = gg1$data[[2]]$y,
                    ymax = gg1$data[[3]]$y) 
  df3 <- data.frame(x = gg1$data[[4]]$x,
                    ymean <- gg1$data[[4]]$y,
                    ymin = gg1$data[[5]]$y,
                    ymax = gg1$data[[6]]$y)
  
  # use the loess data to add the 'ribbon' to plot 
  GG <- g1 +
    geom_ribbon(data = df2, aes(x = x, ymin = ymin, ymax = ymax),
                fill = "grey", alpha = 0.5) +
    geom_ribbon(data = df3, aes(x = x, ymin = ymin, ymax = ymax),
                fill = "darkgreen", alpha = 0.1)
  
  # modify SI names
  
  facet <- unique(GG$data$Index)
  
  GG <- GG + 
    theme_bw(base_size = 12) %+replace%
    theme(panel.grid =  element_blank(),
          strip.text.x = element_text(size = 12, colour = "black", angle = 0)) +
    facet_wrap(~paste(facet))
  
  return(GG)
  
}

plots <- lapply(as.list(names(ALL)), create_gg, data = data_long)

#x axis needs to be removed from top 4 plots
gg_mod1 <- function(plot) {
    plot + theme(axis.text.x = element_blank(),
                 axis.title.x = element_blank())
}

plot_list <- plots[1:4]
red_plots <- lapply(plot_list, gg_mod1)

#reassemble plot list
plots <- c(red_plots, plots[5:6])

#y axis needs to be removed from 3 left right
gg_mod2 <- function(plot) {
  plot + theme(axis.text.y = element_blank(),
               axis.title.y = element_blank())
}

plot_list <- plots[c(2, 4, 6)]
red_plots <- lapply(plot_list, gg_mod2)

pdf("O:/Projects/KP0011/si.pdf", width = 9, height = 8)
gridExtra::grid.arrange(plots[[1]], red_plots[[1]],
                        plots[[3]], red_plots[[2]],
                        plots[[5]], red_plots[[3]], 
                        widths = c(10, 9),
                        heights = c(8, 8, 9),
                        layout_matrix = rbind(c(1, 2),
                                              c(3, 4), 
                                              c(5, 6)))
dev.off()

###################################################################################### -
# RFE SnsCnp ----
###################################################################################### -

# O:\Projects\KP0011\1\Analysis\rfe_cubist\Scripts\rfe_cubist_summary.R

###################################################################################### -
# Multvarmod validation  ----
###################################################################################### -

# > cubist ----

#data for all cubist models
dir <- "O:/Projects/KP0011/1/Analysis/multvarmod/Results/cubist/SnsCnp/fullsample/center_scale"
preProc <- "center_scale"
files <- as.list(list.files(paste(dir), pattern = ".rds"))
setwd(dir)
data <- lapply(files, readRDS)
names(data) <- sapply(strsplit(unlist(files), "\\."), "[[", 1)

  # #ORIGINAL GRAPH
  # 
  #   #only smth avg rflt bin3
  #   d <- data[[2]]
  #   D <- unlist(d, recursive = FALSE)
  # 
  #   data <- D[[7]]
  # 
  #   data$data$obsf <- as.factor(data$data$obs)
  #   
  #   Plot <- ggplot() +
  #     geom_boxplot(aes(x = data$data$obsf, y = data$data$pred),
  #                  size = 0.5, alpha = 1) +
  #     stat_summary(mapping = aes(x = data$data$obs + 1, y = data$data$pred), fun.y = "median", 
  #                  colour = "red", size = 3, geom = "point", shape=18) +
  #     scale_y_continuous(breaks = seq(0,10,2), limits = c(-2, 12)) +
  #     geom_abline(intercept = -1, slope = 1, color = "red", size = 1, linetype="dashed") +
  #     xlab("observed") + ylab("predicted (10-fold CV)") +
  #     annotate("text", x = 9, y = 0, label= paste("RMSE = ", data$RMSE, sep = ""), size = 2.5) +
  #     annotate("text", x = 9, y = -0.5, label = paste("Coms = ", data$bestTune$committees, " ", 
  #                                                     "Nbrs = ", data$bestTune$neighbors, sep = ""), size = 2.5) +
  #     annotate("text", x = 9, y = -1, label = paste("nsamp = ", data$nsamp, sep = ""), size = 2.5) +
  #     ggtitle(paste("training = ", paste(data$training, collapse = ", "),
  #                   "\nvalidation = ", paste(data$validation, collapse = ", "),
  #                   paste("\ndata = ", i, sep = ""), paste("\npreProc = ", preProc, sep = ""), sep = "")) +
  #     theme_bw() +
  #     theme(plot.title = element_text(size = 7, face = "bold"))
  
#only smth avg rflt bin3
d <- data[[2]]
D <- unlist(d, recursive = FALSE)

#FPWW018 -> FPWW018 AND FPWW018 -> FPWW022
data1 <- D[[5]]$data %>% dplyr::select(-rowIndex)
data1$obsf <- as.factor(data1$obs)

data2 <- D[[7]]$data
data2$obsf <- as.factor(data2$obs)

data <- rbind(data1, data2)

plot_cubist <- ggplot(data) +
  ggtitle("B") +
  geom_boxplot(aes(x = obsf, y = pred, group = interaction(obsf, Exp), fill = Exp),
               size = 0.1, alpha = 1, width = 0.4) +
  scale_y_continuous(breaks = seq(0,10,2), limits = c(0, 10)) +
  geom_abline(intercept = -1, slope = 1, color = "grey", size = 1, linetype="dashed") +
  scale_fill_discrete(name="Validation",
                      breaks=c("FPWW018", "FPWW022"),
                      labels=c("Within-year \n(2017)", "Across-year \n(2018)")) +
  xlab("observed") + ylab("predicted (10-fold CV)") +
  annotate("text", x = 2, y = 8, label= paste("RMSE = 0.68 \n Coms = 5; Nbrs = 0 \n nsamp = 2555")
           , size = 4, col = "#F8766D") +
  annotate("text", x = 10, y = 2, label= paste("RMSE = 1.45 \n Coms = 5; Nbrs = 0 \n nsamp = 638")
           , size = 4, col = "#00BFC4") +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(face = "bold"))

##### -

# > PLSR ----

#data for all cubist models
dir <- "O:/Projects/KP0011/1/Analysis/multvarmod/Results/PLSR/SnsCnp/fullsample/center_scale"
preProc <- "center_scale"
files <- as.list(list.files(paste(dir), pattern = ".rds"))
setwd(dir)
data <- lapply(files, readRDS)
names(data) <- sapply(strsplit(unlist(files), "\\."), "[[", 1)

  # #ORIGINAL GRAPH
  # 
  #   #only smth avg rflt bin3
  #   d <- data[[4]]
  #   D <- unlist(d, recursive = FALSE)
  # 
  #   data <- D[[7]]
  # 
  #   data$data$obsf <- as.factor(data$data$obs)
  #   
  #   Plot <- ggplot() +
  #     geom_boxplot(aes(x = data$data$obsf, y = data$data$pred),
  #                  size = 0.5, alpha = 1) +
  #     stat_summary(mapping = aes(x = data$data$obs + 1, y = data$data$pred), fun.y = "median", colour = "red", size = 3, geom = "point", shape=18) +
  #     scale_y_continuous(breaks = seq(0,10,2), limits = c(-2, 12)) +
  #     geom_abline(intercept = -1, slope = 1, color = "red", size = 1, linetype="dashed") +
  #     xlab("observed") + ylab("predicted (5-times repeated 10-fold CV)") +
  #     annotate("text", x = 9, y = 0, label= paste("RMSE = ", data$RMSE, sep = ""), size = 2.5) +
  #     annotate("text", x = 9, y = -0.5, label = paste("ncomp = ", data$ncomp, sep = ""), size = 2.5) +
  #     annotate("text", x = 9, y = -1, label = paste("nsamp = ", data$nsamp, sep = ""), size = 2.5) +
  #     ggtitle(paste("training = ", paste(data$training, collapse = ", "),
  #                   "\nvalidation = ", paste(data$validation, collapse = ", "),
  #                   paste("\ndata = ", i, sep = ""), paste("\npreProc = ", preProc, sep = ""), sep = "")) +
  #     theme_bw(base_size = 15) +
  #     theme(plot.title = element_text(size = 7, face = "bold"))

#only smth avg rflt
d <- data[[4]]
D <- unlist(d, recursive = FALSE)

#FPWW018 -> FPWW018 AND FPWW018 -> FPWW022
data1 <- D[[5]]$data %>% dplyr::select(-rowIndex)
data1$obsf <- as.factor(data1$obs)

data2 <- D[[7]]$data
data2$obsf <- as.factor(data2$obs)

data <- rbind(data1, data2)

plot_plsr <- ggplot(data) +
  ggtitle("A") +
  geom_boxplot(aes(x = obsf, y = pred, group = interaction(obsf, Exp), fill = Exp),
               size = 0.1, alpha = 1, width = 0.4) +
  scale_y_continuous(breaks = seq(0,10,2), limits = c(0, 10)) +
  geom_abline(intercept = -1, slope = 1, color = "grey", size = 1, linetype="dashed") +
  xlab("observed") + ylab("predicted (5 times repeated 10-fold CV)") +
  annotate("text", x = 2, y = 8, label= paste("RMSE = 0.91 \n ncomp = 8 \n nsamp = 2555")
           , size = 4, col = "#F8766D") +
  annotate("text", x = 10, y = 2, label= paste("RMSE = 1.98 \n ncomp = 8 \n nsamp = 638")
           , size = 4, col = "#00BFC4") +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none")

# > combine ----

tiff("O:/Projects/KP0011/1/Figures/multvarmod_val.tiff", width = 15, height = 5, units = 'in', res = 300)
gridExtra::grid.arrange(plot_plsr, plot_cubist, widths = c(0.9, 1.1),
                        layout_matrix = rbind(c(1, 2)))
dev.off()

###################################################################################### -
# WVLT sensitivity  ----
###################################################################################### -

#load data
result <- readRDS("O:/Projects/KP0011/1/Analysis/Senescence dynamics/Data/Function_output/wvlt_search/SR_complete.rds")

#extract information from result
stat_data_exp_SnsCnp <- extract_inf(result, trait = "SnsCnp", out = "stat_data_exp")

#tidy up output
final <- stat_data_exp_SnsCnp %>% 
  dplyr::select(1, 3, contains("cor")) %>%
  mutate(wvlt = gsub("SI._R|_r", "", Index) %>% as.numeric(),
         Index = str_sub(Index, 1, 3)) %>% 
  tidyr::gather(par, corr, cor_onsen:cor_tsen) %>% 
  mutate(par = gsub("cor_", "", par)) %>% 
  filter(par != "tsen") %>% 
  mutate(SI = ifelse(Index == "SI0", "(R678-R500)/x", 
                     ifelse(Index == "SI1", "R678/x", 
                            ifelse(Index == "SI2", "x/R765", 
                                   "(R678-x)/R765"))))

p <- ggplot(final) +
  geom_line(aes(x = wvlt, y = corr, group = interaction(Exp, par), color = Exp, linetype = par)) +
  facet_wrap(~SI) +
  xlab("Wavelength (nm)") + ylab("Correlation coefficient") +
  # labs(linetype="Parameter", colour="Year") +
  scale_colour_discrete(name  ="Year",
                        breaks=c("FPWW012", "FPWW018", "FPWW022"),
                        labels=c("2016", "2017", "2018")) +
  scale_linetype_discrete(name="Parameter",
                      breaks=c("endsen", "midsen", "onsen"),
                      labels=c("Endsen", "Midsen", "Onsen")) +
  theme_bw(base_size = 15)

tiff("O:/Projects/KP0011/1/Figures/sensitivity.tiff", width = 10, height = 8, units = 'in', res = 300)
print(p)
dev.off()

###################################################################################### -
# RFE primary traits  ----
###################################################################################### -

### > Performance profiles ----

dir <- "O:/Projects/KP0011/1/Analysis/RFE_New"
setwd(dir)

files <- as.list(list.files(dir, pattern = ".rds"))

RMSE_test <- RMSE_train <- ranks <- list()

for(i in files){
  
  print(i)
  
  #load data
  data <- readRDS(i)
  
  #Tidy up
  allranks <- lapply(data, "[[", 1)
  rnk <- allranks %>%
    Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="var"), .)
  colnames(rnk) <- c("subset_size", paste("Resample", 1:30, sep = ""))
  
  ranks[[i]] <- rnk %>% 
    tibble::as.tibble() %>%
    gather(resample, rank, 2:31) %>%
    group_by(subset_size) %>%
    arrange(subset_size) %>%
    summarise_at(vars(rank), funs(mean, sd), na.rm = TRUE) %>%
    arrange(mean) %>%
    mutate(year = unlist(strsplit(i, "\\.|_"))[2],
           trait = unlist(strsplit(i, "\\.|_"))[1]) %>% 
    mutate(order = row_number())
  
  allRMSE <- lapply(data, "[[", 2) %>% lapply(as.data.frame)
  RMSE <- allRMSE %>% 
    lapply(function(d) {
      dd <- apply(d, 2, unlist) %>% as.data.frame()
    }) %>% 
    Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="subset_size"), .)
  
  RMSE_train[[i]] <- RMSE %>% 
    tibble::as.tibble() %>%
    select(1, contains("train")) %>% 
    `colnames<-`(c("subset_size", paste("Resample", 1:30, sep = ""))) %>% 
    tidyr::gather(resample, RMSE, 2:31) %>% 
    group_by(subset_size) %>% 
    arrange(subset_size) %>% 
    summarise_at(vars(RMSE), funs(mean, sd), na.rm = TRUE) %>% 
    mutate(year = unlist(strsplit(i, "\\.|_"))[2],
           trait = unlist(strsplit(i, "\\.|_"))[1],
           data = "train")

  RMSE_test[[i]] <- RMSE %>% 
    tibble::as.tibble() %>%
    select(1, contains("test")) %>% 
    `colnames<-`(c("subset_size", paste("Resample", 1:30, sep = ""))) %>% 
    tidyr::gather(resample, RMSE, 2:31) %>% 
    group_by(subset_size) %>% 
    arrange(subset_size) %>% 
    summarise_at(vars(RMSE), funs(mean, sd), na.rm = TRUE) %>% 
    mutate(year = unlist(strsplit(i, "\\.|_"))[2],
           trait = unlist(strsplit(i, "\\.|_"))[1],
           data = "test")

}

#final dataset for RMSE
all <- c(RMSE_test, RMSE_train)
all <- do.call("rbind", all) %>%  
  filter(data == "test") %>% 
  filter(subset_size < 51)

all$trait <- ifelse(all$trait == "gpc", "GPC", "GY")

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.5) # move them .05 to the left and right

p1 <- ggplot(all, aes(x=subset_size, y=mean, col = year)) + 
  ggtitle("A") +
  geom_point(position = pd) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=2, position = pd) +
  geom_line(position=pd, aes(group = year)) +
  facet_wrap(~trait, ncol = 2, scales = "free") +
  xlab("feature subset size") +
  ylab("mean test RMSE") + 
  theme_bw() %+replace%
  theme(plot.title = element_text(hjust = -0.05, size=15, face="bold"))


###################################################################################### -

### > Feature ranks ----

### > > GY ----

#final dataset for ranks:gy

order <- ranks[[3]][c(1, 6)]
ranks_16 <- ranks[[3]][1:6]
ranks_17 <- ranks[[4]][1:5]
ranks_17 <- full_join(ranks_17, order, by = "subset_size")
ranks_18 <- ranks[[5]][1:5]
ranks_18 <- full_join(ranks_18, order, by = "subset_size")

ranks <- rbind(ranks_16, ranks_17, ranks_18) %>% 
  arrange(year, trait, order) %>% 
  group_by(year) %>%
  #drop middle features for better graphs
  slice(-15:-88) %>% 
  mutate(trait = ifelse(trait == "gy", "GY", "GPC"))

ranks$order <- as.factor(ranks$order)

pd <- position_dodge(0.5) # move them .05 to the left and right
p2a <- ggplot(ranks, aes(x=order, y=mean, group = year, color = year)) + 
  ggtitle("") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1.5, position=pd) +
  geom_vline(aes(xintercept = 14.5), size = 1, linetype="dashed") +
  # geom_line(position=pd) +
  geom_point(position=pd) +
  ylab("Mean feature rank") +
  xlab("Feature")+
  # Add categories to axis
  scale_x_discrete(
    breaks = ranks$order,
    labels = ranks$subset_size,
    expand = c(0,0)
  ) +
  coord_flip() +
  theme_bw() %+replace%
  theme(axis.title.y =  element_blank(),
        plot.title = element_text(size=15, face="bold")) +
  facet_wrap(~trait)

### > > GPC ----

#re-run loop first

order <- ranks[[1]][c(1, 6)]
ranks_16 <- ranks[[1]][1:6]
ranks_17 <- ranks[[2]][1:5]
ranks_17 <- full_join(ranks_17, order, by = "subset_size")

ranks <- rbind(ranks_16, ranks_17) %>% 
  arrange(year, trait, order) %>% 
  group_by(year) %>%
  #drop middle features for better graphs
  slice(-15:-88) %>% 
  mutate(trait = ifelse(trait == "gy", "GY", "GPC"))

ranks$order <- as.factor(ranks$order)

pd <- position_dodge(0.5) # move them .05 to the left and right
p2b <- ggplot(ranks, aes(x=order, y=mean, group = year, color = year)) + 
  ggtitle("B") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1.5, position=pd) +
  geom_point(position=pd) +
  scale_color_manual(values = c("#F8766D", "#00BA38")) +
  ylab("Mean feature rank") +
  xlab("Feature") +
  geom_vline(aes(xintercept = 14.5), size = 1, linetype="dashed") +
  # Add categories to axis
  scale_x_discrete(
    breaks = ranks$order,
    labels = ranks$subset_size,
    expand = c(0,0)
  ) +
  coord_flip() +
  theme_bw() %+replace%
  theme(legend.position = "none",
        plot.title = element_text(hjust = -0.83, size=15, face="bold")) +
  facet_wrap(~trait)

#### arrange all plots

ggsave(
  filename = "O:/Projects/KP0011/1/Figures/rfe_GPC_GY.tiff",
  gridExtra::grid.arrange(p1, p2a, p2b, 
                        widths = c(1, 1.135),
                        heights = c(1, 1.5),
                        layout_matrix = rbind(c(1, 1),
                                              c(3, 2))),
  device = "tiff", 
  dpi = 300)

