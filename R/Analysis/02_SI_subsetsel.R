#====================================================================================== -

#HEADER ----

# AUTHOR: Jonas Anderegg

# Prepare subset of SI based on filtering criteria

#====================================================================================== -

.libPaths("T:/R3UserLibs")

library(dplyr)
library(tidyr)

#====================================================================================== -

# Plot SI dynamics for all calculated SI over all plots,
# Used to identify SI with a more or less monotoneous decreasing pattern over time,
# to help the selection of subset of SI

#load data
dat <- readRDS("O:/Projects/KP0011/1/Senescence-Project/Data/SI/SI_sc_posthead_compl.rds")

#list all SI
ind <- dplyr::select(dat, contains("SI")) %>% names()
dat <- dat[dat$Exp == "FPWW022",]

##CREATE PLOTS

pdf("O:/Projects/KP0011/1/Analysis/Figures/Sendyn_IND_linear_FPWW022_new.pdf", 
    onefile = TRUE, paper = "a4")

op <- par(mfrow = c(5,3),
          oma = c(3,2,0,0) + 1.5,
          mar = c(0,0,1,1) + 1.5)

for (i in ind){
  
  print(i)
  
  xvals <- split(dat$meas_GDDAH, dat$Plot_ID)
  yvals <- split(dat[[i]], dat$Plot_ID)
  
  plot(1:max(unlist(xvals), na.rm = TRUE), xlim = c(-100, 1000), ylim = (c(min(unlist(yvals), na.rm = TRUE),max(unlist(yvals), na.rm = TRUE))), type = "n", 
       xlab = "meas_GDDAH", ylab = i)
  
  mapply(lines, xvals, yvals, type="l")
  
  text(x = 50, y = min(unlist(yvals)) + 1, labels = paste(gsub("SI_", "", gsub("_r", "", i))), col = "red")
  
}

title(xlab = "meas_GDDAH",
      ylab = "Scaled value of SI",
      outer = TRUE, line = 1.5,
      cex.lab = 1.5)

par(op)

dev.off()

#====================================================================================== -

### STARTING FROM HERE: OLD !!!!!!

# load senescence parameter data
  pars_all <- readRDS("T:/PhD/DataAnalysis/FPWW012_FPWW018/ANOVA_ind/Data/pars_lin_Ind_wide.rds")
  
  #select based on FPWW012 data (both replicates measured, allowing to calculate w2)
  pars <- pars_all %>% filter(Exp == "FPWW012")

# select onsen 
  onsen <- select(pars_all, contains("onsen"))

# select based on heritabilities
  h2 <- read.csv("T:/PhD/DataAnalysis/FPWW012_FPWW018/ANOVA_ind/Results/Heritabilities_2016_all.csv") %>%
    filter(grepl('onsen', X)) %>% filter(H2_2016 > 0.6) %>% droplevels()

  traits <- as.character(h2$X)

  onsen <- onsen %>%
    select(one_of(traits))
  
# Drop indices to avoid very high correlations: 
  ## SI are dropped preferring nb over bb, SI with clear physiological meaning,
  ## SI developped in wheat or cereals

  drop <- c("SI_GNDVI_IRIS", "NDVI_bb", "NDVI_MODIS", "NDVI_nb_CASI", 
            "SI_PRI_10", "SI_PRI_6", "SI_PRI570", 
            "VOG2",
            "GLI",
            "VIgreen",
            "lai",
            "PSSR1",
            "RVI1",
            "RVI2",
            "PSND2", #somewhat lower h2
            "PSND3", #somewhat lower h2
            "780_740",
            "760_730",
            "780_550", #maize
            "CHLRE", #trees
            "MTCI", #trees
            "NDTI", #noisy, lower h2
            "CARRE", #trees, maize, soybean
            "CARG", #trees, maize, soybean
            "TGI_bb",
            "NDREI",
            "SR_dm",
            "ARVI",
            "VARI700",
            "PSSR2",
            "TCARI_r",
            "NHI_ALI",
            "RGR2",
            "NDMI",
            "LCI",
            "WI_NDVI",
            "NGRDI",
            "NDVI",
            "RRDI", #not monotoneous
            "NDWI1650",
            "WDRVI", #broad, MSR_rev (narrow) is used
            "CHLG",
            "CIG",
            "CIRE",
            "VOG1", 
            "VOG3",
            "GM",
            "PBI",
            "PSSR",
            "REM",
            "HI", #sugar beet health
            "MCARI_rev", #not monotonesous
            "TCARI_OSAVI_rev" #not monotoneous, noisy at start
)

dropnames <- onsen %>%
  select(matches(paste(drop, collapse = "|"))) %>% names()

cleaned <- onsen %>%
  select(-one_of(dropnames))

M <- cor(cleaned, use = "pairwise.complete.obs")
corrplot(M, method = "color")

onsen_pars <- names(cleaned)

zdf <- as.data.frame(as.table(M))
zdf

subset(zdf, abs(Freq) > 0.96)


##################################################################################################################  

# select midsen 
midsen <- select(pars_all, contains("midsen"))

# select based on heritabilities
h2 <- read.csv("T:/PhD/DataAnalysis/FPWW012_FPWW018/ANOVA_ind/Results/Heritabilities_2016_all.csv") %>%
  filter(grepl('midsen', X)) %>% filter(H2_2016 > 0.7) %>% droplevels()

traits <- as.character(h2$X)

midsen <- midsen %>%
  select(one_of(traits))

drop <- c("SI_GNDVI_IRIS", "NDVI_bb", "NDVI_MODIS", "NDVI_nb_CASI",
          "PSND1", "PSND3", "PSND2", "NDRE", #very high corr with mND705
          "SI_PRI_10", "SI_PRI_6", "SI_PRI570",
          "VOG2",
          "MSR_rev",
          "RGR2", "RGR", #very high corr with VARIgreen
          "lai",
          "RVI1",
          "NDWI2130",
          "780_740",
          "SR_dm",
          "NHI",
          "HI", #sugarbeat health
          "ARVI",
          "NGRDI",
          "CIRE",
          "NDMI",
          "LWVI2",
          "760_730",
          "780_700",
          "PSSR",
          "CHLRE",
          "VARI700",
          "SIWSI",
          "VIgreen",
          "OCAR", # much lower h2 for onsen
          "NDREI",
          "WI_NDVI",
          "REM",
          "SI_SR",
          "LCI2",
          "NDVI",
          "MTCI", #for trees
          "MCARI2",
          "NDWI1650", #slightly lower h2; both in onsen set
          "SI_OSAVI",
          "VOG1",
          "GLI", #highly unstable over time
          "TCARI_OSAVI_rev" #included in onsen dataset
)

dropnames <- midsen %>%
  select(matches(paste(drop, collapse = "|"))) %>% names()

cleaned <- midsen %>%
  select(-one_of(dropnames))

M <- cor(cleaned, use = "pairwise.complete.obs")
corrplot(M, method = "color")

midsen_pars <- names(cleaned)

##################################################################################################################  

# select endsen 
endsen <- select(pars, contains("endsen"))

# select based on heritabilities
h2 <- read.csv("T:/PhD/DataAnalysis/FPWW012_FPWW018/ANOVA_ind/Results/Heritabilities_2016_all.csv") %>%
  filter(grepl('endsen', X)) %>% filter(H2_2016 > 0.8) %>% droplevels()

traits <- as.character(h2$X)

endsen <- endsen %>%
  select(one_of(traits))

drop <- c("SI_GNDVI_IRIS", "NDVI_bb", "NDVI_MODIS", "NDVI_nb_CASI","NDVI",
          "SI_PRI_10", "SI_PRI_6", "SI_PRI570",
          "VOG2",
          "RGR2",
          "lai",
          "RVI1",
          "NDWI2130",
          "780_740",
          "SR_dm",
          "NHI",
          "ARVI",
          "NGRDI",
          "CIRE",
          "NDMI",
          "LWVI2",
          "760_730",
          "780_700",
          "1200", #quite noisy, compared to WI
          "900", #WI instead
          "PSSR",
          "CHLRE",
          "VARI700",
          "SIWSI",
          "VIgreen",
          "OCAR", # much lower h2 for onsen
          "NDREI",
          "WI_NDVI",
          "REM",
          "SI_SR",
          "PSND1", "PSND2", "PSND4", #1 and 4 used for onsen
          "SIPI", #already used, h2 somewhat lower
          "mND705", #already used, h2 somewhat lower
          "LCI2",
          "NDVI",
          "MTCI", #for trees
          "MCARI1",
          "MCARI2",
          "MSAVI",
          "SAVI",
          "NDWI1650", #slightly lower h2; both in onsen set
          "SI_OSAVI",
          "VOG1",
          "EVI", #vegetation(shrubs, grass etc.)
          "TVI", "MTVI1" #noisier at the beginning than MTVI2, otherwise equivalent
)

dropnames <- endsen %>%
  select(matches(paste(drop, collapse = "|"))) %>% names()

cleaned <- endsen %>%
  select(-one_of(dropnames))

M <- cor(cleaned, use = "pairwise.complete.obs")
corrplot(M, method = "color")

endsen_pars <- names(cleaned)


zdf <- as.data.frame(as.table(M))
zdf

subset(zdf, abs(Freq) > 0.96)

##################################################################################################################  

# select tsen 
tsen <- select(pars, contains("tsen"))

# select based on heritabilities
h2 <- read.csv("T:/PhD/DataAnalysis/FPWW012_FPWW018/ANOVA_ind/Results/Heritabilities_2016_all.csv") %>%
  filter(grepl('tsen', X)) %>% filter(H2_2016 > 0.5) %>% droplevels()

traits <- as.character(h2$X)

tsen <- tsen %>%
  select(one_of(traits))

drop <- c("NDVI_bb",
          "NDVI_MODIS",
          "NDVI_nb",
          "PSSR",
          "CIRE",
          "VOG1", "VOG2",
          "RVI1",
          "PSND1",
          "lai",
          "SI_SR",
          "LCI2",
          "REM",
          "ARVI",
          "CHLRE",
          "MTCI",
          "NDREI"
)
          

dropnames <- tsen %>%
  select(matches(paste(drop, collapse = "|"))) %>% names()

cleaned <- tsen %>%
  select(-one_of(dropnames))

M <- cor(cleaned, use = "pairwise.complete.obs")
corrplot(M, method = "color")

tsen_pars <- names(cleaned)

##################################################################################################################  

# Combine

all_pars <- c(onsen_pars, midsen_pars, endsen_pars, tsen_pars)

# select data

data_pars <- select(pars, all_pars)

#save selected dynpars
dynpars_decorr <- names(data_pars)
saveRDS(dynpars_decorr, "T:/PhD/Data_Analysis_3/Trait_pred/Data/dynpars_decorr.rds")

#save selected SI
names(data_pars) <- gsub("onsen_|midsen_|endsen_|tsen_", "", names(data_pars))
names(data_pars) <- gsub("lin_", "", names(data_pars))
SI_decorr <- unique(names(data_pars))
saveRDS(SI_decorr, "T:/PhD/Data_Analysis_3/Trait_pred/Data/SI_decorr.rds")


names(data_pars) <- gsub("lin_SI_", "", names(data_pars))
names(data_pars) <- gsub("_r", "", names(data_pars))
names(data_pars) <- gsub("ev", "", names(data_pars))
names(data_pars) <- gsub("_ALI", "", names(data_pars))

M <- cor(data_pars, use = "pairwise.complete.obs")
corrplot(M, method = "color", tl.cex = 0.65, tl.col = "black")



#count SI
si <- unique(sapply(regmatches(names(data_pars), regexpr("_", names(data_pars)), invert = TRUE), "[[", 2))

zdf <- as.data.frame(as.table(M))
zdf

subset(zdf, abs(Freq) > 0.96)

##################################################################################################################  
