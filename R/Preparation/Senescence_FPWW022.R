
#########################################
# Transform senescence ratings to GDDAH 
#########################################

library(readxl)
library(writexl)
library(tidyr)
library(dplyr)

#read scenescence scorings and heading data
data <- read_excel("O:/Projects/KP0011/1/Data preparation/Data/raw/Senescence_FPWW022.xlsx")
data_heading <- read_excel("O:/Projects/KP0011/1/Data preparation/Data/out/Heading_FPWW022_DAS.xlsx")

###########################################################################################################################

### first, the heading information needs to be completed for the lot which was not scored

  # define mean heading date for checks to attribute to check plots of the lot which was not scored
    d <- data_heading[c(1:378),c(5:8)]
    checks <- c("SURETTA", "CH CLARO", "CH NARA") #define which cultivars represent checks
    d <- d[!d$Gen_Name %in% checks,] #excude the checks from the df
  # set mean headingDAS and mean heading GDDAS based on the scored lot
    d <- add_row(d, Gen_Name = "SURETTA", HeadingDAS = 219, heading_date = as.Date("2018-05-24", "%Y-%m-%d"), heading_GDDAS = 1279.45)
    d <- add_row(d, Gen_Name = "CH NARA", HeadingDAS = 219, heading_date = as.Date("2018-05-24", "%Y-%m-%d"), heading_GDDAS = 1279.45)
    d <- add_row(d, Gen_Name = "CH CLARO", HeadingDAS = 218, heading_date = as.Date("2018-05-23", "%Y-%m-%d"), heading_GDDAS = 1260.91)

  # add this heading information to unscored lot
    data2 <- right_join(d, data, by = "Gen_Name")
   
  # HEADING DATA READY

###########################################################################################################################

### rearrange senescence scoring data
    
  # seperate Fl0 and Cnp rating data
    data_Fl0 <- select(data2, c(1:8, starts_with("Fl0")))
    data_Cnp <- select(data2, c(1:8, starts_with("Can")))
  
  # transform both to long format
    data_long_Fl0 <- data_Fl0 %>% gather(object_rating, SnsFl0, `Fl0 14.06.2018`:`Fl0 11.07.2018`, factor_key = TRUE)
    data_long_Cnp <- data_Cnp %>% gather(object_rating, SnsCnp, `Can 14.06.2018`:`Can 11.07.2018`, factor_key = TRUE)
  
  # bind together
    data_long <- cbind(data_long_Fl0, data_long_Cnp["SnsCnp"])
  
  # get scoring dates
    data_long$grading_date <- sapply(lapply(lapply(as.character(data_long$object_rating), strsplit, " ", fixed = TRUE), unlist), function(x) x=as.character(x[2]))
    data_long$grading_date <- as.Date(data_long$grading_date, "%d.%m.%Y")
  
  # calculate graingDAS, heading DAS and gradingDAH and add to dataframe
    data_long$grading_DAS <- data_long$grading_date - as.Date("17.10.2017", "%d.%m.%Y")
    data_long$heading_DAS <- as.difftime(data_long$HeadingDAS, units = "days")
    data_long$grading_DAH <- data_long$grading_DAS - data_long$heading_DAS
  
  # SENESCENCE SCORING DATA READY

###########################################################################################################################  
  
### Get GDD data

  # load temperature data (from agrometeo)
    data_climate <- read.csv("O:/Projects/KP0011/1/Data preparation/Data/raw/data_climate_2018.csv", sep = ";", na.strings = "?")
    data_climate$time_string <- strptime(data_climate$time_string, "%d.%m.%Y")
    
  # calcualte mean daily temperature, using hourly means 
    day_temp <- aggregate(data_climate$mean_temp, list(day = cut(data_climate$time_string, "days")), mean, na.rm = TRUE)
    names(day_temp)[2] <- "mean"
  
  # get min and max hourly temperatures and calculate meanMM temperature  
    day_temp$min <- aggregate(data_climate$mean_temp, list(day = cut(data_climate$time_string, "days")), min, na.rm = TRUE)$x
    day_temp$max <- aggregate(data_climate$mean_temp, list(day = cut(data_climate$time_string, "days")), max, na.rm = TRUE)$x
    day_temp$meanMM <- (day_temp$max + day_temp$min)/2
    
    day_temp$day <- as.Date(day_temp$day)
    
  # set mean values < 0 to 0 (base temperature)
    day_temp$mean_calc <- ifelse(day_temp$mean > 0, day_temp$mean, 0)
    day_temp$meanMM_calc <- ifelse(day_temp$meanMM > 0, day_temp$meanMM, 0)
  
  # calcualte GDD and GDDMM
    day_temp$GDD <- cumsum(day_temp$mean_calc)
    day_temp$GDDMM <- cumsum(day_temp$meanMM_calc)

  # save GDD data
    write_xlsx(day_temp, "O:/Projects/KP0011/1/Data preparation/Data/raw/GDDAH_2018.xlsx")
    
  # GDD DATA READY
    
###########################################################################################################################  
  
### Calculate GDDAH for each grading date and add to dataframe 
  
  # Calcualte gradingGDDAH as gradingGDD - headingGDD for each plot and each scoring date
  # This takes a while
    grading_GDDAH <- NULL
    for (i in 1:nrow(data_long)) {
        
        if (!is.na(data_long$heading_date)[i]) {
          
          grading_GDDAH[i] <- day_temp[day_temp$day == paste(data_long$grading_date)[i], 9] -
            day_temp[day_temp$day == paste(data_long$heading_date)[i], 9]
          
        }
        
        else paste("NA")
        
        print(grading_GDDAH[i])
  
      }
    data_long$grading_GDDAH <- grading_GDDAH 

  # calcualte grading_  GDDAS
    data_long$grading_GDDAS <- as.numeric(data_long$heading_GDDAS) + as.numeric(data_long$grading_GDDAH)

  # tidy up data
    data_long <- data_long %>% select(Plot_ID, Lot, RangeLot, RowLot, Gen_Name, grading_date, 
                                     heading_date, heading_DAS, heading_GDDAS, grading_DAS, 
                                     grading_DAH, grading_GDDAH, grading_GDDAS, SnsFl0, SnsCnp)
  
  # save 
    write_xlsx(as.data.frame(data_long), "O:/Projects/KP0011/1/Data preparation/Data/out/Senescence_FPWW022_DAS.xlsx")
