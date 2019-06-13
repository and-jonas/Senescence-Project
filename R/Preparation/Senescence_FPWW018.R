
#########################################
# Transform senescence ratings to GDDAH 
#########################################

library(readxl)
library(tidyr)
library(dplyr)
library(xlsx)

#read data

data <- read_excel("O:/Projects/KP0011/1/Data preparation/Data/raw/Senescence_FPWW018.xlsx", 1)

data_heading <- read.csv("O:/Projects/KP0011/1/Data preparation/Data/out/Heading_FPWW018_DAS.csv")
 
#use heading dates of lot 2 on lot 6

##drop check values and replace by mean values 

  d <- data_heading[c(1:378),-c(5:8)]
  checks <- c("SURETTA", "CH CLARO", "NARA")
  d <- d[!d$Gen_Name %in% checks,]
  
  d <- add_row(d, Gen_Name = "SURETTA", HeadingDAS = 209, heading_date = as.factor(as.Date("2017-05-29", "%Y-%m-%d")), heading_GDDAS = 1100.12)
  d <- add_row(d, Gen_Name = "NARA", HeadingDAS = 211, heading_date = as.factor(as.Date("2017-05-31", "%Y-%m-%d")), heading_GDDAS = 1143.56)
  d <- add_row(d, Gen_Name = "CH CLARO", HeadingDAS = 211, heading_date = as.factor(as.Date("2017-05-31", "%Y-%m-%d")), heading_GDDAS = 1143.56)
   
##add values to Lot 6 data
  
  data2 <- right_join(d, data, by = "Gen_Name")

#data tidying
##Fl0 and Cnp in separate columns
    
 data3 <- as_data_frame(data2)
 
  data_Fl0 <- select(data3, c(1:8, starts_with("Fl0")))
  data_Cnp <- select(data3, c(1:8, starts_with("Can")))
  
  data_long_Fl0 <- data_Fl0 %>% gather(object_rating, SnsFl0, `Fl0 18.06.2017`:`Fl0 16.07.2017`, factor_key = TRUE)
  data_long_Cnp <- data_Cnp %>% gather(object_rating, SnsCnp, `Can 18.06.2017`:`Can 16.07.2017`, factor_key = TRUE)
  
  data_long <- bind_cols(data_long_Fl0, data_long_Cnp)
  
  data_long$grading_date <- sapply(lapply(lapply(as.character(data_long$object_rating), strsplit, " ", fixed = TRUE), unlist), function(x) x=as.character(x[2]))
    data_long$grading_date <- as.Date(data_long$grading_date, "%d.%m.%Y")
  
  data_long$grading_DAS <- data_long$grading_date - as.Date("01.11.2016", "%d.%m.%Y")
    data_long$heading_DAS <- as.difftime(data_long$HeadingDAS, units = "days")
  data_long$grading_DAH <- data_long$grading_DAS - data_long$heading_DAS


#calculate GDD
  
  data_climate <- read.csv("O:/Projects/KP0011/1/Data preparation/Data/raw/data_climate_2017.csv", sep = ",")
  data_climate$time_string <- strptime(data_climate$time_string, "%d.%m.%Y %H:%M")
  
  day_temp <- aggregate(data_climate$mean_temp, list(day = cut(data_climate$time_string, "days")), mean)
  names(day_temp)[2] <- "mean"
  
  day_temp$min <- aggregate(data_climate$mean_temp, list(day = cut(data_climate$time_string, "days")), min)$x
  day_temp$max <- aggregate(data_climate$mean_temp, list(day = cut(data_climate$time_string, "days")), max)$x
  day_temp$meanMM <- (day_temp$max + day_temp$min)/2
  
  day_temp$day <- as.Date(day_temp$day)
  
  day_temp$GDD <- cumsum(day_temp$mean)
  day_temp$GDDMM <- cumsum(day_temp$meanMM)

#save GDD data
  
  write.csv(day_temp, "O:/Projects/KP0011/1/Data preparation/Data/raw/GDDAH_2017.csv")
  
#calculate GDDAH and add to data_long
  
  grading_GDDAH <- NULL
  
    for (i in 1:nrow(data_long)) {
      
      if (!is.na(data_long$heading_date)[i]) {
        
        grading_GDDAH[i] <- day_temp[day_temp$day == paste(data_long$grading_date)[i], 7] -
          day_temp[day_temp$day == paste(data_long$heading_date)[i], 7]
        
      }
      
      else paste("NA")
      
      print(grading_GDDAH[i])

    }
  
  data_long$grading_GDDAH <- grading_GDDAH 
  data_long <- as.data.frame(data_long)
  
  data_long$grading_GDDAS <- as.numeric(data_long$heading_GDDAS) + as.numeric(data_long$grading_GDDAH)

#Tidy up data

 data_long <- data_long %>% select(Plot_ID, Lot, RangeLot, RowLot, Gen_Name, grading_date, 
                                   heading_date, heading_DAS, heading_GDDAS, grading_DAS, grading_DAH, grading_GDDAH, grading_GDDAS, SnsFl0, SnsCnp)
  
  
write.csv(as.data.frame(data_long), "O:/Projects/KP0011/1/Data preparation/Data/out/Senescence_FPWW018_DAS.csv", row.names = FALSE)
