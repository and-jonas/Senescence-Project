
##################################################
# Transform senescence ratings to GDDAH and GDDAS 
##################################################

library(readxl)
library(tidyr)
library(dplyr)
library(xlsx)

#read data

data2016 <- read.csv("O:/Projects/KP0011/1/Data preparation/Data/raw/Senescence_ww012.csv", sep = ",", check.names = FALSE)
data_heading2016 <- read.csv("O:/Projects/KP0011/1/Data preparation/Data/out/Heading_FPWW012_DAS.csv")

data2016$Plot_ID == data_heading2016$Plot_ID

#merge data 

data2 <- full_join(data_heading2016, data2016)

#data tidying
##Fl0 and Cnp in separate columns

data3 <- as_data_frame(data2)

data_Fl0 <- select(data3, c(1:8, starts_with("Fl0")))
data_Cnp <- select(data3, c(1:8, starts_with("Cnp")))

data_long_Fl0 <- data_Fl0 %>% gather(object_rating, SnsFl0, `Fl0 28.06.2016`:`Fl0 20.07.2016`, factor_key = TRUE)
data_long_Cnp <- data_Cnp %>% gather(object_rating, SnsCnp, `Cnp 28.06.2016`:`Cnp 20.07.2016`, factor_key = TRUE)

data_long <- bind_cols(data_long_Fl0, data_long_Cnp)

data_long$grading_date <- sapply(lapply(lapply(as.character(data_long$object_rating), strsplit, " ", fixed = TRUE), unlist), function(x) x=as.character(x[2]))
data_long$grading_date <- as.Date(data_long$grading_date, "%d.%m.%Y")

data_long$grading_DAS <- data_long$grading_date - as.Date("13.10.2015", "%d.%m.%Y")
data_long$heading_DAS <- as.difftime(data_long$HeadingDAS, units = "days")
data_long$grading_DAH <- data_long$grading_DAS - data_long$heading_DAS

#calculate GDD

data_climate <- read.csv("O:/Projects/KP0011/1/Data preparation/Data/raw/data_climate_2016.csv", sep = ";")
data_climate$time_string <- strptime(data_climate$time_string, "%d.%m.%Y %H:%M")

day_temp <- aggregate(data_climate$mean_temp, list(day = cut(data_climate$time_string, "days")), mean)
names(day_temp)[2] <- "mean"

day_temp$min <- aggregate(data_climate$mean_temp, list(day = cut(data_climate$time_string, "days")), min)$x
day_temp$max <- aggregate(data_climate$mean_temp, list(day = cut(data_climate$time_string, "days")), max)$x
day_temp$meanMM <- (day_temp$max + day_temp$min)/2

day_temp$day <- as.Date(day_temp$day)

day_temp$GDD <- cumsum(day_temp$mean)
day_temp$GDDMM <- cumsum(day_temp$meanMM)

# save GDD data
write.csv(day_temp[,-c(5,7)], "O:/Projects/KP0011/1/Data preparation/Data/raw/GDDAH_2016.csv")


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
names(data_long)[2] <- paste("Lot")

data_long$grading_GDDAS <- data_long$heading_GDDAS + data_long$grading_GDDAH

#Tidy up data

data_long <- data_long %>% select(Plot_ID, Lot, RangeLot, RowLot, Gen_Name, grading_date, 
                                  heading_date, heading_DAS, heading_GDDAS, grading_DAS, grading_DAH, grading_GDDAH, grading_GDDAS, SnsFl0, SnsCnp)


write.csv(as.data.frame(data_long), "O:/Projects/KP0011/1/Data preparation/Data/out/Senescence_FPWW012_DAS.csv", row.names = FALSE)

