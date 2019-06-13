
#########################################
# Estimate heading date in DAS
#########################################

library("xlsx")

d.head <- read.csv("O:/Projects/KP0011/1/Data preparation/Data/raw/Heading_FPWW018.csv", na.string = "NA", sep = ";")

plots <- d.head[,1]

HeadingDAS<-NULL

for (i in 1:length(plots)){
  x<- d.head[i,-c(1:6)]
  
  if(is.na(x[3])){
    print("NA")
    HeadingDAS<-c(HeadingDAS, NA)
  } else if(x[1] == 60){
    print("205")
    HeadingDAS<-c(HeadingDAS, 205)
  } else if (x[2] == 60){
    print("211")
    HeadingDAS<-c(HeadingDAS, 211)
  } else if (x[3] == 60){
    print("214")
    HeadingDAS<-c(HeadingDAS, 214)    
  } else if (x[1] > 56 & x[1] < 60 & x[2] >= 60){
    print("206")
    HeadingDAS<-c(HeadingDAS, 206)
  } else if (x[1] > 53 & x[1] < 60 & x[2] >= 60){
    print("207")
    HeadingDAS<-c(HeadingDAS, 207)
  } else if (x[1] > 50 & x[1] < 60 & x[2] >= 60){
    print("208")
    HeadingDAS<-c(HeadingDAS, 208)
  } else if (x[2] > 56 & x[2] < 60 & x[3] >= 60){
    print("212")
    HeadingDAS<-c(HeadingDAS, 212)
  } else if (x[2] > 50 & x[2] < 60 & x[3] >= 60){
    print("213")
    HeadingDAS<-c(HeadingDAS, 213)
  } else if (x[3] > 56 & x[3] < 60){
    print("215")
    HeadingDAS<-c(HeadingDAS, 215)
  } else if (x[3] > 53 & x[3] < 60){
    print("216")
    HeadingDAS<-c(HeadingDAS, 216)
  } else if (x[3] > 50 & x[3] < 60){
    print("217")
    HeadingDAS<-c(HeadingDAS, 217)
  } else if (x[1] > 45 & x[2] > 60){
    print("209")
    HeadingDAS<-c(HeadingDAS, 209)
  } else if (x[1] > 39 & x[2] > 60){
    print("210")
    HeadingDAS<-c(HeadingDAS, 210) 
  } else if (x[2] > 39 & x[3] > 60){
    print("213")
    HeadingDAS<-c(HeadingDAS, 213)
  } else{
    print("NA")
    HeadingDAS<-c(HeadingDAS, NA)
  }
  
}

HeadingDAS <- cbind(d.head[,1:5], HeadingDAS)
HeadingDAS

hist(HeadingDAS$HeadingDAS)

HeadingDAS$heading_date <- as.Date("01.11.2016", "%d.%m.%Y") + HeadingDAS$HeadingDAS

#Quality Check

HeadingDAS[HeadingDAS$Gen_Name=="CH CLARO", ] #OK
HeadingDAS[HeadingDAS$Gen_Name=="SURETTA", ] #OK
HeadingDAS[HeadingDAS$Gen_Name=="NARA",] #OK


#GDDAS

heading_GDDAS <- NULL

day_temp <- read.csv("O:/Projects/KP0011/1/Data preparation/Data/raw/GDDAH_2017.csv")

for (i in 1:nrow(HeadingDAS)) {
  
  if (!is.na(HeadingDAS$heading_date)[i]) {
    
    heading_GDDAS[i] <- day_temp[day_temp$day == paste(HeadingDAS$heading_date)[i], 6] -
      day_temp[day_temp$day == paste(as.Date("01.11.2016", "%d.%m.%Y")), 6]
    
  }
  
  else {heading_GDDAS[i] <- NA}
  
  print(heading_GDDAS[i])
  
}

HeadingDAS$heading_GDDAS <- round(heading_GDDAS, 2)

d <- HeadingDAS[c(1:378),c(5:8)]
checks <- c("SURETTA", "CH CLARO", "NARA")
d <- d[!d$Gen_Name %in% checks,]

d <- add_row(d, Gen_Name = "SURETTA", HeadingDAS = 209, heading_date = as.Date("2017-05-29", "%Y-%m-%d"), heading_GDDAS = 1100.12)
d <- add_row(d, Gen_Name = "NARA", HeadingDAS = 211, heading_date = as.Date("2017-05-31", "%Y-%m-%d"), heading_GDDAS = 1143.56)
d <- add_row(d, Gen_Name = "CH CLARO", HeadingDAS = 211, heading_date = as.Date("2017-05-31", "%Y-%m-%d"), heading_GDDAS = 1143.56)

##add values to Lot 6 data

data <- right_join(d, HeadingDAS[1:5])

write.csv(data, "O:/Projects/KP0011/1/Data preparation/Data/out/Heading_FPWW018_DAS.csv", row.names = FALSE)
