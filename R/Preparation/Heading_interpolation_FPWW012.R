
#########################################
# Estimate heading date in DAS and GDDAS
#########################################

#DAS

d.head <- read.csv("O:/Projects/KP0011/1/Data preparation/Data/raw/ww012_Heading_raw.csv")

plots <- d.head[,1]

HeadingDAS<-NULL

for (i in 1:length(plots)){
  x<- d.head[i,-c(1:5)]
  
  if (x[1] == 60){
    print("226")
    HeadingDAS<-c(HeadingDAS, 226)
  } else if (x[1] > 57){
    print("227")
    HeadingDAS<-c(HeadingDAS, 227)
  } else if (x[2] == 60){
    print("228")
    HeadingDAS<-c(HeadingDAS, 228)
  } else if (x[2] > 57){
    print("229")
    HeadingDAS<-c(HeadingDAS, 229)
  } else if (x[3] == 60){
    print("230")
    HeadingDAS<-c(HeadingDAS, 230)
  } else if (x[3] > 57){
    print("231")
    HeadingDAS<-c(HeadingDAS, 231)
  } else if (x[4] == 60){
    print("232")
    HeadingDAS<-c(HeadingDAS, 232)
  } else if (x[4] > 57){
    print("233")
    HeadingDAS<-c(HeadingDAS, 233)
  } else if (x[5] == 60){
    print("234")
    HeadingDAS<-c(HeadingDAS, 234)
  } else if (x[5] > 57){
    print("235")
    HeadingDAS<-c(HeadingDAS, 235)
  } else if (x[6] == 60){
    print("236")
    HeadingDAS<-c(HeadingDAS, 236)
  } else if (x[6] > 57){
    print("237")
    HeadingDAS<-c(HeadingDAS, 237)
  } else if (x[6] > 0){
    print("238")
    HeadingDAS<-c(HeadingDAS, 238)
  } else{
    print("239")
    HeadingDAS<-c(HeadingDAS, 239)
  }
  
}

HeadingDAS <- cbind(d.head[,1:5], HeadingDAS)
HeadingDAS

hist(HeadingDAS$HeadingDAS)

HeadingDAS$heading_date <- as.Date("13.10.2015", "%d.%m.%Y") + HeadingDAS$HeadingDAS

#GDDAS

heading_GDDAS <- NULL

day_temp <- read.csv("O:/Projects/KP0011/1/Data preparation/Data/raw/GDDAH_2016.csv")

for (i in 1:nrow(HeadingDAS)) {
  
  if (!is.na(HeadingDAS$heading_date)[i]) {
    
    heading_GDDAS[i] <- day_temp[day_temp$day == paste(HeadingDAS$heading_date)[i], 6] -
      day_temp[day_temp$day == paste(as.Date("13.10.2015", "%d.%m.%Y")), 6]
    
  }
  
  else paste("NA")
  
  print(heading_GDDAS[i])
  
}

HeadingDAS$heading_GDDAS <- heading_GDDAS


write.csv(HeadingDAS, "O:/Projects/KP0011/1/Data preparation/Data/out/Heading_FPWW012_DAS.csv", row.names = FALSE)
