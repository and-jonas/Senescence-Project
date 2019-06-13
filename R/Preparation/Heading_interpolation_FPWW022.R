
# Function to change certain variables in a dataset from numeric into factors
FunNumerize <- function(DataSetName, VariableNames){
  VariableNumbers <- match(VariableNames, colnames(DataSetName))
  for(CGIndex in VariableNumbers){
    DataSetName[[CGIndex]] <- as.numeric(as.character(DataSetName[[CGIndex]]))
  }
  return(DataSetName)
}

#########################################
# Estimate heading date
#########################################

## Get headingDAS and heading date by interpolation of ratings  
  
  #load rating data
  d.head <- read.csv("O:/Projects/KP0011/1/Data preparation/Data/raw/Heading_FPWW022.csv", na.strings = c("NA",""), sep = ";")
  #numerize ratings
  d.head <- FunNumerize(d.head, names(d.head)[6:13])
  
  #interpolate ratings
  plots <- d.head[,1]
  HeadingDAS<-NULL
  for (i in 1:length(plots)){
    x<- d.head[i,-c(1:6)]
    
    if(x[1] == 60){
      print("213")
      HeadingDAS<-c(HeadingDAS, 213)
    } else if (x[2] == 60){
      print("216")
      HeadingDAS<-c(HeadingDAS, 216)
    } else if (x[3] == 60){
      print("218")
      HeadingDAS<-c(HeadingDAS, 218)    
    } else if (x[4] == 60){
      print("220")
      HeadingDAS<-c(HeadingDAS, 220)    
    } else if (x[5] == 60){
      print("223")
      HeadingDAS<-c(HeadingDAS, 223)    
    } else if (x[6] == 60){
      print("225")
      HeadingDAS<-c(HeadingDAS, 225)    
    } else if (x[1] > 56 & x[1] < 60 & x[2] >= 60){
      print("214")
      HeadingDAS<-c(HeadingDAS, 214)
    } else if (x[1] > 53 & x[1] < 60 & x[2] >= 60){
      print("215")
      HeadingDAS<-c(HeadingDAS, 215)
    } else if (x[2] > 56 & x[2] < 60 & x[3] >= 60){
      print("217")
      HeadingDAS<-c(HeadingDAS, 217)
    } else if (x[2] > 50 & x[2] < 60 & x[3] >= 60){
      print("218")
      HeadingDAS<-c(HeadingDAS, 218)
    } else if (x[3] > 56 & x[3] < 60 & x[4] >= 60){
      print("219")
      HeadingDAS<-c(HeadingDAS, 219)
    } else if (x[3] > 50 & x[3] < 60 & x[4] >= 60){
      print("220")
      HeadingDAS<-c(HeadingDAS, 220)
    } else if (x[4] > 56 & x[4] < 60 & x[5] >= 60){
      print("221")
      HeadingDAS<-c(HeadingDAS, 221)
    } else if (x[4] > 53 & x[4] < 60 & x[5] >= 60){
      print("222")
      HeadingDAS<-c(HeadingDAS, 222)
    } else if (x[4] > 50 & x[4] < 60 & x[5] >= 60){
      print("223")
      HeadingDAS<-c(HeadingDAS, 223)
    } else if (x[5] > 56 & x[5] < 60 & x[6] >= 60){
      print("224")
      HeadingDAS<-c(HeadingDAS, 224) 
    } else if (x[5] > 50 & x[5] < 60 & x[6] >= 60){
      print("225")
      HeadingDAS<-c(HeadingDAS, 225) 
    } else if (x[6] > 56 & x[6] < 60 & x[7] >= 60){
      print("226")
      HeadingDAS<-c(HeadingDAS, 226) 
    } else{
      print("NA")
      HeadingDAS<-c(HeadingDAS, NA)
    }
    
}
  
  #add Plot information to interpolated ratings
  HeadingDAS <- cbind(d.head[,1:5], HeadingDAS)
  
  #check heading date distribution
  hist(HeadingDAS$HeadingDAS)
  
  #add heading date by adding headingDAS to sowing date
  HeadingDAS$heading_date <- as.Date("17.10.2017", "%d.%m.%Y") + HeadingDAS$HeadingDAS
  
  #Quality Check
  HeadingDAS[HeadingDAS$Gen_Name=="CH CLARO", ] #OK
  HeadingDAS[HeadingDAS$Gen_Name=="SURETTA", ] #OK
  HeadingDAS[HeadingDAS$Gen_Name=="CH NARA", ] #OK

##############################################################################################################################

## Get headingGDDAS 
  
  #get GDDAH data for the growing season
  day_temp <- read.csv("O:/Projects/KP0011/1/Data preparation/Data/raw/GDDAH_2018.csv", sep = ",")[-1]
  #convert to dates
  day_temp$day <- as.Date(day_temp$day, "%Y-%m-%d")
  
  #calculate headingGDDAS for each plot
  heading_GDDAS <- NULL
  for (i in 1:nrow(HeadingDAS)) {
    
    if (!is.na(HeadingDAS$heading_date)[i]) {
      
      heading_GDDAS[i] <- day_temp[day_temp$day == paste(HeadingDAS$heading_date)[i], 6] -
        day_temp[day_temp$day == paste(as.Date("17.10.2017", "%d.%m.%Y")), 6]
      
    }
    
    else {heading_GDDAS[i] <- NA}
    
    print(heading_GDDAS[i])
    
  }
  HeadingDAS$heading_GDDAS <- round(heading_GDDAS, 2)

  #define mean heading date for checks to attribute to check plots of the lot which was not scored
  d <- HeadingDAS[c(1:378),c(5:8)]
  checks <- c("SURETTA", "CH CLARO", "CH NARA") #define which cultivars represent checks
  d <- d[!d$Gen_Name %in% checks,] #excude the checks from the df
  #set mean headingDAS and mean heading GDDAS based on the scored lot
  d <- add_row(d, Gen_Name = "SURETTA", HeadingDAS = 219, heading_date = as.Date("2018-05-24", "%Y-%m-%d"), heading_GDDAS = 1279.45)
  d <- add_row(d, Gen_Name = "CH NARA", HeadingDAS = 219, heading_date = as.Date("2018-05-24", "%Y-%m-%d"), heading_GDDAS = 1279.45)
  d <- add_row(d, Gen_Name = "CH CLARO", HeadingDAS = 218, heading_date = as.Date("2018-05-23", "%Y-%m-%d"), heading_GDDAS = 1260.91)

  ##add derived mean values to the df again
  data <- right_join(d, HeadingDAS[1:5]) %>% 
    select(Plot_ID, Lot, RangeLot, RowLot, Gen_Name, heading_date, HeadingDAS, heading_GDDAS)

  write_xlsx(data, "O:/Projects/KP0011/1/Data preparation/Data/out/Heading_FPWW022_DAS.xlsx")
