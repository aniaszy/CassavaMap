# CassavaMap v1 15 May 2020

# This code illustrates the methodology of disaggregating cassava harvested
# area data between rural population.
# Various shapefiles of different formats were used for the continental scale map.
# This example shows calculations for a chosen country (Uganda)
# and can be upscaled to the entire continent

options(stringsAsFactors=FALSE)

# Load required libraries
library(sp)
library(raster)
library(rgdal)

# Load 2014 LandScan population raster obtained from https://landscan.ornl.gov/
# This resource is not free and cannot be provided in this repository
PopRaster <- raster("Landscan_2014_pop.tif")

# Convert population raster into rural-only raster with a chosen threshold:

# Here the chosen threshold is 5000 people per pixel and above:
ThresPop <- 5000 # threshold value #

# Remove pixels with high pop density from the raster
PopRaster[PopRaster > ThresPop] <- NA

# Import country shapefile:
Shape_Adm <- readOGR("Admin_shape/UGA_shape.shp")

# Import the harvested area csv data file for the country:
Data_Adm_HA <- read.csv("AgroMaps_data/Afr-HA-Sub1.csv")

# Extract HA for each country of interest
FAOSTAT14 <- read.csv("FAO/FAOSTAT_data_2014.csv")
FAOSTAT14 <- FAOSTAT14[FAOSTAT14$Element == "Area harvested",]

# List countries that you want to work on - for example we have here Uganda only
CountryList <- "Uganda"

# Function to extract harvested area for respective country
# from the 2014 FAOSTAT dataset:
getval<-function(country,element) {
  as.numeric(FAOSTAT14$Value[FAOSTAT14$Area==country & 
                               FAOSTAT14$Element==element])
}

# Add harvested area values per country 
FAOSTAT14HA <- vector("list",length(CountryList))
for (i in 1:length(CountryList)) {
  FAOSTAT14HA <- getval(CountryList,"Area harvested")
}
FAOSTAT14HA <- data.frame(CountryList,FAOSTAT14HA)
colnames(FAOSTAT14HA) <- c("Country","FAOSTAT14HA")

# Tidy up country shapefile
Data_Adm_HA <- Data_Adm_HA[!duplicated(Data_Adm_HA), ] # Remove duplicates
colnames(Data_Adm_HA)[3] = "code" # Change column name in the data.frame to match the one in the polygon
colnames(Data_Adm_HA)[9] = "HA_AgroMaps_ADM"
colnames(Data_Adm_HA)[8] = "YEAR_AgroMaps_ADM"
Data_Adm_HA$COUNTRY_NAME <- as.character(Data_Adm_HA$COUNTRY_NAME)

## Add harvested area info to admin units
HA_SH <- merge(Shape_Adm, Data_Adm_HA, by='code')

# Visualise the production data added to the polygon
spplot(HA_SH,zcol="HA_AgroMaps_ADM")

# Function to calculate per capita harvested area values:
HarvArea_PerCapita <- function(pop_data,poly_data,fao_data){
  # Extract rural population totals per admin unit:
  AllAdmHAPopRural <- extract(pop_data, poly_data)
  AllAdmHAPopRural_ul <- sapply(AllAdmHAPopRural, 
                                 function (x) sum(unlist(x), na.rm=TRUE))
  
  # Convert to dataframe
  AllAdmHAPopRural_df <- data.frame(AllAdmHAPopRural_ul, 
                                     poly_data$HA_AgroMaps_ADM,
                                     poly_data$adm0_name,
                                     as.character(poly_data$code))
  colnames(AllAdmHAPopRural_df) <- c("PopRural","AM_HA","Country","code")
  
  # Sum HA per country
  CountryHA <- aggregate(poly_data$HA_AgroMaps_ADM, 
                         by=list(Category=poly_data$adm0_name), 
                         FUN=sum, na.rm=TRUE)
  colnames(CountryHA) <- c("Country","CountryHA_AM")
  
  
  # Add country total column
  CountryHATotal   <- merge(AllAdmHAPopRural_df,CountryHA,"Country")
  CountryHATotal14 <- merge(CountryHATotal,fao_data, by="Country")
  
  # New column with a ratio of country total in adm unit
  CountryHATotal14$Ratio <- CountryHATotal14$AM_HA / 
    CountryHATotal14$CountryHA_AM
  CountryHATotal14$HAFAO14adj <- CountryHATotal14$Ratio *
    CountryHATotal14$FAOSTAT14HA
  
  # Calculate HA per pixel
  CountryHATotal14$HARate <- CountryHATotal14$AM_HA /
    CountryHATotal14$PopRural
  CountryHATotal14$HARateFAO14Adj <- CountryHATotal14$HAFAO14adj /
    CountryHATotal14$PopRural
  
  FinalPoly <- merge(poly_data,CountryHATotal14,by="code")
  return(list(FinalPoly,CountryHATotal14))
}

Output <- HarvArea_PerCapita(PopRaster, HA_SH, FAOSTAT14HA)
FinalPolygon <- Output[[1]]
CountryHATotal14 <- Output[[2]]


# Set the maximum value of harv area hectares per pixel allowed
MPPPV <- 50 

# Calculate host density per pixel
CalcHost <- function (x_poly,x_data,x_LS){
  x_finalpoly <- merge(x_poly,x_data,by="code")
  rp <- rasterize(x_finalpoly,x_LS,'HARateFAO14Adj')
  x <- rp*x_LS
  x
}


# Reallocate harvested area if it exceeds the threshold MPPPV
ReallocateHA <- function (x_poly,x_data,x_LS,x) {
  x_Exc <- reclassify(x, matrix(c(0,MPPPV,0),
                                ncol=3, byrow=TRUE))
  x_Exc <- x_Exc-MPPPV
  x_Exc <- reclassify(x_Exc, matrix(c(-1*MPPPV-1,0.00000001,0),
                                    ncol=3, byrow=TRUE))
  y_data <- x_data
  
  x_01 <- reclassify(x, matrix(c(0,MPPPV,0,
                                 MPPPV+0.00000001,10000000,1),
                               ncol=3, byrow=TRUE))
  x_01sum <- cellStats(x_01, sum)
  y_data$HAFAO14adj[1] <- x_data$HAFAO14adj[1]-cellStats(x_01,sum)*MPPPV
  x_LSExc <- x_01*x_LS
  y_data$PopRural[1] <- x_data$PopRural[1]-cellStats(x_LSExc, sum)
  y_data$HARateFAO14Adj[1] <- y_data$HAFAO14adj[1]/y_data$PopRural[1]
  
  # Remove pixels with excess values
  x_Exc0 <- reclassify(x, matrix(c(0,MPPPV,1,
                                   MPPPV+0.00000001,10000000,0),
                                 ncol=3, byrow=TRUE))
  x_LSExc0 <- x_Exc0*x_LS # Create raster with LandScan population outside exceeded limits
  
  x_finalpoly2 <- merge(x_poly,y_data,by="code")
  rp <- rasterize(x_finalpoly2,x_LSExc0,'HARateFAO14Adj')
  resultd <- rp*x_LSExc0
  result <- list(resultd,y_data,x_LSExc0,x_01*MPPPV,x_01sum)
  result
  
}

# Supplementary function to reallocate harvested area if it exceeds the threshold MPPPV
ReallocateHA2 <- function (x_poly,x_data,x_LS,x,oldx01) {
  x_Exc <- reclassify(x, matrix(c(0,MPPPV,0),
                                ncol=3, byrow=TRUE))
  x_Exc <- x_Exc-MPPPV
  x_Exc <- reclassify(x_Exc, matrix(c(-1*MPPPV-1,0.00000001,0),
                                    ncol=3, byrow=TRUE))
  y_data <- x_data
  x_01 <- reclassify(x, matrix(c(0,MPPPV,0,
                                 MPPPV+0.00000001,10000000,1),
                               ncol=3, byrow=TRUE))
  y_data$HAFAO14adj[1] <- x_data$HAFAO14adj[1]-cellStats(x_01,sum)*MPPPV
  x_01sum <- cellStats(x_01, sum) # sum of pixels exceeded value
  x_LSExc <- x_01*x_LS # count no of people in pixels with exceeded values
  y_data$PopRural[1] <- x_data$PopRural[1]-cellStats(x_LSExc, sum)
  y_data$HARateFAO14Adj[1] <- y_data$HAFAO14adj[1]/y_data$PopRural[1]
  
  # Remove pixels with excess values
  x_Exc0 <- reclassify(x, matrix(c(0,MPPPV,1,
                                   MPPPV+0.00000001,10000000,0),
                                 ncol=3, byrow=TRUE))
  x_LSExc0 <- x_Exc0*x_LS # Create raster with LandScan population outside exceeded limits
  
  x_finalpoly2 <- merge(x_poly,y_data,by="code")
  rp <- rasterize(x_finalpoly2,x_LSExc0,'HARateFAO14Adj')
  resultd <- rp*x_LSExc0
  result <- list(resultd,y_data,x_LSExc0,oldx01+x_01*MPPPV,x_01sum)
  result
  
}

# Function redistributing HA until result has no values exceeding MPPPV
HostCtd <- function (a,b,c,d,e){
  success <- FALSE
  atmp = a
  btmp = b
  ctmp = c
  dtmp = d
  etmp = e
  while (!success){
    resulttmp <- ReallocateHA2(atmp,btmp,ctmp,dtmp,etmp)
    atmp = x_poly
    btmp = resulttmp[[2]]
    ctmp = resulttmp[[3]]
    dtmp = resulttmp[[1]]
    etmp = resulttmp[[4]]
    success <- cellStats(resulttmp[[1]], max) <= MPPPV
    print(resulttmp[[5]])
  }
  return(resulttmp)
}


### Run host calculation and redistribution:
codes <- HA_SH$code
resultFinal <- list() 
cnt = 0
for (i in 1:length(HA_SH$code)){
  cnt = cnt +1
  print(paste0("Iteration: ",cnt," in progress"))
  x_poly <- HA_SH[HA_SH$code == codes[i],]
  x_data <- CountryHATotal14[CountryHATotal14$code == codes[i],]
  x_LS <- crop(PopRaster,extent(x_poly))
  x_LS <- mask(x_LS, x_poly)
  x <- CalcHost(x_poly,x_data,x_LS)
  x
  
  x_LS01c <- reclassify(x_LS, matrix(c(0.000000001,8000000,1),
                                     ncol=3, byrow=TRUE))
  if(is.na(x_data$HARateFAO14Adj[1])){
    resultFinal[[i]] <- x
  }else{
    if(is.na(x_data$HAFAO14adj)){
      resultFinal[[i]] <- reclassify(x,matrix(c(MPPPV+0.000001,100000000,MPPPV),
                                              ncol=3,byrow=TRUE))
    }else{
      if(cellStats(x_LS01c,sum)*MPPPV<x_data$HAFAO14adj){
        resultFinal[[i]] <- reclassify(x, matrix(c(0.000000001,8000000,MPPPV),
                                                 ncol=3, byrow=TRUE))
      }else{
        if(cellStats(x, max) <= MPPPV){
          resultFinal[[i]] <- x
        }else{
          result <- ReallocateHA(x_poly,x_data,x_LS,x)
          if(cellStats(result[[1]],max) <= MPPPV){
            resultFinal[[i]] <- result[[1]]+result[[4]]
          }else{
            HostFinali <- HostCtd(x_poly,result[[2]],result[[3]],result[[1]],result[[4]])
            resultFinal[[i]] <- HostFinali[[4]]+HostFinali[[1]]  # Final host distribution result
          }
        }
      }
    }
    
  }
}

# Create the final raster with harvested area:
HA_raster <- do.call(merge,resultFinal)
plot(HA_raster)

