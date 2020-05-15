# CassavaMap v1 15 May 2020

# This code illustrates the methodology of disaggregating cassava production
# data between rural population.
# Various sProdpefiles of different formats were used for the continental scale map.
# This example shows calculations for a chosen country (Uganda)
# and can be upscaled to the entire continent

options(stringsAsFactors=FALSE)

# Load required libraries
library(sp)
library(raster)
library(rgdal)

# Load 2014 Population raster obtained from https://landscan.ornl.gov/
PopRaster <- raster("C:/Repos/CassavaDensity/Population/LandScanRural5000maxThreshold.tif")

# Convert population raster into rural-only raster with a chosen threshold:

# Here the chosen threshold is 5000 people per pixel and above:
ThresPop <- 5000 # threshold value #

# Remove pixels with high pop density from the raster
PopRaster[PopRaster > ThresPop] <- NA

# Import country shapefile:
Shape_Adm <- readOGR("Admin_shape/UGA_shape.shp")

# Import the production csv data file:
Data_Adm_Prod <- read.csv("AgroMaps_data/Afr-Prod-Sub1.csv")

# Extract production for each country of interest
FAOSTAT14 <- read.csv("FAO/FAOSTAT_data_2014.csv")
FAOSTAT14 <- FAOSTAT14[FAOSTAT14$Element == "Production",]

# List countries that you want to work on - as example Uganda is given but can be a longer list
CountryList <- "Uganda"

# Function to extract production for respective countries
# from the 2014 FAOSTAT dataset:
getval<-function(country,element) {
  as.numeric(FAOSTAT14$Value[FAOSTAT14$Area==country & 
                               FAOSTAT14$Element==element])
}

# Add Prod area values per country 
FAOSTAT14Prod <- vector("list",length(CountryList))
for (i in 1:length(CountryList)) {
  FAOSTAT14Prod <- getval(CountryList,"Production")
}
FAOSTAT14Prod <- data.frame(CountryList,FAOSTAT14Prod)
colnames(FAOSTAT14Prod) <- c("Country","FAOSTAT14Prod")

# Tidy up country sProdpefile
Data_Adm_Prod <- Data_Adm_Prod[!duplicated(Data_Adm_Prod), ] # Remove duplicates
colnames(Data_Adm_Prod)[5] = "code" # CProdnge column name in the data.frame to match the one in the polygon
colnames(Data_Adm_Prod)[11] = "Prod_AgroMaps_ADM"
colnames(Data_Adm_Prod)[10] = "YEAR_AgroMaps_ADM"
Data_Adm_Prod$COUNTRY_NAME <- as.character(Data_Adm_Prod$COUNTRY_NAME)

## Add Prod area info to admin units
Prod_SH <- merge(Shape_Adm, Data_Adm_Prod, by='code')

# Visualise the production data added to the polygon
spplot(Prod_SH,zcol="Prod_AgroMaps_ADM")

# Function to calculate per capita Prod area values:
Prod_PerCapita <- function(pop_data,poly_data,fao_data){
  # Extract rural population totals per admin unit:
  AllAdmProdPopRural <- extract(pop_data, poly_data)
  AllAdmProdPopRural_ul <- sapply(AllAdmProdPopRural, 
                                 function (x) sum(unlist(x), na.rm=TRUE))
  
  # Convert to dataframe
  AllAdmProdPopRural_df <- data.frame(AllAdmProdPopRural_ul, 
                                     poly_data$Prod_AgroMaps_ADM,
                                     poly_data$adm0_name,
                                     as.character(poly_data$code))
  colnames(AllAdmProdPopRural_df) <- c("PopRural","AM_Prod","Country","code")
  
  # Sum Prod per country
  CountryProd <- aggregate(poly_data$Prod_AgroMaps_ADM, 
                         by=list(Category=poly_data$adm0_name), 
                         FUN=sum, na.rm=TRUE)
  colnames(CountryProd) <- c("Country","CountryProd_AM")
  
  
  # Add country total column
  CountryProdTotal   <- merge(AllAdmProdPopRural_df,CountryProd,"Country")
  CountryProdTotal14 <- merge(CountryProdTotal,fao_data, by="Country")
  
  # New column with a ratio of country total in adm unit
  CountryProdTotal14$Ratio <- CountryProdTotal14$AM_Prod / 
    CountryProdTotal14$CountryProd_AM
  CountryProdTotal14$ProdFAO14adj <- CountryProdTotal14$Ratio *
    CountryProdTotal14$FAOSTAT14Prod
  
  # Calculate Prod per pixel
  CountryProdTotal14$ProdRate <- CountryProdTotal14$AM_Prod /
    CountryProdTotal14$PopRural
  CountryProdTotal14$ProdRateFAO14Adj <- CountryProdTotal14$ProdFAO14adj /
    CountryProdTotal14$PopRural
  
  FinalPoly <- merge(poly_data,CountryProdTotal14,by="code")
  return(list(FinalPoly,CountryProdTotal14))
}

Output <- Prod_PerCapita(PopRaster, Prod_SH, FAOSTAT14Prod)
FinalPolygon <- Output[[1]]
CountryProdTotal14 <- Output[[2]]


# Set the maximum value of cassava production per pixel allowed
MPPPV <- 1000

# Calculate host density per pixel
CalcHost <- function (x_poly,x_data,x_LS){
  x_finalpoly <- merge(x_poly,x_data,by="code")
  rp <- rasterize(x_finalpoly,x_LS,'ProdRateFAO14Adj')
  x <- rp*x_LS
  x
}


# Reallocate Prod area if it exceeds the threshold MPPPV
ReallocateProd <- function (x_poly,x_data,x_LS,x) {
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
  y_data$ProdFAO14adj[1] <- x_data$ProdFAO14adj[1]-cellStats(x_01,sum)*MPPPV
  x_LSExc <- x_01*x_LS
  y_data$PopRural[1] <- x_data$PopRural[1]-cellStats(x_LSExc, sum)
  y_data$ProdRateFAO14Adj[1] <- y_data$ProdFAO14adj[1]/y_data$PopRural[1]
  
  # Remove pixels with excess values
  x_Exc0 <- reclassify(x, matrix(c(0,MPPPV,1,
                                   MPPPV+0.00000001,10000000,0),
                                 ncol=3, byrow=TRUE))
  x_LSExc0 <- x_Exc0*x_LS # Create raster with LandScan population outside exceeded limits
  
  x_finalpoly2 <- merge(x_poly,y_data,by="code")
  rp <- rasterize(x_finalpoly2,x_LSExc0,'ProdRateFAO14Adj')
  resultd <- rp*x_LSExc0
  result <- list(resultd,y_data,x_LSExc0,x_01*MPPPV,x_01sum)
  result
  
}

# Supplementary function to reallocate Prod area if it exceeds the threshold MPPPV
ReallocateProd2 <- function (x_poly,x_data,x_LS,x,oldx01) {
  x_Exc <- reclassify(x, matrix(c(0,MPPPV,0),
                                ncol=3, byrow=TRUE))
  x_Exc <- x_Exc-MPPPV
  x_Exc <- reclassify(x_Exc, matrix(c(-1*MPPPV-1,0.00000001,0),
                                    ncol=3, byrow=TRUE))
  y_data <- x_data
  x_01 <- reclassify(x, matrix(c(0,MPPPV,0,
                                 MPPPV+0.00000001,10000000,1),
                               ncol=3, byrow=TRUE))
  y_data$ProdFAO14adj[1] <- x_data$ProdFAO14adj[1]-cellStats(x_01,sum)*MPPPV
  x_01sum <- cellStats(x_01, sum) # sum of pixels with exceeded value
  x_LSExc <- x_01*x_LS # count population in pixels with exceeded values
  y_data$PopRural[1] <- x_data$PopRural[1]-cellStats(x_LSExc, sum)
  y_data$ProdRateFAO14Adj[1] <- y_data$ProdFAO14adj[1]/y_data$PopRural[1]
  
  # Remove pixels with excess values
  x_Exc0 <- reclassify(x, matrix(c(0,MPPPV,1,
                                   MPPPV+0.00000001,10000000,0),
                                 ncol=3, byrow=TRUE))
  x_LSExc0 <- x_Exc0*x_LS # Create raster with LandScan population outside exceeded limits
  
  x_finalpoly2 <- merge(x_poly,y_data,by="code")
  rp <- rasterize(x_finalpoly2,x_LSExc0,'ProdRateFAO14Adj')
  resultd <- rp*x_LSExc0
  result <- list(resultd,y_data,x_LSExc0,oldx01+x_01*MPPPV,x_01sum)
  result
  
}

# Function redistributing Prod until result has no values exceeding MPPPV
HostCtd <- function (a,b,c,d,e){
  success <- FALSE
  atmp = a
  btmp = b
  ctmp = c
  dtmp = d
  etmp = e
  while (!success){
    resulttmp <- ReallocateProd2(atmp,btmp,ctmp,dtmp,etmp)
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
codes <- Prod_SH$code
resultFinal <- list() 
cnt = 0
for (i in 1:length(Prod_SH$code)){
  cnt = cnt +1
  print(paste0("Iteration: ",cnt," in progress"))
  x_poly <- Prod_SH[Prod_SH$code == codes[i],]
  x_data <- CountryProdTotal14[CountryProdTotal14$code == codes[i],]
  x_LS <- crop(PopRaster,extent(x_poly))
  x_LS <- mask(x_LS, x_poly)
  x <- CalcHost(x_poly,x_data,x_LS)
  x
  
  x_LS01c <- reclassify(x_LS, matrix(c(0.000000001,8000000,1),
                                     ncol=3, byrow=TRUE))
  if(is.na(x_data$ProdRateFAO14Adj[1])){
    resultFinal[[i]] <- x
  }else{
    if(is.na(x_data$ProdFAO14adj)){
      resultFinal[[i]] <- reclassify(x,matrix(c(MPPPV+0.000001,100000000,MPPPV),
                                              ncol=3,byrow=TRUE))
    }else{
      if(cellStats(x_LS01c,sum)*MPPPV<x_data$ProdFAO14adj){
        resultFinal[[i]] <- reclassify(x, matrix(c(0.000000001,8000000,MPPPV),
                                                 ncol=3, byrow=TRUE))
      }else{
        if(cellStats(x, max) <= MPPPV){
          resultFinal[[i]] <- x
        }else{
          result <- ReallocateProd(x_poly,x_data,x_LS,x)
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

# Create the final raster with Prod area:
Prod_raster <- do.call(merge,resultFinal)
plot(Prod_raster)

