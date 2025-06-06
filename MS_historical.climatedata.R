#Purpose: Pulling climate date from ClimateNa for each location
  #Using ClimateNA because the data goes back to 1901 while terra only goes back to 1950
#How output is used: 
  #Climate data is added to the full dataset for each catalogNumber (i.e., herbarium sheets)
#Author: LAS
#Date: 04/19/2024
#R version: 4.3.3

#load package
require(ClimateNAr) #see documentation in project folder on how I got this to install but this is stupid
require(tidyverse)
#install.packages("elevatr")
require(elevatr) #need this package to get elevation points for each lat/long
#install.packages("sf")
require(sf) #need to create projection of spatial data. 

#load data
x <- read_csv("Cleaned_data/QL_fulldf.csv")

#need to subset data by the various "historical" bins for the ClimateNA calculations. Can only go back to 1901
  #Data has to be in a very specific order
  #ID1, ID2, lat, long, elevation

#1901-1930
dat.1901_30 <- x %>% 
  select(c(catalogNumber,decimalLatitude,decimalLongitude,decade)) %>% #pulling only the location data for each herbarium sheet
  group_by(catalogNumber) %>% #grouping by herbarium sheet
  summarise_all(mean) %>%  #averaging the lat long (they don't change. I double checked)
  rename(lat = "decimalLatitude") %>% 
  rename(lon = "decimalLongitude") %>% 
  rename(ID1 = "catalogNumber") 

dat.1901_30 <- dat.1901_30 %>% 
  filter(decade %in% c(1900,1910,1920)) %>%  #filtering out only the rows that contain the data provided within the vector
  add_column(ID2 = ".", .after = "ID1", ) %>% #adding a dummy column 
  relocate(lat, .before = lon) %>% 
  select(!decade)

#getting elevation points for each lat/long point
#need projection point for data first 
sf_points <- st_as_sf(dat.1901_30, coords = c("lon","lat"), crs=4326) #4326 refers to the WGS 84 
  #geographic coordinate system, which is commonly used for GPS coordinates and other global positioning systems

eldat <-  get_elev_point(sf_points) #pulls elevation points for each lat/long in dataframe 
dat.1901_30 <- cbind(dat.1901_30, eldat) #binding the elevation data with the latlong data. Using cbind because they are in order 

dat.1901_30 <-  dat.1901_30 %>% 
  select(c("ID1", "ID2", "lat", "lon", "elevation")) %>% 
  rename(el = "elevation")

#1931-1960
dat.1931_60 <- x %>% 
  select(c(catalogNumber,decimalLatitude,decimalLongitude,decade)) %>% 
  group_by(catalogNumber) %>% 
  summarise_all(mean) %>%  
  rename(lat = "decimalLatitude") %>% 
  rename(lon = "decimalLongitude") %>% 
  rename(ID1 = "catalogNumber") 

dat.1931_60 <- dat.1931_60 %>% 
  filter(decade %in% c(1930,1940,1950)) %>% 
  add_column(ID2 = ".", .after = "ID1", ) %>% 
  relocate(lat, .before = lon) %>% 
  select(!decade)

#getting elevation points for each lat/long point
#need projection point for data first 
sf_points <- st_as_sf(dat.1931_60, coords = c("lon","lat"), crs=4326) 

eldat <-  get_elev_point(sf_points) 
dat.1931_60 <- cbind(dat.1931_60, eldat) 

dat.1931_60 <-  dat.1931_60 %>% 
  select(c("ID1", "ID2", "lat", "lon", "elevation")) %>% 
  rename(el = "elevation")

#1951-1980
dat.1951_80 <- x %>% 
  select(c(catalogNumber,decimalLatitude,decimalLongitude,decade)) %>% 
  group_by(catalogNumber) %>% 
  summarise_all(mean) %>%  
  rename(lat = "decimalLatitude") %>% 
  rename(lon = "decimalLongitude") %>% 
  rename(ID1 = "catalogNumber") 

dat.1951_80 <- dat.1951_80 %>% 
  filter(decade %in% c(1950,1960,1970)) %>% 
  add_column(ID2 = ".", .after = "ID1", ) %>% 
  relocate(lat, .before = lon) %>% 
  select(!decade)

#getting elevation points for each lat/long point
#need projection point for data first 
sf_points <- st_as_sf(dat.1951_80, coords = c("lon","lat"), crs=4326) 

eldat <-  get_elev_point(sf_points) 
dat.1951_80 <- cbind(dat.1951_80, eldat) 

dat.1951_80 <-  dat.1951_80 %>% 
  select(c("ID1", "ID2", "lat", "lon", "elevation")) %>% 
  rename(el = "elevation")

#1981 -2010
dat.1981_2010 <- x %>% 
  select(c(catalogNumber,decimalLatitude,decimalLongitude,decade)) %>% 
  group_by(catalogNumber) %>% 
  summarise_all(mean) %>%  
  rename(lat = "decimalLatitude") %>% 
  rename(lon = "decimalLongitude") %>% 
  rename(ID1 = "catalogNumber") 

dat.1981_2010 <- dat.1981_2010 %>% 
  filter(decade %in% c(1980,1990,2000))%>% 
  add_column(ID2 = ".", .after = "ID1", ) %>% 
  relocate(lat, .before = lon) %>% 
  select(!decade) 

#getting elevation points for each lat/long point
#need projection point for data first 
sf_points <- st_as_sf(dat.1981_2010, coords = c("lon","lat"), crs=4326) 

eldat <-  get_elev_point(sf_points) 
dat.1981_2010 <- cbind(dat.1981_2010, eldat) 

dat.1981_2010 <-  dat.1981_2010 %>% 
  select(c("ID1", "ID2", "lat", "lon", "elevation")) %>% 
  rename(el = "elevation")

#1991 - 2020
dat.1991_2020 <- x %>% 
  select(c(catalogNumber,decimalLatitude,decimalLongitude,decade)) %>% 
  group_by(catalogNumber) %>% 
  summarise_all(mean) %>% 
  rename(lat = "decimalLatitude") %>% 
  rename(lon = "decimalLongitude") %>% 
  rename(ID1 = "catalogNumber") 

dat.1991_2020 <- dat.1991_2020 %>% 
  filter(decade %in% c(1990,2000,2010)) %>% 
  add_column(ID2 = ".", .after = "ID1", ) %>%
  relocate(lat, .before = lon) %>% 
  select(!decade)

#getting elevation points for each lat/long point
#need projection point for data first 
sf_points <- st_as_sf(dat.1991_2020, coords = c("lon","lat"), crs=4326) 

eldat <-  get_elev_point(sf_points) 
dat.1991_2020 <- cbind(dat.1991_2020, eldat) 

dat.1991_2020 <-  dat.1991_2020 %>% 
  select(c("ID1", "ID2", "lat", "lon", "elevation")) %>% 
  rename(el = "elevation")

#no data on 2021
#2022
dat.2022 <- x %>% 
  select(c(catalogNumber,decimalLatitude,decimalLongitude,year.x)) %>% 
  group_by(catalogNumber) %>% 
  summarise_all(mean) %>%  
  rename(lat = "decimalLatitude") %>% 
  rename(lon = "decimalLongitude") %>% 
  rename(ID1 = "catalogNumber") 

dat.2022 <- dat.2022 %>% 
  filter(year.x == 2022) %>% 
  add_column(ID2 = ".", .after = "ID1", ) %>% 
  relocate(lat, .before = lon) %>% 
  select(!year.x)

#getting elevation points for each lat/long point
#need projection point for data first 
sf_points <- st_as_sf(dat.2022, coords = c("lon","lat"), crs=4326) 

eldat <-  get_elev_point(sf_points) 
dat.2022 <- cbind(dat.2022, eldat) 

dat.2022 <-  dat.2022 %>% 
  select(c("ID1", "ID2", "lat", "lon", "elevation")) %>% 
  rename(el = "elevation")

#Climate NA data for each historical bin----
  #change each period argument to match the correct time interval
  #MSY is the timescale argument here. We can use "Y" for yearly/annual, M' for monthly, 'S' for seasonal, 'SY' for annual and seasonal, or 'MSY' for all
clm1901 <- ClimateNAr::ClimateNA_API2(ClimateBC_NA='NA', inputFile=dat.1901_30, period='Normal_1901_1930.nrm',MSY='MSY')
clm1931 <- ClimateNAr::ClimateNA_API2(ClimateBC_NA='NA', inputFile=dat.1931_60, period='Normal_1931_1960.nrm',MSY='MSY')
clm1951 <- ClimateNAr::ClimateNA_API2(ClimateBC_NA='NA', inputFile=dat.1951_80, period='Normal_1951_1980.nrm',MSY='MSY')
clm1981 <- ClimateNAr::ClimateNA_API2(ClimateBC_NA='NA', inputFile=dat.1981_2010, period='Normal_1981_2010.nrm',MSY='MSY')
clm1991 <- ClimateNAr::ClimateNA_API2(ClimateBC_NA='NA', inputFile=dat.1991_2020, period='Normal_1991_2020.nrm',MSY='MSY')

#Can only do 5 pulls per hour so set a 60min timer and then run the next line 
clm2022 <- ClimateNAr::ClimateNA_API2(ClimateBC_NA='NA', inputFile=dat.2022, period='Year_2022.ann',MSY='MSY')


#bind all dataframes together
clmdata <- bind_rows(clm1901, clm1931, clm1951, clm1981, clm1991, clm2022) %>% 
  #distinct(ID1, .keep_all = TRUE) %>%  #removing duplicate rows
  filter(!ID1 == "NA") #for some reason there in an NA catalogNumber that is showing up even though each individual clmX dataset doesn't have it. 
  #It appears to be a duplicate of CHSC075376 with less complete climate data

#checking that 84 observations is correct given the climate data (89 herb. sheets overall but 4 are before 1901)
  #counting the number of sheets within each decade to check the 84 observations
testcount <- x %>% 
  filter(year.x >= 1901) %>% #filtering it by years that are greater or equal to 1901
  distinct(catalogNumber) %>% #grabbing the distinct herbarium sheets
  summarise(n = n()) #counting to double check that the "distinct" on line 197 isn't removing something that we should keep 
testcount


#saving .csv data ----
write.csv(clmdata, file = "Cleaned_data/climate.data.csv")



#Variable description----
# Directly calculated annual variables:
#       MAT	mean annual temperature (°C)
#       MWMT	mean warmest month temperature (°C)
#       MCMT	mean coldest month temperature (°C)
#       TD	temperature difference between MWMT and MCMT, or continentality (°C)
#       MAP	mean annual precipitation (mm)
#       AHM	annual heat-moisture index (MAT+10)/(MAP/1000))
#       SHM	summer heat-moisture index ((MWMT)/(MSP/1000))

# Derived annual variables:
      # DD<0 (or DD_0) 	degree-days below 0°C, chilling degree-days
      # DD>5 (or DD5)	degree-days above 5°C, growing degree-days
      # DD<18 (or DD_18)	degree-days below 18°C, heating degree-days
      # DD>18 (or DD18)	degree-days above 18°C, cooling degree-days
      # NFFD	the number of frost-free days
      # FFP	frost-free period
      # bFFP	the day of the year on which FFP begins
      # eFFP	the day of the year on which FFP ends
      # PAS	precipitation as snow (mm). For individual years, it covers the period between August in the previous year and July in the current year
      # EMT	extreme minimum temperature over 30 years (°C)
      # EXT	extreme maximum temperature over 30 years (°C)
      # Eref	Hargreaves reference evaporation (mm)
      # CMD	Hargreaves climatic moisture deficit (mm)
      # MAR	mean annual solar radiation (MJ m‐2 d‐1)
      # RH	mean annual relative humidity (%)
      # CMI	Hogg’s climate moisture index (mm)
      # DD1040	degree-days above 10°C and below 40°C