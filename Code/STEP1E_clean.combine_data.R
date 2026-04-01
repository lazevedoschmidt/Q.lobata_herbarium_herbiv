#Purpose: To clean and combine datasets
  #I want to combine the occurence data with the data I collected and the 
  #presence/absence dataset
  #It's going to be a massive dataset but it'll have all the info necessary for each
  #analysis down the road
#How the output is used: 
  #Combined dataset for all future uses
#Author: LAS
#Dave: 04/11/2024
  #Date last revised: 08/11/2025
    #Climate data was updated to be the actual annual values as noted below in the code
    #In previous version, the 30-year normals were used and that was a mistake. This has now been corrected. 
  #Date last revised: 10/22/2024
    #notes on revision: Added a step to add in climate variables so that the redundant code can be removed from
    #other scrips
#R version: 4.3.3.

#load packages
require(tidyverse)

#load in datasets
QL_occurence <- read_csv("Raw_data/QL_occurrences.csv") #full (not the subsampled) downloaded data from CCH2 

QL_abund <- read_csv("Cleaned_data/cleaned.abundance_QLmatrix.csv") #abundance matrix from the cleaned folder 
  #can "easily" convert this to pres/ab if needed later

QL_rare <- read_csv("Rare_Output/Per_sheet/Quercus_l. _summary.csv") #there is always a space before the "_X" because of TWs code which is lame but whaterver

QL_rich <- read_csv("Cleaned_data/herb.richness.csv") #this data had to be cleaned once prior in STEP1B.2

ffg_comp <- read_csv("Cleaned_data/ffg_comp.csv") #This is the frequency (%) output from the TW code. In here we now have all the frequency of each FFG + a "chewing" column

##ALL four above datasets were double checked and varified to be correct/accurate on October 3rd 2024
##FFG_comp was created on 10/10/2024

#cleaning data----
#Occurrence CCH2 data
# QL_occurence <-  QL_occurence %>% 
#   select(!"...1")  #removing random index column 

#Abundance data
QL_abund <- QL_abund %>%
  select(!c("...1", "None"))

#Rarefaction data
QL_rare <-  QL_rare %>% 
  select(!"...1") %>% #removing random index column 
  filter(Species != "Total") %>%  #remove the total row from the dataset otherwise it messes everything up
  rename(catalogNumber = Species) #renaming the species column to the catalogNumber so we can join properly below 
    

#richness data
QL_rich <- QL_rich %>% 
  select(!"...1") #removing random index column 

#ffg_comp
ffg_comp <- ffg_comp %>% 
  rename(unique_leaf_ID = ID) %>% 
  select(!"...1")
ffg_comp$unique_leaf_ID <- as.character(ffg_comp$unique_leaf_ID)

#convert percentages of FFG to 0:1 and also convert chewing to 0:1 rather than a sum
ffg_comp2 <- ffg_comp %>% 
  mutate(across(c(Skeletonization, `Surface feeding`, Gall, `Hole Feeding`, 
                  Mining, `Margin Feeding`, Piercing, Chewing), ~ . / 100)) %>% #converting percentages into ratio (0:1)
  mutate(Chewing = if_else(Chewing > 1, 1, Chewing)) #converting all chewing damage into a 0:1

#combining all datasheets----
#joining the df's that are lumped by the catalogNumber column already rather than individual leaves
catNumjoin <- left_join(QL_occurence, QL_rare, by="catalogNumber")  
catNumjoin <- left_join(catNumjoin, QL_rich, by = "catalogNumber")

# 
# QL_fulldf <- QL_herb %>% 
#   #inner_join(QL_occurence, by = c("catalogNumber","unique_leaf_ID","image_leaf_ID","")) %>% 
#   inner_join(QL_abund, by = c("catalogNumber","unique_leaf_ID","image_leaf_ID",""), multiple = "all") %>% 
#   inner_join(QL_rare, by = "catalogNumber", multiple = "all")

#Adding decade column to dataframe
#splitting "date_collected" column into month, day, year columns
QL_abund$date <- as.Date(QL_abund$date_collected,format = "%m/%d/%Y") #creating "date" column in the needed format for splitting
QL_abund$year <- year(ymd(QL_abund$date)) #New year column
QL_abund$month <- month(ymd(QL_abund$date)) #New month column
QL_abund$day <- day(ymd(QL_abund$date)) #New day column 

QL_abund <- QL_abund %>% 
  mutate(decade = floor(year/10)*10) %>%   #creating a "decade" column based on the year
  relocate(year, .after  = date_collected) %>% 
  relocate(month, .after = year) %>% 
  relocate(day, .after = month) %>% 
  relocate(decade, .after = day) %>% 
  relocate(date, .before = year)

#Combining all dataframes together
QL_fulldf <-  left_join(QL_abund, catNumjoin, by = "catalogNumber")
names(QL_fulldf)[names(QL_fulldf) == "%DMG"] <- "perc.dam" #base way of renaming a column since dplyer is hating the "%"
names(QL_fulldf)[names(QL_fulldf) == "%Spec"] <- "perc.spec"
names(QL_fulldf)[names(QL_fulldf) == "%Gall"] <- "perc.gall"
names(QL_fulldf)[names(QL_fulldf) == "%Mine"] <- "perc.mine"


#saving cleaned data files ----
write.csv(QL_abund, "Cleaned_data/QL.cleaneddata.csv", row.names = FALSE)
write.csv(catNumjoin, "Cleaned_data/catNumjoin.csv", row.names = FALSE)
write.csv(QL_fulldf, "Cleaned_data/QL_fulldf.csv", row.names = FALSE)

#Need to add leaf level data to the cleaned data above
leaflevel <- read_csv("Rare_Output/Per_leaf/Quercus_l. _summary.csv") #manually move over the column headers in excel before loading in this data
QL_plantfull <- read_csv("Cleaned_data/QL_fulldf.csv")

#leaf level data
leaflevel <- leaflevel %>% 
  select(!"...1") %>% 
  filter(Species != "Total") %>% 
  rename(unique_leaf_ID = "Species") 

#leaf level
names(leaflevel)[names(leaflevel) == "%DMG"] <- "leaftotal.freq" #base way of renaming a column since dplyer is hating the "%" and "#"
names(leaflevel)[names(leaflevel) == "%Spec"] <- "leafspec.freq"
names(leaflevel)[names(leaflevel) == "%Gall"] <- "leafgall.freq"
names(leaflevel)[names(leaflevel) == "%Mine"] <- "leafmine.freq"
names(leaflevel)[names(leaflevel) == "#leaves"] <- "num.leaves"
names(leaflevel)[names(leaflevel) == "#FFGs"] <- "leafffg.rich"

leaflevel <- leaflevel %>% 
  rename(leaftotal.rich = DTs,
         leafspec.rich = SpecDTs,
         leafgall.rich = GDTs,
         leafmine.rich = MDTs) 

#need to clean the plant data for only the specific columns we need
  #catalog number, unique_leaf_ID, all perc_area, and climate variables 
# QL_plantfull2 <- QL_plantfull %>% 
#   select(c(catalogNumber:decade))
QL_plantfull$unique_leaf_ID <- as.character(QL_plantfull$unique_leaf_ID)

#join the leaf level data with the full dataset above but without the perc.X columns at the end
leaflevel2 <- full_join(QL_plantfull,leaflevel, by="unique_leaf_ID") %>% 
  rename(leafperc_area_chew = perc_area_mandib, #renaming these to double check myself during analysis but ALL perc_area_ are at the leaf level because I took each measurement on individual leaves
         leafperc_area_mine = perc_area_mine,
         leafperc_area_gall = perc_area_gall)

leaflevel3 <- full_join(leaflevel2, ffg_comp2, by = "unique_leaf_ID")
#renaming FFG columns before saving
names(leaflevel3)[names(leaflevel3) == "Skeletonization"] <- "ffg_skel"
names(leaflevel3)[names(leaflevel3) == "Surface feeding"] <- "ffg_surf"
names(leaflevel3)[names(leaflevel3) == "Gall"] <- "ffg_gall"
names(leaflevel3)[names(leaflevel3) == "Hole Feeding"] <- "ffg_hole"
names(leaflevel3)[names(leaflevel3) == "Mining"] <- "ffg_mine"
names(leaflevel3)[names(leaflevel3) == "Margin Feeding"] <- "ffg_marg"
names(leaflevel3)[names(leaflevel3) == "Piercing"] <- "ffg_p"
names(leaflevel3)[names(leaflevel3) == "Chewing"] <- "ffg_chew"

#Double checking that there are no herbarium sheets before 1901
leaflevel4 <- leaflevel3 %>% 
  filter(decade != "1890") %>% 
  filter(decade != "1880") %>% 
  filter(decade != "1870") %>% 
  filter(decade != "1860")


#Adding climate data and decade information to cleaned data
  #This step needs to come after leaflevel3 has been created
clmdata <- read_csv("Cleaned_data/1901_2022_herbariumClimateData.csv") #climate data
  #Updated on 08.11.2025
  #ClimateNA data was pulled from the web version for each individual herbarium sheet
    #The data is NOW the annual variables and not the 30-year normals. 
    #That was a mistake from previous code and has now been updated. As annoying as this is, I'm glad we figured this out now
    #rather than later. 

clmdata_clean <- clmdata %>% distinct()

leaffulldata2 <- clmdata_clean %>%
  #select(!c("...1", "ID2", "ID1...272", "ID1...273")) %>% #removing unneeded columns
  #rename(catalogNumber = "ID1...2") %>% #renaming so we can join the two datasets by catalogNumber
  left_join(leaflevel4, by = "catalogNumber", multiple = "all") #combining climate data with herbivory data

names(leaffulldata2)[names(leaffulldata2) == "#leaves"] <- "sheet.num.leaves"

#Handling leapyears----
# Step 1: 
is_leap_year <- function(year) {
  (year %% 4 == 0 & year %% 100 != 0) | (year %% 400 == 0)
}
#step 2: standarized DOY to 365 days for leap years
standardize_doy <- function(doy, is_leap) {
  ifelse(is_leap, doy * 365/366, doy)
}

# Apply to your data
leaffulldata2$is_leap <- is_leap_year(leaffulldata2$year.x)
leaffulldata2$doy.clean <- standardize_doy(leaffulldata2$startDayOfYear, leaffulldata2$is_leap)

# Verification
cat("Leap years found:", sum(is_leap_year(unique(leaffulldata2$year.x))), "\n")
cat("DOY range - Original:", range(leaffulldata2$startDayOfYear), "\n")
cat("DOY range - Clean:", range(leaffulldata2$doy.clean), "\n")

#saving files
write.csv(leaffulldata2, "Cleaned_data/QL_leaflevel.df.csv")

#Reclassifiying relative leaf age----
#Purpose: To create old and young leaf datasets that are binned by the relative age
#column (leaf_rel_age) and then subset further by the month
#this method should clean our data to use both the relative age but also deal with the
#actual age of the leaf. EX: a leaf that I catagorized as "young" but is from October
#is a much older leaf that a "young" leaf found April. Leaves with a relative age
#of "young" but later in the season will be reclassified as "old"
#How output is used: 
#These cleaned datasets will be used in PCA, and model analysis, as well as any visualization based on age
#Date started: 10.22.2024
#Date last revised: 08/12/2025
#Dataset was updated with the correct climate data and needed to be rerun here to account for the corrections made in the 
#previous scrips. 
#Author: LAS
#R version: 4.3.3.

#load in data
data <- read_csv("Cleaned_data/QL_leaflevel.df.csv") #cleaned dataset at the leaf level

#Subset data by the relative age column
old <- data %>% 
  filter(leaf_rel_age == "old")

young <- data %>% 
  filter(leaf_rel_age == "young")

#subset old data by month
#old == August - December
old2 <- old %>% 
  filter(month.x %in% c("8", "9", "10", "11", "12"))
#removes 60 leaves from the dataset

#young = April (4) - July (7)
young2 <- young %>% 
  filter(month.x %in% c("4", "5", "6", "7"))
#removes 114 leaves but we will move those to the "old" dataset

#Reclassify the lost 114 young leaves to "old"
newold <- young %>% 
  filter(month.x %in% c("1", "2", "3", "8", "9", "10", "11", "12")) %>% #subsetting the young dataset by the old months
  bind_rows(old2) %>% #binding with the subset old data
  mutate(reclas_leaf_rel_age = "old")  #creating a new column where all leaves are labeled as "old"

newold2 <- newold %>% 
  select(!leaf_rel_age) %>% #removing old relative age column
  rename(leaf_rel_age = reclas_leaf_rel_age) #renaming the newly classified column so that future scripts run

#saving output files
write.csv(newold2, "Cleaned_data/reclass.old.leaves.csv")
write.csv(young2, "Cleaned_data/reclass.young.leaves.csv")

#Cleaning all data so this can be removed from the scripts that follow----
oldleaves2 <- newold2%>%
  select(!c("perc_area_PS", "perc_area_ovi")) %>% 
  mutate_at(c("leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall",
              "leaftotal.freq", "leafspec.freq", "leafgall.freq",
              "leafmine.freq"), as.numeric) %>% #R is importing these as characters which isn't correct so updating them here
  mutate(leaf_rel_age = as.factor(leaf_rel_age)) %>% 
  select(c("MAT":"CMI_12", "doy.clean","catalogNumber", "Lat", "long","leaf_rel_age","year.x","leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall",
           "leaftotal.freq", "leafspec.freq", "leafgall.freq",
           "leafmine.freq", "leaftotal.rich", "leafspec.rich", "leafgall.rich", "leafmine.rich")) %>% 
  mutate(across(starts_with("leafperc_area_"), ~ . / 100)) %>% #converting percent leaf area measurements to proportions (0:1)
  mutate(leafgen.freq = leaftotal.freq - leafspec.freq, #creating generalist columns for frequency and richness by supstracting specialist from total
         leafgen.rich = leaftotal.rich - leafspec.rich) 

#creating a new column of data that is 0:1 for chewing damage to then run models on
oldleaves3 <- oldleaves2
oldleaves3$chewbinary <- oldleaves3$leafperc_area_chew 
oldleaves3 <- oldleaves3 %>% 
  mutate(value = if_else(chewbinary < 0.0001, 0, 1)) %>% 
  rename(chewvalue = value) %>% 
  select(!chewbinary)
oldleaves3$chewvalue <- pmax(0.0001, pmin(0.9999, oldleaves3$chewvalue))

# #converting all 0 = 0.0001 and 1 = 0.9999
oldleaves3$leaftotal.freq <- pmax(0.0001, pmin(0.9999, oldleaves3$leaftotal.freq))
oldleaves3$leafgen.freq <- pmax(0.0001, pmin(0.9999, oldleaves3$leafgen.freq))
oldleaves3$leafspec.freq <- pmax(0.0001, pmin(0.9999, oldleaves3$leafspec.freq))
oldleaves3$leafmine.freq <- pmax(0.0001, pmin(0.9999, oldleaves3$leafmine.freq))
oldleaves3$leafgall.freq <- pmax(0.0001, pmin(0.9999, oldleaves3$leafgall.freq))
oldleaves3$leafperc_area_chew <- pmax(0.0001, pmin(0.9999, oldleaves3$leafperc_area_chew))
oldleaves3$leafperc_area_mine <- pmax(0.0001, pmin(0.9999, oldleaves3$leafperc_area_mine))
oldleaves3$leafperc_area_gall <- pmax(0.0001, pmin(0.9999, oldleaves3$leafperc_area_gall))

oldleaves3$leafperc_area_chew <- round(oldleaves3$leafperc_area_chew, 4) #rounding percent area damaged columns to only have 4 decimal places
oldleaves3$leafperc_area_gall <- round(oldleaves3$leafperc_area_gall, 4) #this will make all columns standarized
oldleaves3$leafperc_area_mine <- round(oldleaves3$leafperc_area_mine, 4)
oldleaves3$doy.clean <- round(oldleaves3$doy.clean,2) #Rounding the DOY column since it's been standarized now to deal with leap years

#removing columns with incorrect climate data (i.e., -9999.0)
oldleaves3 <- oldleaves3%>% 
  select(where(~if_else(any(. == -9999.0, na.rm = TRUE), FALSE, TRUE)))

youngleaves2 <- young2%>%
  select(!c("perc_area_PS", "perc_area_ovi")) %>% 
  mutate_at(c("leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall",
              "leaftotal.freq", "leafspec.freq", "leafgall.freq",
              "leafmine.freq"), as.numeric) %>% #R is importing these as characters which isn't correct so updating them here
  mutate(leaf_rel_age = as.factor(leaf_rel_age)) %>% 
  select(c("MAT":"CMI_12","doy.clean","catalogNumber", "Lat", "long","leaf_rel_age","year.x","leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall",
           "leaftotal.freq", "leafspec.freq", "leafgall.freq",
           "leafmine.freq", "leaftotal.rich", "leafspec.rich", "leafgall.rich", "leafmine.rich")) %>% 
  mutate(across(starts_with("leafperc_area_"), ~ . / 100)) %>% #converting percent leaf area measurements to proportions (0:1)
  mutate(leafgen.freq = leaftotal.freq - leafspec.freq, #creating generalist columns for frequency and richness by supstracting specialist from total
         leafgen.rich = leaftotal.rich - leafspec.rich) 

youngleaves3 <- youngleaves2
youngleaves3$chewbinary <- youngleaves3$leafperc_area_chew 
youngleaves3 <- youngleaves3 %>% 
  mutate(value = if_else(chewbinary < 0.0001, 0, 1)) %>% 
  rename(chewvalue = value) %>% 
  select(!chewbinary)
youngleaves3$chewvalue <- pmax(0.0001, pmin(0.9999, youngleaves3$chewvalue))

# #converting all 0 = 0.0001 and 1 = 0.9999
youngleaves3$leaftotal.freq <- pmax(0.0001, pmin(0.9999, youngleaves3$leaftotal.freq))
youngleaves3$leafgen.freq <- pmax(0.0001, pmin(0.9999, youngleaves3$leafgen.freq))
youngleaves3$leafspec.freq <- pmax(0.0001, pmin(0.9999, youngleaves3$leafspec.freq))
youngleaves3$leafmine.freq <- pmax(0.0001, pmin(0.9999, youngleaves3$leafmine.freq))
youngleaves3$leafgall.freq <- pmax(0.0001, pmin(0.9999, youngleaves3$leafgall.freq))
youngleaves3$leafperc_area_chew <- pmax(0.0001, pmin(0.9999, youngleaves3$leafperc_area_chew))
youngleaves3$leafperc_area_mine <- pmax(0.0001, pmin(0.9999, youngleaves3$leafperc_area_mine))
youngleaves3$leafperc_area_gall <- pmax(0.0001, pmin(0.9999, youngleaves3$leafperc_area_gall))

youngleaves3$leafperc_area_chew <- round(youngleaves3$leafperc_area_chew, 4)
youngleaves3$leafperc_area_gall <- round(youngleaves3$leafperc_area_gall, 4)
youngleaves3$leafperc_area_mine <- round(youngleaves3$leafperc_area_mine, 4)
youngleaves3$doy.clean <- round(youngleaves3$doy.clean,2)

youngleaves3<- youngleaves3%>% 
  select(where(~if_else(any(. == -9999.0, na.rm = TRUE), FALSE, TRUE)))

#saving output files
write.csv(oldleaves3, "Cleaned_data/cleaned_reclass.old.leaves.csv")
write.csv(youngleaves3, "Cleaned_data/cleaned_reclass.young.leaves.csv")

