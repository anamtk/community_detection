# Ana Miller-ter Kuile
# January 25, 2023
# fish MSOM data prep

# this script preps a single site of repeat fish surveys for 
# a preliminary MSOM analysis


# Load packages -----------------------------------------------------------


package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

fish <- read.csv(here("data_raw",
                      "Monthly_Fish_All_Years_20221018.csv"))

bs <- read.csv(here("data_raw",
                    "Annual_fish_comb_20220809.csv"))



# Explore fish ------------------------------------------------------------

colnames(fish)


#important variables
# YEAR, MONTH, SITE, TRANSECT, VIS, SP_CODE,
# SIZE, COUNT, AREA

# are all areas the same
unique(fish$AREA) #YES

unique(fish$SP_CODE)


# Clean dataset -----------------------------------------------------------

#get months of resurvey to match yearly all-site survey months
# get just our practice site
# fill visibility across a month for each site
fish1 <- fish %>%
  #get months from yearly surveys only
  filter(MONTH %in% c(7, 8, 9, 10)) %>%
  #select our initial test site
  filter(SITE == "AQUE" & TRANSECT == 1) %>%
  #make NA values for visibility
  mutate(VIS = case_when(VIS == -99999 ~ NA_real_,
                         TRUE ~ VIS)) %>%
  group_by(YEAR, SITE, MONTH, TRANSECT) %>%
  #fill missing VIS values for all species on a transect-year-month
  fill(VIS, .direction = "updown") %>%
  ungroup() %>%
  #set rid of NA counts
  mutate(COUNT = case_when(COUNT == -99999 ~ NA_integer_,
                           TRUE ~ COUNT)) %>%
  group_by(YEAR, SITE, MONTH, TRANSECT, VIS, SP_CODE) %>%
  #get total count per species for all year-months
  summarise(COUNT = sum(COUNT, na.rm = T)) %>%
  #ungroup that dataset
  ungroup() 
  

groups <- fish1 %>%
  distinct(YEAR, MONTH) %>%
  complete(YEAR, MONTH) %>%
  group_by(YEAR) %>%
  arrange(MONTH) %>%
  mutate(REP = 1:n())

fish2 <- fish1 %>%
  full_join(groups, by = c('YEAR', 'MONTH')) %>%
  #make numreic variables for year and species
  mutate(specID = as.numeric(as.factor(SP_CODE)),
         yrID = as.numeric(as.factor(YEAR))) %>%
  #make the dataset occupancy
  mutate(OCC = case_when(COUNT > 0~ 1,
                         COUNT == 0 ~ 0,
                         TRUE ~ NA_real_)) %>%
  dplyr::select(-YEAR, -SP_CODE, -SITE, -TRANSECT,
                -VIS, -COUNT, -MONTH) 

t <- fish2 %>%
  dplyr::select(yrID, specID, REP) %>%
  filter(!is.na(specID)) %>%
  complete(yrID, specID, REP)

fish3 <- fish2 %>%
  full_join(t, by = c("yrID", "specID", "REP")) %>%
  filter(!is.na(specID))

#LOTS O ZEROS YO 
#4182-268 #= 3914 0's' - 94%!!


# Get covariates in order -------------------------------------------------
#Visibility covariate - by year-month

#408 observations have NA values for vis
vis <- fish2 %>%
  distinct(YEAR, MONTH, VIS) %>%
  dplyr::select(VIS) %>%
  mutate(VIS = scale(VIS)) %>%
  as_vector()

hist(vis)

vis[which(is.na(vis))] # 8/82 months have NA for visibility

#408/4402 #NA for visibility


species <- unique(fish2$SP_CODE)

#Average sizes - by species
sizesa <- bs %>%
  filter(SIZE != -99999) %>%
  dplyr::select(SIZE, SP_CODE, COUNT)

sizes <- fish %>%
  filter(SIZE != -99999) %>%
  dplyr::select(SIZE, SP_CODE, COUNT) %>%
  bind_rows(sizesa) %>%
  group_by(SP_CODE) %>%
  summarise(AVG_SIZE = weighted.mean(SIZE, COUNT)) %>%
  filter(SP_CODE %in% species) %>%
  mutate(AVG_SIZE = scale(AVG_SIZE)) %>%
  dplyr::select(AVG_SIZE) %>%
  as_vector()

hist(sizes)


# Prep observed data ------------------------------------------------------

nspecs <- max(fish3$specID)
nyrs <- max(fish3$yrID)
nreps <- 4

#need to convert these to numeric in the dataframe
yr <- fish3$yrID
spec <- fish3$specID

#make a blank array
y <- array(NA, dim = c(nspecs, nyrs, nreps))

#STILL NEED TO SOLVE THiS - it's NA for all reps
# for 2002 and 2009 and should be for just one of the reps - not
# sure which one - but i'm closer!

for(i in 1:dim(fish3)[1]){
  y[spec[i], yr[i], 1:4] <- as.numeric(fish3[i,4])
}

fish2 %>%
  filter(REP == 4) %>%
  group_by(YEAR, REP) %>%
  summarise(total = sum(OCC))

write.csv(y, 
          here("data_outputs",
               "raw_community",
               "single_site_MSOM_matrix.csv"))

z <- (y>0)*1
z[z == 0] <- NA
# Get covariates and list elements ----------------------------------------

n.species <- ncol(y)
n.years <- nrow(y)
n.reps <- 4

