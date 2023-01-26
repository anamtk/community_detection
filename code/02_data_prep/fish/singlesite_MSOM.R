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
  mutate(REP = 1:n()) %>%
  ungroup() 

fish2 <- groups %>%
  full_join(fish1, by = c('YEAR', 'MONTH')) %>%
  #make numreic variables for year and species
  mutate(specID = as.numeric(as.factor(SP_CODE)),
         yrID = as.numeric(as.factor(YEAR))) %>%
  #make the dataset occupancy
  mutate(OCC = case_when(COUNT > 0~ 1,
                         COUNT == 0 ~ 0,
                         TRUE ~ NA_real_)) %>%
  dplyr::select(-SITE, -TRANSECT, -COUNT) 

t <- fish2 %>%
  dplyr::select(YEAR, yrID, specID, REP) %>%
  group_by(YEAR, yrID) %>%
  complete(specID, REP)

fish3 <- t %>%
  full_join(fish2, by = c("yrID", "specID", "REP", "YEAR")) %>%
  filter(!is.na(specID)) %>%
  group_by(specID) %>%
  fill(SP_CODE, .direction = "updown") %>%
  ungroup() %>%
  mutate(MONTH = replace_na(MONTH, 7))

#LOTS O ZEROS YO 
#4182-268 #= 3914 0's' - 94%!!

# Get covariates in order -------------------------------------------------
#Visibility covariate - by year-month

#408 observations have NA values for vis
vis <- fish3 %>%
  distinct(YEAR, REP, VIS) %>%
  mutate(VIS = scale(VIS)) %>%
  pivot_wider(names_from = "REP",
              values_from = "VIS") %>%
  column_to_rownames(var = "YEAR") %>%
  as.matrix()

hist(vis)

vis[which(is.na(vis))] # 8/82 months have NA for visibility

species <- unique(fish3$SP_CODE)

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
rep <- fish3$REP

#make a blank array
y <- array(NA, dim = c(nspecs, nyrs, nreps))

#STILL NEED TO SOLVE THiS - it's NA for all reps
# for 2002 and 2009 and should be for just one of the reps - not
# sure which one - but i'm closer!

for(i in 1:dim(fish3)[1]){
  y[spec[i], yr[i], rep[i]] <- as.numeric(fish3[i,7])
}


write.csv(y, 
          here("data_outputs",
               "raw_community",
               "single_site_MSOM_matrix.csv"))


#z matrix is year and species summed over all reps
z <- fish3 %>% 
  group_by(yrID, specID) %>%
  summarise(tot_occ = sum(OCC, na.rm = T)) %>%
  mutate(tot_occ = case_when(tot_occ > 0 ~ 1,
                             tot_occ == 0 ~ 0,
                             TRUE ~ NA_real_)) %>%
  pivot_wider(names_from = 'yrID',
              values_from = "tot_occ") %>%
  column_to_rownames(var = "specID") %>%
  as.matrix()

z[z == 0] <- NA
# Get covariates and list elements ----------------------------------------

n.species <- nrow(y)
n.years <- ncol(y)
n.reps <- 4

data <- list(y = y,
             z = z,
             vis = vis,
             size = sizes,
             n.species = n.species,
             n.years = n.years,
             n.reps = n.reps)

saveRDS(data, here("data_outputs",
                   "model_inputs",
                   "fish_data_singlesite.RDS"))
