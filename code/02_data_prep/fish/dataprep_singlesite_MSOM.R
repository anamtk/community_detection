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

#fish survey data
fish <- read.csv(here("data_raw",
                      "Monthly_Fish_All_Years_20221018.csv"))

#data on bdoy size distributions for species
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
  
#there are two years with missing months from the survey
#that we want to populate with NA values for the model
# we also want to give repeat months 1:4 IDs so we can
# model them later in JAGS (jags likes numbers)
groups <- fish1 %>%
  #get the distinct combos of year and month
  distinct(YEAR, MONTH) %>%
  #add new months based on the ones missing
  complete(YEAR, MONTH) %>%
  #group by year
  group_by(YEAR) %>%
  #arrange in order by month each year
  arrange(MONTH) %>%
  #give each month within a year a 1:4 count
  mutate(REP = 1:n()) %>%
  #ungroup
  ungroup() 

#combine those month IDs and missing months with the 
#fish observation dataset
fish2 <- groups %>%
  full_join(fish1, by = c('YEAR', 'MONTH')) %>%
  #make numreic variables for year and species
  mutate(specID = as.numeric(as.factor(SP_CODE)),
         yrID = as.numeric(as.factor(YEAR))) %>%
  #make the dataset occupancy
  mutate(OCC = case_when(COUNT > 0~ 1,
                         COUNT == 0 ~ 0,
                         TRUE ~ NA_real_)) %>%
  #get rid of columns we don't need for the 
  #single site 
  dplyr::select(-SITE, -TRANSECT, -COUNT) 

#for the replicates that don't have any data, we want
# to give those NA values, that's what this data pull does
t <- fish2 %>%
  #select variables of interest
  dplyr::select(YEAR, yrID, specID, REP) %>%
  #group by the year
  group_by(YEAR, yrID) %>%
  #make sure that the missing replicates for 2002 and 2009
  # get NA values for all species
  complete(specID, REP)

#now combine that back wit hteh fisth dataset
fish3 <- t %>%
  full_join(fish2, by = c("yrID", "specID", "REP", "YEAR")) %>%
  #remove any rows where speciesID is not defined
  filter(!is.na(specID)) %>%
  group_by(specID) %>%
  #fill the species code for later metadata matching
  fill(SP_CODE, .direction = "updown") %>%
  ungroup() %>%
  #missing months were for the years where one survey
  # wasn't completed, both years in month 7, so 
  # we'll add in that data
  mutate(MONTH = replace_na(MONTH, 7))

# Get covariates in order -------------------------------------------------

####Visibility covariate - by year-month##

#408 observations have NA values for vis
vis <- fish3 %>%
  #get distinct visibility values for each
  #survey month in each year
  distinct(YEAR, REP, VIS) %>%
  #scale this variable
  mutate(VIS = scale(VIS)) %>%
  #in the model, I have the visibility data called 
  # as a matrix with rows of year, columns of replicate
  #within year
  pivot_wider(names_from = "REP",
              values_from = "VIS") %>%
  column_to_rownames(var = "YEAR") %>%
  #make that a matrix
  as.matrix()

#look at the data
hist(vis)

#do we have missing data?
vis[which(is.na(vis))] # 8/82 months have NA for visibility
#we can impute the missing data in the model

### BODY SIZE covariate ###
#get unique species codes to match to the body size 
# dataset
species <- unique(fish3$SP_CODE)

#Select size columns in yearly surveys - by species
sizesa <- bs %>%
  #remove all the NA size values
  filter(SIZE != -99999) %>%
  #select only the variables of interest
  dplyr::select(SIZE, SP_CODE, COUNT)

#again, with the monthly data, do the same thing
#select the size data
sizes <- fish %>%
  #remove the missing values for size
  filter(SIZE != -99999) %>%
  #select only variables of interest
  dplyr::select(SIZE, SP_CODE, COUNT) %>%
  #add in the yearly surveys (They seem like
  # they're on different days)
  bind_rows(sizesa) %>%
  #group by species ID
  group_by(SP_CODE) %>%
  #get a weighted average based on the distributions
  # of sizes of individuals observed in taht size
  # class
  summarise(AVG_SIZE = weighted.mean(SIZE, COUNT)) %>%
  #select only the species in the species codes in the dataset
  filter(SP_CODE %in% species) %>%
  #scle the variable so the model plays nice
  mutate(AVG_SIZE = scale(AVG_SIZE)) %>%
  dplyr::select(AVG_SIZE) %>%
  #make this variable a vector
  as_vector()

#look at that variable's distribution
hist(sizes)


# Prep observed data ------------------------------------------------------

#Now we need to make an array of the observed
#data with rows for species, columns for years,
# and matrix elements for each replicate
nspecs <- max(fish3$specID) #get the dimension of rows
nyrs <- max(fish3$yrID) #get the dimension of columns
nreps <- 4 #get the dimension of matrices

#now, generate IDs for the for loop where 
# we will populate the matrix
yr <- fish3$yrID #get a yearID for each iteration of the loop
spec <- fish3$specID #get a species ID for each iteration of the loop
rep <- fish3$REP #get a replicate for each iteration of the loop

#make a blank array with dims of species x years x reps
y <- array(NA, dim = c(nspecs, nyrs, nreps))

#fill taht array based on the values in those columns
# for occupancy
for(i in 1:dim(fish3)[1]){ #dim[1] = n.speies
  #using info from the dataframe on the species of row i,
  # the year of row i and the replicate of row i,
  # populate that space in the array with the column in
  # the dataframe that corresponds to the 1-0 occupancy
  # for that speciesxyearxreplicate combo
  y[spec[i], yr[i], rep[i]] <- as.numeric(fish3[i,7])
}

#generate the z-matrix of speciesxyear which assumes
# no false positives
#z matrix is year and species summed over all reps
z <- fish3 %>% 
  group_by(yrID, specID) %>%
  #get total occupancy for each speciesxyear
  summarise(tot_occ = sum(OCC, na.rm = T)) %>%
  #set that to 1-0 values again
  mutate(tot_occ = case_when(tot_occ > 0 ~ 1,
                             tot_occ == 0 ~ 0,
                             TRUE ~ NA_real_)) %>%
  #make this into a matrix that is species x year
  #filled with the 1-0 occupancies
  pivot_wider(names_from = 'yrID',
              values_from = "tot_occ") %>%
  #set the rownames to be species
  column_to_rownames(var = "specID") %>%
  #make it a matrix
  as.matrix()

#set all zeros (which could be true or false) to NA
z[z == 0] <- NA
# Get covariates and list elements ----------------------------------------

#values for for-loops in the model
n.species <- nrow(y)
n.years <- ncol(y)
n.reps <- 4

#make the data list
data <- list(y = y,
             z = z,
             vis = vis,
             size = sizes,
             n.species = n.species,
             n.years = n.years,
             n.reps = n.reps)

#export that for using with the model
saveRDS(data, here("data_outputs",
                   "model_inputs",
                   "fish_data_singlesite.RDS"))


# Raw community matrix ----------------------------------------------------

#for downstream analyses, we also want the 1-0 matrix for 
# occupancy of speciesxyear - which we can generate and export
matrix <- fish3 %>%
  group_by(YEAR, SP_CODE) %>%
  summarise(OCC = sum(OCC, na.rm = T)) %>%
  ungroup() %>%
  mutate(OCC = case_when(OCC > 0 ~ 1,
                         OCC ==0 ~ 0,
                         TRUE ~ NA_real_)) %>%
  pivot_wider(names_from = "YEAR",
              values_from = "OCC") %>%
  column_to_rownames(var = 'SP_CODE') %>%
  as.matrix()

#Export that matrix to a central location for all matrices
saveRDS(matrix, here("data_outputs",
                       "community_matrices",
                       "fish_AQUE1_raw_matrix.RDS"))





