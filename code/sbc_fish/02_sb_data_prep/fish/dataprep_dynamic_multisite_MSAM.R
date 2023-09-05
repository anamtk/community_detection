# Ana Miller-ter Kuile
# January 25, 2023
# fish MSOM data prep

# this script preps a multiple sites for a dynamic occupancy model


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

colnames(bs)

#ABUR, AQUE, MOHK

#important variables
# YEAR, MONTH, SITE, TRANSECT, VIS, SP_CODE,
# SIZE, COUNT, AREA

# are all areas the same
unique(fish$AREA) #YES

unique(fish$SP_CODE)

# Clean dataset -----------------------------------------------------------

allfsh <- bs %>%
  filter(!(SITE == "ABUR" & TRANSECT %in% c(2) & YEAR %in% c(2002,
                                                               2003,
                                                               2004,
                                                               2005))) %>%
  dplyr::select(-QUAD, -SIDE)

#get months of resurvey to match yearly all-site survey months
# get just our practice site
# fill visibility across a month for each site
fish1 <- fish %>%
  #get months from yearly surveys only
  filter(MONTH %in% c(7, 8, 9, 10)) %>%
  bind_rows(allfsh) %>%
  #make NA values for visibility
  mutate(VIS = case_when(VIS == -99999 ~ NA_real_,
                         TRUE ~ VIS)) %>%
  group_by(YEAR, SITE,MONTH, TRANSECT) %>%
  #fill missing VIS values for all species on a transect-year-month
  fill(VIS, .direction = "updown") %>%
  ungroup() %>%
  #set rid of NA counts
  mutate(COUNT = case_when(COUNT == -99999 ~ NA_integer_,
                           TRUE ~ COUNT)) %>%
  group_by(YEAR, SITE,  MONTH, TRANSECT, VIS, SP_CODE) %>%
  #get total count per species for all year-months
  summarise(COUNT = sum(COUNT, na.rm = T)) %>%
  #ungroup that dataset
  ungroup() %>%
  unite(SITE_TRANS,
        c("SITE", "TRANSECT"),
        sep = "_",
        remove = F) %>%
  filter(SITE_TRANS != "ABUR_3") %>%
  filter(SITE %in% c("ABUR", "AQUE", "MOHK", "NAPL",
                     "BULL", "GOLB", "IVEE",
                     "SCTW", "SCDI"))

#WORKING: ABUR, AQUE, MOHK, NAPL, BULL, GOLB, IVEE, SCTW, SCDI
#NOT WORKING: CARP, AHND
# 

#some fish are never observed on transects, which in a standard
#MSAM, I would decide to keep them in. However, with the correlation
#structure in the model, this breaks things 
fish1 <- fish1 %>%
  group_by(SP_CODE) %>%
  mutate(max = max(COUNT, na.rm =T)) %>%
  filter(max != 0) %>%
  ungroup()


fish1 %>%
  distinct(SITE_TRANS, YEAR, VIS) %>%
  filter(is.na(VIS)) %>%
  tally()

#29 out of 1135 have missing vis (9%)

#33 out of 1145 have missing vis (9%)

#there are two years with missing months from the survey
#that we want to populate with NA values for the model
# we also want to give repeat months 1:4 IDs so we can
# model them later in JAGS (jags likes numbers)
groups <- fish1 %>%
  #get the distinct combos of year and month
  distinct(SITE_TRANS, YEAR, MONTH) %>%
  #group by year
  group_by(SITE_TRANS, YEAR) %>%
  #arrange in order by month each year
  arrange(MONTH) %>%
  #give each month within a year a 1:4 count
  mutate(REP = 1:n()) %>%
  #ungroup
  ungroup() 

#combine those month IDs and missing months with the 
#fish observation dataset
fish2 <- groups %>%
  full_join(fish1, by = c('SITE_TRANS', 'YEAR', 'MONTH'),
            multiple = "all") %>%
  #make numreic variables for year and species
  mutate(specID = as.numeric(as.factor(SP_CODE)),
         yrID = as.numeric(as.factor(YEAR))) %>%
  #get rid of columns we don't need for the 
  #single site 
  dplyr::select(-SITE, -TRANSECT) 

#for the replicates that don't have any data, we want
# # to give those NA values, that's what this data pull does
# t <- fish2 %>%
#   #select variables of interest
#   dplyr::select(SITE_TRANS, YEAR, yrID, specID, REP) %>%
#   #group by the year
#   group_by(SITE_TRANS, YEAR, yrID) %>%
#   #make sure that the missing replicates for 2002 and 2009
#   # get NA values for all species
#   complete(specID, REP)

#now combine that back wit hteh fisth dataset
fish3 <- fish2 %>%
  #full_join(fish2, by = c("SITE_TRANS", "yrID", "specID", "REP", "YEAR"),
  #          multiple = "all") %>%
  #remove any rows where speciesID is not defined
  filter(!is.na(specID)) %>%
  #group_by(specID) %>%
  #fill the species code for later metadata matching
  #fill(SP_CODE, .direction = "updown") %>%
  #ungroup() %>%
  #missing months were for the years where one survey
  # wasn't completed, both years in month 7, so 
  # we'll add in that data
  #mutate(MONTH = replace_na(MONTH, 7)) %>%
  #get a site id that is numerical
  mutate(siteID = as.numeric(as.factor(SITE_TRANS))) %>%
  #scale visibility covaraite
  mutate(VIS = scale(VIS))  


missing <- fish3 %>%
  dplyr::select(specID, siteID, REP, yrID) %>%
  mutate(specID = as.factor(specID)) %>%
  group_by(siteID,  REP, yrID) %>%
  complete(specID) %>%
  ungroup() %>%
  mutate(specID = as.numeric(specID)) %>%
  filter(!is.na(specID))

new <- missing %>% 
  anti_join(fish3)

fish4 <- fish3 %>%
  full_join(new, by = c("siteID", "REP", "yrID", "specID")) %>%
  group_by(siteID) %>%
  fill(SITE_TRANS) %>%
  ungroup() %>%
  group_by(specID) %>%
  fill(SP_CODE) %>%
  ungroup() %>%
  group_by(yrID) %>%
  fill(YEAR) %>%
  ungroup() %>%
  group_by(SITE_TRANS, REP, YEAR) %>%
  fill(MONTH) %>%
  ungroup() %>%
  group_by(SITE_TRANS, REP, YEAR) %>%
  fill(VIS) %>%
  ungroup()

fish4 %>%
  filter(COUNT ==0) %>%
  tally()

fish4 %>%
  filter(COUNT >0) %>%
  tally()

# Get covariates in order -------------------------------------------------

####Visibility covariate - site, year, visit month

#Now we need to make an array of the observed
#vis data with rows for sites, columns for years,
# and matrix elements for each replicate
nsites <- max(fish4$siteID)#get the dimension of rows 45
nyrs <- max(fish4$yrID) #get the dimension of columns 23
nreps <- max(fish4$REP) #get the dimension of matrices 4

#now, generate IDs for the for loop where 
# we will populate the matrix
yr <- fish4$yrID #get a yearID for each iteration of the loop
site <- fish4$siteID #get a site ID for each iteration of the loop
rep <- fish4$REP #get a replicate for each iteration of the loop

#make a blank array with dims of sites x years x reps
vis <- array(NA, dim = c(nsites, nyrs, nreps))

#fill taht array based on the values in those columns
# for occupancy
for(i in 1:dim(fish4)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the site of row i,
  # the year of row i and the replicate of row i,
  # populate that space in the array with the column in
  # the dataframe that corresponds to the scaled vis data
  # for that sitexyearxreplicate combo
  vis[site[i], yr[i], rep[i]] <- as.numeric(fish4[i,5])
}

### BODY SIZE covariate ###
#get unique species codes to match to the body size 
# dataset
species <- unique(fish4$SP_CODE)

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
#data with rows for species, columns for sites,
# third dimension for years, and fourth dimension for reps in 
#each year
nspecs <- max(fish4$specID) #get the dimension of rows
nsites <- max(fish4$siteID) #get dimension for columns
nyrs <- max(fish4$yrID) #get the dimension of 3rd dimension
nreps <- max(fish4$REP) #get the dimension of  4th dimension

#now, generate IDs for the for loop where 
# we will populate the matrix
yr <- fish4$yrID #get a yearID for each iteration of the loop
site <- fish4$siteID #site ID for each iteration fo the loop
spec <- fish4$specID #get a species ID for each iteration of the loop
rep <- fish4$REP #get a replicate for each iteration of the loop

#make a blank array with dims of species x years x reps
y <- array(NA, dim = c(nspecs, nsites, nyrs, nreps))

#fill taht array based on the values in those columns
# for occupancy
for(i in 1:dim(fish4)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the species of row i,
  #the site of row i, 
  # the year of row i and the replicate of row i,
  # populate that space in the array with the column in
  # the dataframe that corresponds to the 1-0 occupancy
  # for that speciesxyearxreplicate combo
  y[spec[i], site[i], yr[i], rep[i]] <- as.numeric(fish4[i,7])
}

#generate the z-matrix of speciesxsitexyear which assumes
# no false positives
#z matrix is speciesxsitexyear summed over all reps
Ndf <- fish4 %>%
  group_by(siteID, yrID, specID) %>%
  #get total abundance for each sitexspeciesxyear
  summarise(tot = sum(COUNT, na.rm = T)) 


nspecs <- max(Ndf$specID) #get the dimension of rows
nsites <- max(Ndf$siteID) #get dimension for columns
nyrs <- max(Ndf$yrID) #get the dimension of 3rd dimension

#now, generate IDs for the for loop where
# we will populate the matrix
yr <- Ndf$yrID #get a yearID for each iteration of the loop
site <-Ndf$siteID #site ID for each iteration fo the loop
spec <- Ndf$specID #get a species ID for each iteration of the loop

#make a blank array with dims of species x years
ymax <- array(NA, dim = c(nspecs, nsites, nyrs))

#fill taht array based on the values in those columns
# for occupancy
for(i in 1:dim(Ndf)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the species of row i,
  #the site of row i,
  # the year of row i and the replicate of row i,
  # populate that space in the array with the column in
  # the dataframe that corresponds to the 1-0 occupancy
  # for that speciesxyearxreplicate combo
  ymax[spec[i], site[i], yr[i]] <- as.numeric(Ndf[i,4])
}

#set all zeros (which could be true or false) to NA
ymax[ymax == 0] <- NA

# Get covariates and list elements ----------------------------------------

#values for for-loops in the model
n.species <- nrow(y)
n.years <- length(unique(fish4$yrID))
n.transects <- length(unique(fish4$SITE_TRANS))


n.start <- fish4 %>%
  distinct(siteID, yrID) %>%
  group_by(siteID) %>%
  filter(yrID == min(yrID)) %>%
  arrange(siteID) %>%
  ungroup() %>%
  dplyr::select(yrID) %>%
  as_vector()



n.end <- fish4 %>%
  distinct(siteID, yrID) %>%
  group_by(siteID) %>%
  filter(yrID == max(yrID)) %>%
  arrange(siteID) %>%
  ungroup() %>%
  dplyr::select(yrID) %>%
  as_vector()

#site x year matrix
n.rep <- fish4 %>%
  distinct(siteID, yrID, REP) %>%
  group_by(siteID, yrID) %>%
  tally(name = "REP") %>%
  pivot_wider(names_from = yrID,
              values_from = REP) %>%
  column_to_rownames(var = "siteID") %>%
  dplyr::select("1", "2", "3", "4", '5',
                "6", '7', '8', '9', '10',
                '11','12','13','14','15',
                '16','17','18','19','20',
                '21','22', '23') %>%
  as.matrix()


n.rep[which(is.na(n.rep))] <- 1


# Make R covariance matrix ------------------------------------------------

#n.species x n.species matrix of covariance between species abundances
#for the omega parameter prior in the multivariate normal distribution
# this omega will be somewhat of the covariance matrix similar to a
# JSDM 

#R needs to be positive definite, so i think > 0

#will take an average of abundances across all sites and years and 
#then get the covariance matrix from that


# t <- fish4 %>%
#   group_by(yrID, siteID, specID) %>%
#   summarise(COUNT = mean(COUNT, na.rm = T)) %>%
#   ungroup() %>%
#   unite("site_year", c("yrID", "siteID"),
#         sep = "_") %>%
#   dplyr::select(specID, COUNT, site_year) %>%
#   pivot_wider(names_from = specID, 
#               values_from = COUNT,
#               values_fill = 0) %>%
#   column_to_rownames(var = "site_year") %>%
#   mutate(across(everything(), ~replace_na(.x, 0)))
# 
# R <- cor(t)
# 
# library(ggcorrplot)
# 
# ggcorrplot(cor(t), type = "lower",
#            lab = FALSE)

#trying shelby's code
R<-diag(x=0.1, n.species, n.species)

# Make data list to export ------------------------------------------------


#make the data list
data <- list(y = y,
             vis = vis,
             size = sizes,
             n.species = n.species,
             n.years = n.years,
             n.start = n.start,
             n.end = n.end,
             n.transects = n.transects,
             n.rep = n.rep,
             #for initials
             ymax = ymax,
             R = R)

#export that for using with the model
saveRDS(data, here("data_outputs",
                   'sbc_fish',
                   "model_inputs",
                   "fish_msam_dynmultisite.RDS"))



# Export metadata for post summaries --------------------------------------

fish5 <- fish4 %>%
  ungroup() %>%
  distinct(SITE_TRANS, YEAR, siteID, yrID)

write.csv(fish5, here("data_outputs",
                      'sbc_fish',
                      "metadata",
                      "site_year_IDs.csv"),
          row.names = F)

fish6 <- fish4 %>%
  ungroup() %>%
  distinct(SP_CODE, specID)

write.csv(fish6, here("data_outputs",
                      'sbc_fish',
                      "metadata",
                      "species_IDs.csv"),
          row.names = F)

# Raw community matrix ----------------------------------------------------
# 
# 
# 
# #for downstream analyses, we also want the 1-0 matrix for 
# # occupancy of speciesxyear - which we can generate and export
# matrix <- fish4 %>%
#   group_by(YEAR, SP_CODE) %>%
#   summarise(OCC = sum(OCC, na.rm = T)) %>%
#   ungroup() %>%
#   mutate(OCC = case_when(OCC > 0 ~ 1,
#                          OCC ==0 ~ 0,
#                          TRUE ~ NA_real_)) %>%
#   pivot_wider(names_from = "YEAR",
#               values_from = "OCC") %>%
#   column_to_rownames(var = 'SP_CODE') %>%
#   as.matrix()
# 
# #Export that matrix to a central location for all matrices
# saveRDS(matrix, here("data_outputs",
#                      "community_matrices",
#                      "fish_AQUE1_raw_matrix.RDS"))






