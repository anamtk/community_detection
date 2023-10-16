# Ana Miller-ter Kuile
# July 27, 2023
# konza bird data prep 

# this script preps a multiple sites for a dynamic abundance model


# Load packages -----------------------------------------------------------


package.list <- c("here", "tidyverse",
                  'ggcorrplot')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

birds <- read.csv(here("02_konza_birds",
                       "data_raw",
                       "CBP011.csv"))

str(birds)
#for body size
#Supplementary Data 1 here: https://figshare.com/s/b990722d72a26b5bfead
avonet <- read.csv(here('02_konza_birds',
                        'data_raw',
                        'AVONET_BirdLife.csv'))

avonet2 <- read.csv(here('02_konza_birds',
                         'data_raw',
                        'AVONET_BirdTree.csv'))

avonet3 <- read.csv(here('02_konza_birds',
                         'data_raw',
                         'AVONET_eBird.csv'))

#to link the konza birds via codes and scientific names to the avonet data
#from here: https://www.birdpop.org/pages/birdSpeciesCodes.php
codes <- read.csv(here('02_konza_birds',
                       'data_raw',
                       'IBP-AOS-LIST23.csv'))

# Filter out watersheds of interest ---------------------------------------

#only going to do watersheds where we also have NPP data
#which are 001D, 004B, and 020B

birds1 <- birds %>%
  filter(WATERSHED %in% c("001D", "004B", "020B"))


# Match up species codes to scientific names ------------------------------

#bird size
konz_codes <- unique(birds$AOUCODE)

#get the codes from IBP that only 
#match Konza birds
codes1 <- codes %>%
  filter(SPEC %in% konz_codes) %>%
  dplyr::select(-SP, -CONF, -SPEC6, -CONF6)

#what's funky?
bird_id <- birds1 %>%
  distinct(AOUCODE, SPECNAME, COMMONNAME) %>%
  filter(SPECNAME != "NONE") %>%
  left_join(codes1, by = c("AOUCODE" = "SPEC"))

#some things not IDed to species, so we want to remove those
#either / or Empidonax
#unique ID is actually common name, not SPECNAME

#clean up so that things will mesh well
birds2 <- birds1 %>%
  mutate(COMMONNAME = case_when(COMMONNAME == "Le Conte's Sparrow" ~ "LeConte's Sparrow",
                              TRUE ~ COMMONNAME)) %>%
  filter(COMMONNAME != "Empidonax sp.") %>%
  filter(!str_detect(COMMONNAME, "/"))  %>%
  filter(COMMONNAME != "sparrow sp.")

#Get unique Bird IDS with common name, AOU code and scientific names
bird_id2 <- birds2 %>%
  distinct(COMMONNAME, SPECNAME, AOUCODE) %>%
  left_join(codes1, by = c("COMMONNAME")) %>%
  #American Goldfinch is acting up and I don't know why....!!! GAH
  mutate(SPEC = case_when(COMMONNAME == "American Goldfinch" ~ "AGOL",
                          TRUE ~ SPEC)) %>%
  mutate(SCINAME = case_when(COMMONNAME == "American Goldfinch" ~ "Spinus tristis",
                             TRUE ~ SCINAME)) %>%
  filter(SPECNAME != "NONE") %>%
  distinct(COMMONNAME, SPEC, SCINAME) 
#71 total bird species


# Manupulate data to abundance structure ----------------------------------

#Data are individuals seen at different distances from a line transect

#need to summarise data as count data for each recording
#period in each year for each species

birds3 <- birds2 %>%
  group_by(RECYEAR, RECMONTH, RECDAY, TRANSNUM,
           WATERSHED, COMMONNAME) %>%
  #get number observed
  summarise(NOBS = sum(COUNT)) %>%
  ungroup()
  
#i'm not sure if every birds is documented in each survey period, so
#might need to fill this out with the zeros

birds4 <- birds3 %>%
  #unique identifier is common name
  mutate(COMMONNAME = as.factor(COMMONNAME)) %>%
  #group by all the survey ID info
  group_by(RECYEAR, RECMONTH, RECDAY, TRANSNUM, WATERSHED) %>%
  #make sure each species is in each survey
  complete(COMMONNAME) %>%
  ungroup() %>%
  #get rid of "NONE" species category
  filter(COMMONNAME != "No birds detected") %>%
  #set all NA counts to 0s
  mutate(NOBS = case_when(is.na(NOBS) ~ 0,
                          TRUE ~ NOBS)) %>%
  #re-order the factor of species now that NONE is gone
  mutate(COMMONNAME = as.character(COMMONNAME))%>%
  mutate(TransID = as.numeric(as.factor(TRANSNUM))) %>%
  mutate(yrID = as.numeric(as.factor(RECYEAR))) %>%
  mutate(SpecID = as.numeric(as.factor(COMMONNAME))) 


#no missing years in the time series
birds4 %>%
  distinct(RECYEAR, TRANSNUM) %>%
  group_by(TRANSNUM) %>%
  tally()



# Covariates for detection ------------------------------------------------

#get info from bird dataset to blend in here:
size_meta <- birds4 %>%
  distinct(COMMONNAME, SpecID)
#body size
sizes <- bird_id2 %>%
  #combining with EBird data because the sceitnfiic names
  #match up beautifully and no having to change them
  left_join(avonet3, by = c("SCINAME" = "Species2")) %>%
  dplyr::select(COMMONNAME, SCINAME, Mass) %>%
  #get species ID for modeling
  left_join(size_meta, by = "COMMONNAME")

#get info from bird dataset on ids of different surveys to blend to
#survey dataset
survey_meta <- birds4 %>%
  distinct(RECYEAR, RECMONTH, RECDAY, TRANSNUM, WATERSHED,
           TransID, yrID)

#survey length
survey <- birds %>%
  filter(WATERSHED %in% c("001D", "004B", "020B")) %>%
  distinct(RECYEAR, RECMONTH, RECDAY, TRANSNUM, 
           WATERSHED, DURATION) %>%
  filter(!is.na(DURATION)) %>%
  #zero durations are duplicates of other surveys, so removing them
  filter(DURATION > 0)


#add survey lengths to bird dataset
survey2 <- birds4 %>%
  left_join(survey, by = c("RECYEAR", "RECMONTH", "RECDAY", "TRANSNUM",
                           "WATERSHED")) %>%
  distinct(RECYEAR, RECMONTH, RECDAY, TRANSNUM, WATERSHED, DURATION,
           TransID, yrID) %>% 
  group_by(TRANSNUM, RECYEAR) %>%
  arrange(RECMONTH) %>%
  mutate(REP = 1:n()) %>%
  ungroup() 


# Get Rep data in birds dataset -------------------------------------------

birds5 <- birds4 %>%
  left_join(survey2, by = c("RECYEAR", "RECMONTH", "RECDAY", 
                            "TRANSNUM", "WATERSHED", "TransID",
                            "yrID"))

# Prep data structure for JAGS --------------------------------------------

n.species <- length(unique(birds5$COMMONNAME))

n.transects <- length(unique(birds5$TRANSNUM))

years <- birds5 %>%
  distinct(RECYEAR) %>%
  mutate(yearnum = 1:n())

n.start <- birds5 %>%
  ungroup() %>%
  distinct(TRANSNUM, RECYEAR) %>%
  left_join(years, by = "RECYEAR") %>%
  arrange(TRANSNUM, yearnum) %>%
  group_by(TRANSNUM) %>%
  filter(yearnum == min(yearnum)) %>%
  ungroup() %>%
  dplyr::select(yearnum) %>%
  as_vector()

n.end <- birds5 %>%
  ungroup() %>%
  distinct(TRANSNUM, RECYEAR) %>%
  left_join(years, by = "RECYEAR") %>%
  arrange(TRANSNUM, yearnum) %>%
  group_by(TRANSNUM) %>%
  filter(yearnum == max(yearnum)) %>%
  ungroup() %>%
  dplyr::select(yearnum) %>%
  as_vector()

n.rep <- birds5 %>%
  distinct(TRANSNUM, RECYEAR, RECMONTH, RECDAY) %>%
  group_by(TRANSNUM, RECYEAR) %>%
  tally() %>%
  pivot_wider(names_from = RECYEAR,
              values_from = n) %>%
  column_to_rownames(var = "TRANSNUM") %>%
  as.matrix()

n.years <- length(unique(birds5$RECYEAR))



# Get covariates in order -------------------------------------------------

#Size Covariate
size <- sizes %>%
  mutate(Mass = scale(Mass)) %>%
  dplyr::select(Mass, SpecID) %>%
  arrange(SpecID) %>%
  column_to_rownames(var = "SpecID") %>%
  as_vector()

####Effort covariate - site, year, visit month
#Now we need to make an array of the observed
#effort data with rows for sites, columns for years,
# and matrix elements for each replicate

#now, generate IDs for the for loop where 
# we will populate the matrix
yr <- survey2$yrID
site <- survey2$TransID  
rep <- survey2$REP

survey3 <- survey2 %>%
  mutate(effort = scale(DURATION))

#make a blank array with dims of sites x years x reps
effort <- array(NA, dim = c(n.transects, #rows
                         n.years, #columns
                         4 #array elements - number of replicates
                         ))

#fill taht array based on the values in those columns
# for occupancy
for(i in 1:dim(survey2)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the site of row i,
  # the year of row i and the replicate of row i,
  # populate that space in the array with the column in
  # the dataframe that corresponds to the scaled vis data
  # for that sitexyearxreplicate combo
  effort[site[i], yr[i], rep[i]] <- as.numeric(survey3[i,10])
}

#now, generate IDs for the for loop where 
# we will populate the matrix
yr <- birds5$yrID #get a yearID for each iteration of the loop
site <- birds5$TransID #site ID for each iteration fo the loop
spec <- birds5$SpecID #get a species ID for each iteration of the loop
rep <- birds5$REP #get a replicate for each iteration of the loop

y <- array(NA, dim = c(n.species, #rows
                       n.transects, #column
                       n.years, #first array level
                       4 #second array level
))

#fill taht array based on the values in those columns
# for occupancy
for(i in 1:dim(birds5)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the species of row i,
  #the site of row i, 
  # the year of row i and the replicate of row i,
  # populate that space in the array with the column in
  # the dataframe that corresponds to the 1-0 occupancy
  # for that speciesxyearxreplicate combo
  y[spec[i], site[i], yr[i], rep[i]] <- as.numeric(birds5[i,7])
}



#generate the y-max values
Ndf <- birds5 %>%
  group_by(TransID, yrID, SpecID) %>%
  #get total abundance for each sitexspeciesxyear
  summarise(tot = sum(NOBS, na.rm = T)) 


nspecs <- max(Ndf$SpecID) #get the dimension of rows
nsites <- max(Ndf$TransID) #get dimension for columns
nyrs <- max(Ndf$yrID) #get the dimension of 3rd dimension

#now, generate IDs for the for loop where
# we will populate the matrix
yr <- Ndf$yrID #get a yearID for each iteration of the loop
site <-Ndf$TransID #site ID for each iteration fo the loop
spec <- Ndf$SpecID #get a species ID for each iteration of the loop

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


# Make R covariance matrix ------------------------------------------------

#n.species x n.species matrix of covariance between species abundances
#for the omega parameter prior in the multivariate normal distribution
# this omega will be somewhat of the covariance matrix similar to a
# JSDM 

#R needs to be positive definite,

#trying shelby's code from Kiona/Jessica - need to ask what this means
R<-diag(x=0.1, n.species, n.species)

#omega also needs priors, which I'm going to attempt to define using
#covariance among species abundances, we'll see how it goes

t <- birds5 %>%
  group_by(yrID, TransID, SpecID) %>%
  summarise(COUNT = mean(NOBS, na.rm = T)) %>%
  ungroup() %>%
  unite("site_year", c("yrID", "TransID"),
        sep = "_") %>%
  dplyr::select(SpecID, COUNT, site_year) %>%
  pivot_wider(names_from = SpecID,
              values_from = COUNT,
              values_fill = 0) %>%
  column_to_rownames(var = "site_year") %>%
  mutate(across(everything(), ~replace_na(.x, 0)))
# 
# 
# ggcorrplot(cor(t), type = "lower",
#            lab = FALSE)

#set omega init to this - not sure if it will work with the NA values
#or if i will need to define those as a value?? we can try it...
omega.init <- cor(t)

# Prep list for JAGS ------------------------------------------------------

data <- list(n.species = n.species,
             n.transects = n.transects,
             n.start = n.start,
             n.end = n.end,
             n.rep = n.rep,
             n.years = n.years,
             effort = effort,
             size = size,
             y = y,
             #omega prior:
             R = R,
             #initials
             ymax = ymax,
             omega.init = omega.init)


#export that for using with the model
saveRDS(data, here('02_konza_birds',
                   "data_outputs",
                   "MSAM",
                   "model_inputs",
                   "bird_msam_dynmultisite.RDS"))

# Uncorrected bray for all sites-years ------------------------------------

#pulling out and calculating bray-curtis for all 
#communities so we can compare to the "corrected" version
#after the model

diss_fun <- function(site){
  
  matrix <- birds5 %>%
    filter(TransID == site) %>%
    group_by(yrID, SpecID) %>%
    summarise(COUNT = max(NOBS, na.rm = T)) %>%
    pivot_wider(names_from = yrID,
                values_from = COUNT) %>%
    column_to_rownames(var = 'SpecID')
  
  a <- matrix(NA, nrow = nrow(matrix),
              ncol = ncol(matrix))
  
  b <- matrix(NA, nrow = nrow(matrix),
              ncol = ncol(matrix))
  
  c <- matrix(NA, nrow = nrow(matrix),
              ncol = ncol(matrix))
  
  for(r in 1:nrow(matrix)){
    for(t in 2:ncol(matrix)){
      a[r, t] <- min(c(matrix[r,t-1], matrix[r,t]))
      b[r,t] <- matrix[r,t-1] - a[r,t]
      c[r,t] <- matrix[r,t] - a[r,t]
    }
  }
  
  A <- colSums(a)
  B <- colSums(b)
  C <- colSums(c)
  
  bray <- (B + C)/(2*A+B+C)
  
  bray_df <- birds5 %>%
    filter(TransID == site) %>%
    distinct(yrID, TransID) %>%
    cbind(bray) %>%
    mutate(type = "observed")
  
  return(bray_df)
}

sites <- unique(birds5$TransID)
results <- lapply(sites, FUN = diss_fun)

results_df <- do.call(rbind, results)

saveRDS(results_df, here('05_visualizations',
                         'viz_data',
                         'konza_observed_bray.RDS'))


# Export metadata ---------------------------------------------------------

sites <- birds5 %>%
  distinct(RECYEAR, TRANSNUM,
           WATERSHED,
           TransID, yrID)

write.csv(sites, here('02_konza_birds',
               'data_outputs',
               'metadata',
               'site_year_IDs.csv'))
