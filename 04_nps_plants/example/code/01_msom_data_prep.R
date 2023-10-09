# Ana Miller-ter Kuile
# July 27, 2023
# konza bird data prep 

# this script preps a multiple sites for a dynamic occupancy model


# Load packages -----------------------------------------------------------


package.list <- c("here", "tidyverse",
                  'ggcorrplot')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

#Get data from the internet (don't commit this to Github, it's too big!)
url <- 'http://lter.konza.ksu.edu/sites/default/files/data/CBP011.csv'

dest <- here('04_nps_plants',
             'example',
             'data_raw',
             'raw_Konza_birds.csv')

download.file(url = url, destfile = dest)

#read that data into R
birds <- read.csv(here('04_nps_plants',
                       'example',
                       'data_raw',
                       'raw_Konza_birds.csv'))

str(birds)


# Manipulate data to occupancy structure ----------------------------------

#Data are individuals seen at different distances from a line 
#transect. Might be worth getting some size data from AVONET for this
#dataset, making it consistent with the fish dataset? 

#need to summarise data as 1-0 for each
#recording period in each year for each species.
occ <- birds %>%
  group_by(RECYEAR, RECMONTH, RECDAY, TRANSNUM,
           WATERSHED, SPECNAME, COMMONNAME) %>%
  summarise(total = sum(COUNT)) %>%
  ungroup()

#Transect is the lowest ID, NO1B is both Transect 6 and 10
occ %>%
  distinct(TRANSNUM, WATERSHED)

#Get survey covariate of effort
survey <- birds %>%
  distinct(RECYEAR, RECMONTH, RECDAY, TRANSNUM, OBSERVER, DURATION,
           TIME) %>%
  filter(!is.na(DURATION)) %>%
  mutate(DURATION = case_when(DURATION == 0 ~ NA_real_,
                              TRUE ~ DURATION)) %>%
  group_by(TRANSNUM, RECYEAR) %>%
  arrange(RECMONTH) %>%
  mutate(REP = 1:n()) %>%
  ungroup() %>%
  mutate(yrID = as.numeric(as.factor(RECYEAR))) %>%
  mutate(TransID = as.numeric(as.factor(TRANSNUM)))

rep <- survey %>%
  distinct(yrID, TransID, RECMONTH, REP)

occ2 <- occ %>%
  mutate(SPECNAME = as.factor(SPECNAME)) %>%
  #factor species name so we can complete it by survey interval
  dplyr::select(-COMMONNAME) %>%
  #group by all the survey ID info
  group_by(RECYEAR, RECMONTH, RECDAY, TRANSNUM, WATERSHED) %>%
  #make sure each species is in each survey
  complete(SPECNAME) %>%
  #presence: 1 when there was at least one detected, 0 if none detected
  mutate(presence = case_when(!is.na(total) ~ 1,
                              TRUE ~ 0)) %>%
  ungroup() %>%
  #get rid of "NONE" species category
  filter(SPECNAME != "NONE") %>%
  #remove unecessary columns
  dplyr::select(-total) %>%
  #re-order the factor of species now that NONE is gone
  mutate(SPECNAME = as.factor(as.character(SPECNAME)))%>%
  mutate(TransID = as.numeric(as.factor(TRANSNUM))) %>%
  mutate(yrID = as.numeric(as.factor(RECYEAR))) %>%
  mutate(SpecID = as.numeric(as.factor(SPECNAME))) %>%
  left_join(rep, by = c("yrID", "TransID", "RECMONTH"))
  
         
#no missing years in the time series
occ2 %>%
  distinct(RECYEAR, TRANSNUM) %>%
  group_by(TRANSNUM) %>%
  tally()



# Prep data structure for JAGS --------------------------------------------

n.species <- length(unique(occ2$SPECNAME))

n.transects <- length(unique(occ2$TRANSNUM))

years <- occ2 %>%
  distinct(RECYEAR) %>%
  mutate(yearnum = 1:n())

n.start <- occ2 %>%
  ungroup() %>%
  distinct(TRANSNUM, RECYEAR) %>%
  left_join(years, by = "RECYEAR") %>%
  arrange(TRANSNUM, yearnum) %>%
  group_by(TRANSNUM) %>%
  filter(yearnum == min(yearnum)) %>%
  ungroup() %>%
  dplyr::select(yearnum) %>%
  as_vector()

n.end <- occ2 %>%
  ungroup() %>%
  distinct(TRANSNUM, RECYEAR) %>%
  left_join(years, by = "RECYEAR") %>%
  arrange(TRANSNUM, yearnum) %>%
  group_by(TRANSNUM) %>%
  filter(yearnum == max(yearnum)) %>%
  ungroup() %>%
  dplyr::select(yearnum) %>%
  as_vector()

n.rep <- occ2 %>%
  distinct(TRANSNUM, RECYEAR, RECMONTH, RECDAY) %>%
  group_by(TRANSNUM, RECYEAR) %>%
  tally() %>%
  pivot_wider(names_from = RECYEAR,
              values_from = n) %>%
  column_to_rownames(var = "TRANSNUM") %>%
  as.matrix()

n.years <- length(unique(occ2$RECYEAR))



# Get covariates in order -------------------------------------------------

####Effort covariate - site, year, visit month
#Now we need to make an array of the observed
#effort data with rows for sites, columns for years,
# and matrix elements for each replicate

#now, generate IDs for the for loop where 
# we will populate the matrix
yr <- survey$yrID
site <- survey$TransID  
rep <- survey$REP

#make a blank array with dims of sites x years x reps
effort <- array(NA, dim = c(n.transects, #rows
                         n.years, #columns
                         4 #array elements - number of replicates
                         ))

#fill that array based on the values in those columns
# for occupancy
for(i in 1:dim(survey)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the site of row i,
  # the year of row i and the replicate of row i,
  # populate that space in the array with the column in
  # the dataframe that corresponds to the scaled vis data
  # for that sitexyearxreplicate combo
  effort[site[i], yr[i], rep[i]] <- as.numeric(survey[i,6])
}

#now, generate IDs for the for loop where 
# we will populate the matrix
yr <- occ2$yrID #get a yearID for each iteration of the loop
site <- occ2$TransID #site ID for each iteration fo the loop
spec <- occ2$SpecID #get a species ID for each iteration of the loop
rep <- occ2$REP #get a replicate for each iteration of the loop

y <- array(NA, dim = c(n.species, #rows
                       n.transects, #column
                       n.years, #first array level
                       4 #second array level
))

#fill taht array based on the values in those columns
# for occupancy
for(i in 1:dim(occ2)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the species of row i,
  #the site of row i, 
  # the year of row i and the replicate of row i,
  # populate that space in the array with the column in
  # the dataframe that corresponds to the 1-0 occupancy
  # for that speciesxyearxreplicate combo
  y[spec[i], site[i], yr[i], rep[i]] <- as.numeric(occ2[i,7])
}


#generate the z-matrix of speciesxsitexyear which assumes
# no false positives
#z matrix is speciesxsitexyear summed over all reps
zdf <- occ2 %>% 
  group_by(TransID, yrID, SpecID) %>%
  #get total occupancy for each speciesxyear
  summarise(tot_occ = sum(presence, na.rm = T)) %>%
  #set that to 1-0 values again
  mutate(tot_occ = case_when(tot_occ > 0 ~ 1,
                             tot_occ == 0 ~ 0,
                             TRUE ~ NA_real_)) 



nspecs <- max(zdf$SpecID) #get the dimension of rows
nsites <- max(zdf$TransID) #get dimension for columns
nyrs <- max(zdf$yrID) #get the dimension of 3rd dimension

#now, generate IDs for the for loop where 
# we will populate the matrix
yr <- zdf$yrID #get a yearID for each iteration of the loop
site <-zdf$TransID #site ID for each iteration fo the loop
spec <- zdf$SpecID #get a species ID for each iteration of the loop

#make a blank array with dims of species x years 
z <- array(NA, dim = c(nspecs, nsites, nyrs))

#fill taht array based on the values in those columns
# for occupancy
for(i in 1:dim(zdf)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the species of row i,
  #the site of row i, 
  # the year of row i and the replicate of row i,
  # populate that space in the array with the column in
  # the dataframe that corresponds to the 1-0 occupancy
  # for that speciesxyearxreplicate combo
  z[spec[i], site[i], yr[i]] <- as.numeric(zdf[i,4])
}

#set all zeros (which could be true or false) to NA
z[z == 0] <- NA


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

t <- occ2 %>%
  group_by(yrID, TransID, SpecID) %>%
  summarise(presence = mean(presence, na.rm = T)) %>%
  ungroup() %>%
  unite("site_year", c("yrID", "TransID"),
        sep = "_") %>%
  dplyr::select(SpecID, presence, site_year) %>%
  pivot_wider(names_from = SpecID,
              values_from = presence,
              values_fill = 0) %>%
  column_to_rownames(var = "site_year") %>%
  mutate(across(everything(), ~replace_na(.x, 0)))
# 

ggcorrplot(cor(t), type = "lower",
           lab = FALSE)

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
             y = y,
             z = z,
             R = R,
             omega.init = omega.init)

saveRDS(data, here('04_nps_plants',
                   'example',
                   'data_outputs',
                   'model_inputs',
                   'JAGS_data_list.RDS'))
