# Shelby Lamm
# September 15, 2023
# NPS plant data prep 


# this script preps a multiple plots for a dynamic occupancy model


# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse",
                  'ggcorrplot')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

#read data into R
plants <- read.csv(here('04_nps_plants',
                        'data_raw',
                        'NPS_veg_data_PEFO_S.csv'))

str(plants)


# Manipulate data structure ----------------------------------

# ignoring species:
years <- plants %>%
  distinct(EventYear, Plot, Transect, Quadrat) %>%
  group_by(Plot, Transect, Quadrat) %>%
  mutate(quadID = cur_group_id()) %>%
  arrange(EventYear) %>%
  mutate(yrID = 1:n()) %>%
  ungroup() %>%
  distinct(quadID, yrID, EventYear)


survey <- plants %>%
  distinct(EventYear, Plot, Transect, Quadrat, Obs_type) %>%
  group_by(Plot, Transect, Quadrat) %>%
  mutate(quadID = cur_group_id()) %>%
  ungroup() %>%
  mutate(REP = case_when(Obs_type == "Regular" ~ 1,
                         TRUE ~ 2)) %>%
  left_join(years, by = c("quadID", "EventYear"))


rep <- survey %>%
  distinct(yrID, quadID, REP)



# including species:
occ2 <- plants %>%
  #factor species name so we can complete it by survey interval
  mutate(CurrentSpecies = as.factor(CurrentSpecies)) %>%
  #group by all the survey ID info
  group_by(EventYear, Plot, Transect, Quadrat) %>%
  #make sure each species is in each survey
  complete(CurrentSpecies) %>%
  ungroup() %>%
  #presence: 1 when cover class does not equal 0, 0 if cover class = 0
  mutate(presence = case_when(CoverClass == 0 ~ 0,
                              TRUE ~ 1)) %>%
  mutate(CurrentSpecies = as.factor(as.character(CurrentSpecies))) %>%
  #mutate(yrID = as.numeric(as.factor(EventYear))) %>%
  mutate(SpecID = as.numeric(as.factor(CurrentSpecies))) %>%
  group_by(Plot, Transect, Quadrat) %>%
  mutate(quadID = cur_group_id()) %>%
  ungroup() %>%
  mutate(REP = case_when(Obs_type == "Regular" ~ 1,
                         TRUE ~ 2)) %>%
  left_join(years, by = c("quadID", "EventYear"))



# see if missing years in the time series
check_missing <- occ2 %>%
  filter(Obs_type == "Regular") %>%
  distinct(yrID, quadID) %>%
  group_by(quadID) %>%
  tally()


# Prep data structure for JAGS --------------------------------------------

n.species <- length(unique(occ2$CurrentSpecies))

n.quads <- length(unique(occ2$quadID))

n.yr <- years %>%
  select(c(quadID, yrID)) %>%
  arrange(quadID, yrID) %>%
  group_by(quadID) %>%
  filter(yrID == max(yrID)) %>%
  ungroup() %>%
  dplyr::select(yrID) %>%
  as_vector()

n.rep <- rep %>%
  group_by(quadID, yrID) %>%
  pivot_wider(names_from = yrID, 
              values_from = REP,
              values_fn = length) %>%
  column_to_rownames(var = "quadID") %>%
  as.matrix()

n.years <- length(unique(occ2$yrID))


# Get covariates in order -------------------------------------------------


#now, generate IDs for the for loop where 
# we will populate the matrix
yr <- occ2$yrID #get a yearID for each iteration of the loop
site <- occ2$quadID #site ID for each iteration fo the loop
spec <- occ2$SpecID #get a species ID for each iteration of the loop
rep <- occ2$REP #get a replicate for each iteration of the loop

y <- array(NA, dim = c(n.species, #rows
                       n.quads, #column
                       n.years #first array level
))

#fill that array based on the values in those columns
# for [occupancy] presence
for(i in 1:dim(occ2)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the species of row i,
  #the site of row i, 
  # the year of row i and the replicate of row i,
  # populate that space in the array with the column in
  # the dataframe that corresponds to the 1-0 occupancy
  # for that speciesxyearxreplicate combo
  y[spec[i], site[i], yr[i]] <- as.numeric(occ2[i,11])
}


#generate the z-matrix of speciesxsitexyear which assumes
# no false positives
#z matrix is speciesxsitexyear summed over all reps
zdf <- occ2 %>% 
  group_by(quadID, yrID, SpecID) %>%
  #get total occupancy for each speciesxyear
  summarise(tot_occ = sum(presence, na.rm = T)) %>%
  #set that to 1-0 values again
  mutate(tot_occ = case_when(tot_occ > 0 ~ 1,
                             tot_occ == 0 ~ 0,
                             TRUE ~ NA_real_)) 


nspecs <- max(zdf$SpecID) #get the dimension of rows
nsites <- max(zdf$quadID) #get dimension for columns
nyrs <- max(zdf$yrID) #get the dimension of 3rd dimension

#now, generate IDs for the for loop where 
# we will populate the matrix
yr <- zdf$yrID #get a yearID for each iteration of the loop
site <-zdf$quadID #site ID for each iteration fo the loop
spec <- zdf$SpecID #get a species ID for each iteration of the loop

#make a blank array with dims of species x years 
z <- array(NA, dim = c(nspecs, nsites, nyrs))

#fill that array based on the values in those columns
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
  group_by(yrID, quadID, SpecID) %>%
  summarise(presence = mean(presence, na.rm = T)) %>%
  ungroup() %>%
  unite("site_year", c("yrID", "quadID"),
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
             n.quads = n.quads,
             n.yr = n.yr,
             n.rep = n.rep,
             y = y,
             z = z,
             R = R,
             omega.init = omega.init)

saveRDS(data, here('04_nps_plants',
                   'data_outputs',
                   'MSAM',
                   'model_inputs',
                   'nps_msam_dynmultisite.RDS'))

