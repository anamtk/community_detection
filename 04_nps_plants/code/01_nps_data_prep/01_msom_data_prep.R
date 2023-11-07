# Shelby Lamm
# September 15, 2023
# NPS plant data prep 


# this script preps a multiple plots for a dynamic occupancy model


# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", "MASS",
                  'ggcorrplot')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

#read data into R
# plants <- read.csv(here('04_nps_plants',
#                         'data_raw',
#                         'NPS_veg_data_PEFO_S.csv'))
# 
# str(plants)

#this is in the new email from Megan, and includes the lifegroups and 
#durations for all the plants
plants <- read.csv(here('04_nps_plants',
                         'data_raw',
                         'NPS_veg_data.PEFO_S_v2.csv'))

str(plants)


# Remove unknown species --------------------------------------------------

plants %>%
  distinct(CurrentSpecies)

#these are the species we can remove from the dataset because we have
#no info (or little) on their growth forms and duration
#10052015-1
#PEFO_S08_20221025_unk1
#PEFO_S09_20221020_unk2
#PEFO_S10_20221022_unk1
#PEFO20101001_1
#PEFO20111012_1
#UnknownForb

plants <- plants %>%
  filter(!CurrentSpecies %in% c('10052015-1',
                                'PEFO_S08_20221025_unk1',
                                'PEFO_S09_20221020_unk2',
                                'PEFO_S10_20221022_unk1',
                                'PEFO20101001_1',
                                'PEFO20111012_1',
                                'UnknownForb'))

# want to subset the data by grabbing 10 of the 30 plots
# prioritizing plots with repeated measurements
subsample <- plants %>%
  distinct(Plot, Obs_type) %>%
  mutate(REP = case_when(Obs_type == "Regular" ~ 1,
                         TRUE ~ 2)) %>%
  filter(REP == 2) %>%
  # separate letter and number in plot column
  mutate(across(c('Plot'), substr, 2, nchar(Plot)))


subset_plots = as.numeric(sort(subsample$Plot))
# randomly choose 10 plots
#samp_plots = sort(sample(1:14, 10, replace=FALSE))
samp_plots = c(1,3,5,6,7,8,10,12,13,14)
subset_plots <- sort(subset_plots[c(samp_plots)])
subset_plots <- sub("^","S",subset_plots)

plants <- plants %>%
  filter(Plot %in% c(subset_plots))


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

# Ma

# including species:
occ2 <- plants %>%
  #factor species name so we can complete it by survey interval
  mutate(CurrentSpecies = as.factor(CurrentSpecies)) %>%
  #group by all the survey ID info
  group_by(EventYear, Plot, Transect, Quadrat) %>%
  #make sure each species is in each survey
  complete(CurrentSpecies) %>%
  ungroup() %>%
  group_by(CurrentSpecies) %>%
  fill(Lifeform, .direction = "updown") %>%
  fill(Duration, .direction = "updown") %>%
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


check <- occ2 %>%
  #filter(Obs_type == "Regular") %>%
  #distinct(yrID, quadID, SpecID, REP) %>%
  group_by(SpecID) %>%
  mutate(n = sum(presence))


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
  dplyr::select(quadID, yrID) %>%
  arrange(quadID, yrID) %>%
  group_by(quadID) %>%
  filter(yrID == max(yrID)) %>%
  ungroup() %>%
  dplyr::select(yrID) %>%
  as_vector()
n.yr <- unname(n.yr)

n.rep <- rep %>%
  group_by(quadID, yrID) %>%
  pivot_wider(names_from = yrID, 
              values_from = REP,
              values_fn = length) %>%
  column_to_rownames(var = "quadID") %>%
  as.matrix()

n.years <- length(unique(occ2$yrID))


# Get covariates in order -------------------------------------------------

#lifeform_duration, cover classes here as covariates

#LIFE FORMS
#some species identified with just a number or with just a genus,
#so these look like they often have an unknown duration
#I think we can just reclassify these as "unknown grass", "unknown forb"

lifeforms <- occ2 %>%
  unite(lifegroup,
        c(Lifeform, Duration),
        sep = "_",
        remove = F)

t <- lifeforms %>%
  distinct(CurrentSpecies, lifegroup)

#for cell-referencing, you will want to make the group with the 
#most species in it the "baseline"
lifeforms %>%
  distinct(CurrentSpecies, lifegroup) %>%
  group_by(lifegroup) %>%
  tally() %>%
  arrange(desc(n))
#looks like it is forb_annual

lifeforms2 <- lifeforms %>%
  distinct(SpecID, lifegroup) %>%
  #factor lifegroup to put in a logical order, with highest number first
  mutate(lifegroup2 = factor(lifegroup, levels = c('forb_annual', "forb_biennial",
                                                  'forb_perennial','forb_NA',
                                                  'graminoid_annual','graminoid_perennial',
                                                  'graminoid_NA',
                                                  'shrub_perennial', 
                                                  'cactus_perennial', 
                                                   'succulent_perennial'))) %>%
  #make this numeric because JAGS wants numeric factor variables
  mutate(lifegroup2 = as.numeric(lifegroup2))

#pull out as a vector for the model
lifegroup <- lifeforms2 %>%
  dplyr::select(lifegroup2) %>%
  as_vector()
lifegroup <- unname(lifegroup)

#get a number for indexing in the model
n.groups <- length(unique(lifeforms2$lifegroup2))

#COVER CLASSES
#SHELBY - I've set this up so that it should work for you, but I think we 
#should make this actually a scaled value based on the median of each cover
#class
#What that means is to basically just set each cover class to = it's median 
#continuous percentage value
#then use the scale function to get the scaled value of this
#code that should work for this (thoough you'll have to change all the numbers
#after the ~s in the case_when:
# covers <- occ2 %>%
#   mutate(cover = case_when(CoverClass == 0 ~ 0,
#                            CoverClass == 1 ~ 0.05,
#                            CoverClass == 2 ~ 0.3,
#                            CoverClass == 3 ~ 0.75,
#                            CoverClass == 4 ~ 1.5,
#                            CoverClass == 5 ~ 3.5,
#                            CoverClass == 6 ~ 7.5,
#                            CoverClass == 7 ~ 12.5,
#                            CoverClass == 8 ~ 20,
#                            CoverClass == 9 ~ 30,
#                            CoverClass == 10 ~ 42.5,
#                            CoverClass == 11 ~ 62.5,
#                            CoverClass == 12 ~ 87.5,
#                            #double check that the NAs came from "completing"
#                            #species in the pipe that created "occ2" above
##AMtK - did you change this because the NA's were actually 0s above?
## I changed back, but left this here just in case I misunderstood
#                            is.na(CoverClass) ~ 0,
#                            TRUE ~ 0)) %>% #NA_real_)) %>%
#   mutate(cover = scale(cover))

covers <- occ2 %>%
  mutate(cover = case_when(CoverClass == 0 ~ 0,
                           CoverClass == 1 ~ 0.05,
                           CoverClass == 2 ~ 0.3,
                           CoverClass == 3 ~ 0.75,
                           CoverClass == 4 ~ 1.5,
                           CoverClass == 5 ~ 3.5,
                           CoverClass == 6 ~ 7.5,
                           CoverClass == 7 ~ 12.5,
                           CoverClass == 8 ~ 20,
                           CoverClass == 9 ~ 30,
                           CoverClass == 10 ~ 42.5,
                           CoverClass == 11 ~ 62.5,
                           CoverClass == 12 ~ 87.5,
                           TRUE ~ 0)) %>%
                           #TRUE ~ NA_real_)) %>%
  mutate(cover = scale(cover))

yr <- covers$yrID #get a yearID for each iteration of the loop
site <- covers$quadID #site ID for each iteration fo the loop
spec <- covers$SpecID #get a species ID for each iteration of the loop
rep <- covers$REP #get a replicate for each iteration of the loop

cover <- array(NA, dim = c(n.species,
                           n.quads,
                           n.years,
                           2))

for(i in 1:dim(covers)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the species of row i,
  #the site of row i, 
  # the year of row i and the replicate of row i,
  # populate that space in the array with the column in
  # the dataframe that corresponds to the cover class group
  #for that combo
  #(NEED TO UPDATE - probably to median %cover for the cover class)
  cover[spec[i], site[i], yr[i], rep[i]] <- as.numeric(covers[i,18])
}

##AMtK - with imputing code in model - you shouldn't need to do this
#cover[is.na(cover)] <-0

# Get response data  ------------------------------------------------------



#now, generate IDs for the for loop where 
# we will populate the matrix
yr <- occ2$yrID #get a yearID for each iteration of the loop
site <- occ2$quadID #site ID for each iteration fo the loop
spec <- occ2$SpecID #get a species ID for each iteration of the loop
rep <- occ2$REP #get a replicate for each iteration of the loop

y <- array(NA, dim = c(n.species, #rows
                       n.quads, #column
                       n.years, #first array level
                       2))

#fill that array based on the values in those columns
# for [occupancy] presence
for(i in 1:dim(occ2)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the species of row i,
  #the site of row i, 
  # the year of row i and the replicate of row i,
  # populate that space in the array with the column in
  # the dataframe that corresponds to the 1-0 occupancy
  # for that speciesxyearxreplicate combo
  y[spec[i], site[i], yr[i], rep[i]] <- as.numeric(occ2[i,13])
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
#R<-diag(x=0.1, n.species, n.species)

#omega also needs priors, which I'm going to attempt to define using
#covariance among species abundances, we'll see how it goes

# t <- occ2 %>%
#   group_by(yrID, quadID, SpecID) %>%
#   summarise(presence = mean(presence, na.rm = T)) %>%
#   ungroup() %>%
#   unite("site_year", c("yrID", "quadID"),
#         sep = "_") %>%
#   dplyr::select(SpecID, presence, site_year) %>%
#   pivot_wider(names_from = SpecID,
#               values_from = presence,
#               values_fill = 0) %>%
#   column_to_rownames(var = "site_year") %>%
#   mutate(across(everything(), ~replace_na(.x, 0)))
# # 
# 
# ggcorrplot(cov(t), type = "lower",
#            lab = FALSE)
# 
# t[colSums(t, na.rm = T) == 0]
# t[rowSums(t, na.rm = T) == 0]
# 
# #set omega init to this - not sure if it will work with the NA values
# #or if i will need to define those as a value?? we can try it...
# t1 <- t[1:770,]
# t1[colSums(t1, na.rm = T) == 0]
# t2 <- t[771:1540,]
# t3 <- t[1541:2310,]
# 
# cov1 <- cov(t1, use = "complete.obs")
# cov1[cov1 == 0] <- 0.001
# cov2 <- cov(t2, use = "complete.obs")
# cov2[cov2 == 0] <- 0.001
# cov3 <- cov(t3, use = "complete.obs")
# cov3[cov3 == 0] <- 0.001
# 
# omega.init1 <- ginv(cov1)
# omega.init2 <- ginv(cov2)
# omega.init3 <- ginv(cov3)

# Prep list for JAGS ------------------------------------------------------

data <- list(n.species = n.species,
             n.quads = n.quads,
             n.yr = n.yr,
             n.rep = n.rep,
             cover = cover,
             lifegroup = lifegroup,
             n.groups = n.groups,
             y = y,
             z = z)

saveRDS(data, here('04_nps_plants',
                   'data_outputs',
                   'MSAM',
                   'model_inputs',
                   'nps_msam_multisite_subset.RDS'))

