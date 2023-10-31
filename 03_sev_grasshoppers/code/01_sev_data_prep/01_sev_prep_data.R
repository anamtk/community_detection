#Sevilleta grasshopper MSAM data prep
#Ana Miller-ter Kuile
#September 11, 2023

#this script preps the sevilleta grasshopper dataset for an MSAM

#Sites for Sevilleta grasshoppers, based on temporal coverage:

#NPP is from 1999 for two sites: 
#Core-Blak gramma and Core-Creosote
#these sites are also called "Five Points" in the datasets
#I think for meteorological station this includes station:
#49 (Five points) 

#NPP are seasonal - spring and fall of each year
#looks like grouping them by "web" and then taking an average of 
#the plots on that web for all species biomass in each season
#would be the way to combine with grasshopper data

#Climate are hourly - only one site near both the creosote
#and black gramma sites

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse",
                  'ggcorrplot', 'readxl')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load data ---------------------------------------------------------------

hopper <- read.csv(here('03_sev_grasshoppers',
                        'data_raw',
                        'sev106_grasshopper_counts.csv'))

life <- read.csv(here('03_sev_grasshoppers',
                      'data_raw',
                      'SEV_grasshopper_species_list2020.csv'))

# Subset data -------------------------------------------------------------

#only have NPP for 20 years from core_black and core_creosote, so going
#to just look at those two sites

hopper2 <- hopper %>%
  filter(ecosystem %in% c("core_black", "core_creosote"))

#47 species in these two systems
hopper2 %>%
  distinct(SPECIES) %>%
  tally()

# Get unique IDs for species, sites, years --------------------------------

#first, get complete species x survey dataset - set all unobserved to 0
#remove NONE species from all
hopper3 <- hopper2 %>%
  #get a unique site_web_transect ID column
  unite(site_web_trans, 
        c("SITE", "WEB", "TRN"),
        sep = "_",
        remove = F) %>%
  dplyr::select(YEAR, season, site_web_trans,
                SPECIES, CNT) %>%
  group_by(YEAR, season, site_web_trans, SPECIES) %>%
  summarise(CNT = sum(CNT, na.rm = T)) %>%
  ungroup() %>%
  group_by(YEAR, season, site_web_trans) %>%
  mutate(SPECIES = as.factor(SPECIES)) %>%
  complete(SPECIES) %>%
  filter(SPECIES != "NONE") %>%
  mutate(CNT = case_when(is.na(CNT) ~ 0,
                         TRUE ~ CNT)) %>%
  ungroup() %>%
  mutate(SPECIES = as.character(SPECIES)) %>%
  #get ID columns lined up
  mutate(yrID = as.numeric(as.factor(YEAR)),
         siteID = as.numeric(as.factor(site_web_trans)),
         speciesID = as.numeric(as.factor(SPECIES)),
         rep = as.numeric(as.factor(season)))

# Link grasshoppers to their reproductive strategy ------------------------
  
hopper4 <- hopper3 %>%
  left_join(life, by = c("SPECIES" = "CODE")) %>%
  dplyr::select(-FAMILY, -SUBFAMILY,
                -Lifeform, -MacrohabitatMicrohabitat,
                -Diet, -Comments) %>%
    #MESP2 missing from datatable, so i'll set to MESP
    mutate(LifeHistory = case_when(SPECIES == "MESP2" ~ "SU",
                                   TRUE ~ LifeHistory))

# Data objects for jags ---------------------------------------------------

n.species <- length(unique(hopper4$SPECIES))
n.transects <- length(unique(hopper4$site_web_trans))
n.years <- length(unique(hopper4$yrID))
n.start <- hopper4 %>%
  distinct(siteID, yrID) %>%
  group_by(siteID) %>%
  filter(yrID == min(yrID)) %>%
  arrange(siteID) %>%
  ungroup() %>%
  dplyr::select(yrID) %>%
  as_vector()

n.end <- hopper4 %>%
  distinct(siteID, yrID) %>%
  group_by(siteID) %>%
  filter(yrID == max(yrID)) %>%
  arrange(siteID) %>%
  ungroup() %>%
  dplyr::select(yrID) %>%
  as_vector()

#site x year matrix
n.rep <- hopper4 %>%
  distinct(siteID, yrID, rep) %>%
  group_by(siteID, yrID) %>%
  tally(name = "REP") %>%
  pivot_wider(names_from = yrID,
              values_from = REP) %>%
  column_to_rownames(var = "siteID") %>%
  as.matrix()

#no missing
n.rep[which(is.na(n.rep))] 

#which has the most?
hopper4 %>%
  distinct(SPECIES, LifeHistory) %>%
  group_by(LifeHistory) %>%
  tally()
#SU is the highest number

#for sevilleta
#maybe consider removing categorical variable - it likely evens out given
#that the data are collected in two seasons
reprod <- hopper4 %>%
  mutate(LifeHistory = factor(LifeHistory, levels = c("SU", "SP", "MV"))) %>%
  mutate(LifeHistory2 = as.numeric(LifeHistory)) %>%
  distinct(speciesID, LifeHistory2) %>%
  arrange(speciesID) %>%
  column_to_rownames(var = "speciesID") %>%
  as_vector()

# Prep observed data ------------------------------------------------------

#Now we need to make an array of the observed
#data with rows for species, columns for sites,
# third dimension for years, and fourth dimension for reps in 
#each year
nspecs <- max(hopper4$speciesID) #get the dimension of rows
nsites <- max(hopper4$siteID) #get dimension for columns
nyrs <- max(hopper4$yrID) #get the dimension of 3rd dimension
nreps <- max(hopper4$rep) #get the dimension of  4th dimension

#now, generate IDs for the for loop where 
# we will populate the matrix
yr <- hopper4$yrID #get a yearID for each iteration of the loop
site <- hopper4$siteID #site ID for each iteration fo the loop
spec <- hopper4$speciesID #get a species ID for each iteration of the loop
rep <- hopper4$rep #get a replicate for each iteration of the loop

#make a blank array with dims of species x years x reps
y <- array(NA, dim = c(nspecs, nsites, nyrs, nreps))

#fill taht array based on the values in those columns
# for occupancy
for(i in 1:dim(hopper4)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the species of row i,
  #the site of row i, 
  # the year of row i and the replicate of row i,
  # populate that space in the array with the column in
  # the dataframe that corresponds to the 1-0 occupancy
  # for that speciesxyearxreplicate combo
  y[spec[i], site[i], yr[i], rep[i]] <- as.numeric(hopper4[i,5])
}

#generate the total counts by site and year (regardless of replicate)
#for initial ymax value
Ndf <- hopper4 %>%
  group_by(siteID, yrID, speciesID) %>%
  #get total abundance for each sitexspeciesxyear
  summarise(tot = sum(CNT, na.rm = T)) 

nspecs <- max(Ndf$speciesID) #get the dimension of rows
nsites <- max(Ndf$siteID) #get dimension for columns
nyrs <- max(Ndf$yrID) #get the dimension of 3rd dimension

#now, generate IDs for the for loop where
# we will populate the matrix
yr <- Ndf$yrID #get a yearID for each iteration of the loop
site <-Ndf$siteID #site ID for each iteration fo the loop
spec <- Ndf$speciesID #get a species ID for each iteration of the loop

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

t <- hopper4 %>%
  group_by(yrID, siteID, speciesID) %>%
  summarise(COUNT = mean(CNT, na.rm = T)) %>%
  ungroup() %>%
  unite("site_year", c("yrID", "siteID"),
        sep = "_") %>%
  dplyr::select(speciesID, COUNT, site_year) %>%
  pivot_wider(names_from = speciesID,
              values_from = COUNT,
              values_fill = 0) %>%
  column_to_rownames(var = "site_year") %>%
  mutate(across(everything(), ~replace_na(.x, 0)))
# 

ggcorrplot(cor(t), type = "lower",
           lab = FALSE)

#set omega init to this - not sure if it will work with the NA values
#or if i will need to define those as a value?? we can try it...
omega.init <- cor(t)
mean(omega.init, na.rm = T)
omega.init[which(is.na(omega.init))] <- 0.04

# Data list ---------------------------------------------------------------


#make the data list
data <- list(y = y,
             reprod = reprod,
             n.species = n.species,
             n.years = n.years,
             n.start = n.start,
             n.end = n.end,
             n.transects = n.transects,
             n.rep = n.rep,
             #for initials
             ymax = ymax,
             omega.init = omega.init,
             #for omega prior
             R = R)

#export that for using with the model
saveRDS(data, here('sev_grasshopper',
                   "data_outputs",
                   "model_inputs",
                   "sev_msam_dynmultisite.RDS"))


# Extract info on sites and years -----------------------------------------

str(hopper4)

ids <- hopper4 %>%
  distinct(YEAR, site_web_trans, yrID, siteID) %>%
  separate(site_web_trans,
           into = c('site', 'web', 'transect'),
           remove = F)

write.csv(ids, here('03_sev_grasshoppers',
               'data_outputs',
               'metadata',
               'sev_site_year_IDs.csv'))

# Get bray for site one ---------------------------------------------------

#get bray for site 1 
# BOER_1_108 is SiteID 1, so I'll just use it for now

#Just for visualization purposes to compare "raw" vs "corrected" bray
#we will want to get a community matrix for one site across years

#we'll do site one, and do whta folks used to do and just take the 
#maximum number of individuals observed per species per year across
#repeat surveys

matrix <- hopper4 %>%
  filter(site_web_trans == "BOER_1_108") %>%
  group_by(YEAR, speciesID) %>%
  summarise(COUNT = max(CNT, na.rm = T)) %>%
  pivot_wider(names_from = YEAR,
              values_from = COUNT) %>%
  column_to_rownames(var = 'speciesID')

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
years <- 1992:2019

raw_bray <- as.data.frame(cbind(raw_bray = bray,
                                year = years))

saveRDS(raw_bray, here("05_visualizations",
                       "viz_data",
                       "sev_BOER_1_108_raw_bray.RDS"))


# Uncorrected bray for all sites-years ------------------------------------

#pulling out and calculating bray-curtis for all 
#communities so we can compare to the "corrected" version
#after the model

diss_fun <- function(site){
  
  matrix <- hopper4 %>%
    filter(siteID == site) %>%
    group_by(YEAR, speciesID) %>%
    summarise(COUNT = max(CNT, na.rm = T)) %>%
    pivot_wider(names_from = YEAR,
                values_from = COUNT) %>%
    column_to_rownames(var = 'speciesID')
  
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
  
  bray_df <- hopper4 %>%
    filter(siteID == site) %>%
    distinct(yrID, siteID) %>%
    cbind(bray) %>%
    mutate(type = "observed")
  
  return(bray_df)
}

sites <- unique(hopper4$siteID)
results <- lapply(sites, FUN = diss_fun)

results_df <- do.call(rbind, results)

saveRDS(results_df, here('05_visualizations',
                         'viz_data',
                         'sev_observed_bray.RDS'))


