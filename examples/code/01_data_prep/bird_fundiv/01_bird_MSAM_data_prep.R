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

birds <- read.csv(here('examples',
                       'data_raw',
                       'bird_fundiv',
                       'CBP011.csv'))
                       


str(birds)
#for body size
#Supplementary Data 1 here: https://figshare.com/s/b990722d72a26b5bfead
avonet3 <- read.csv(here('examples',
                         'data_raw',
                         'bird_fundiv',
                         'AVONET_eBird.csv'))

#to link the konza birds via codes and scientific names to the avonet data
#from here: https://www.birdpop.org/pages/birdSpeciesCodes.php
codes <- read.csv(here('examples',
                       'data_raw',
                       'bird_fundiv',
                       'IBP-AOS-LIST23.csv'))

# Filter out watersheds of interest ---------------------------------------

birds %>%
  distinct(WATERSHED)
#OG: only going to do watersheds where we also have NPP data
#which are 001D, 004B, and 020B

#going to do all watersheds for now and can subset these other 
#ones with only NPP later - maybe do a SAM with and a SAM without NPP data

#update: now going to only do watersheds in the burning+grazing experiment
#4 yr burn, grazed: N04D, N04B
#4 yr burn, ungrazed: 004A, 004B
#1 year burn, grazed: N01B (two transects)
#1 year burn, ungrazed: 001D, 001A/R20A
#no burn, grazed: N20B
#no burn, ungrazed: 020B, 020C, 020D (can't find on map - removing)

#update update:
#going to pick one transect in each treatment that has higher detections,
#and or primary productivity datae??
#code below to figure this out:
#N04D, 004B, N01B-10, 001D, N20B, 020B

birds1 <- birds %>%
  filter(WATERSHED %in% c("N04D", "N04B",
                          "004A", "004B",
                          "N01B",
                          "001D", "R20A",
                          "N20B",
                          "020B", "020C")) %>%
  #filter(TRANSNUM != 6) %>%
  filter(COMMONNAME !=  "Transect not run" )

# birds1 <- birds %>% 
#   filter(WATERSHED %in% c('001D', '004B','020B'))%>%
#   filter(COMMONNAME !=  "Transect not run" )


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
#sparrow sp.
#Spotted/Eastern Towhee
#Western/Eastern Meadowlark
#unique ID is actually common name, not SPECNAME

#clean up so that things will mesh well
birds2 <- birds1 %>%
  mutate(COMMONNAME = case_when(COMMONNAME == "Le Conte's Sparrow" ~ "LeConte's Sparrow",
                                COMMONNAME == "Greater Prairie-Chicken" ~ "Greater Prairie-chicken",
                              TRUE ~ COMMONNAME)) %>%
  mutate(AOUCODE = case_when(AOUCODE == "eawp" ~ "EAWP",
                             AOUCODE == "CAGO" ~ "CANG",
                             TRUE ~ AOUCODE)) %>%
  filter(COMMONNAME != "Empidonax sp.") %>%
  filter(!str_detect(COMMONNAME, "/"))  %>%
  filter(COMMONNAME != "sparrow sp.")

#Get unique Bird IDS with common name, AOU code and scientific names
bird_id2 <- birds2 %>%
  distinct(COMMONNAME, SPECNAME, AOUCODE) %>%
  left_join(codes1, by = c("COMMONNAME")) %>%
  #American Goldfinch is acting up and I don't know why....!!! GAH
  mutate(SPEC = case_when(COMMONNAME == "American Goldfinch" ~ "AGOL",
                          COMMONNAME == "Barred Owl" ~ "BADO",
                          COMMONNAME == "Greater Prairie-chicken" ~ 'GRPC',
                          COMMONNAME == "Black-and-White Warbler" ~ "BAWW",
                          COMMONNAME == "Canada Goose" ~ "CANG",
                          TRUE ~ SPEC)) %>%
  mutate(SCINAME = case_when(COMMONNAME == "American Goldfinch" ~ "Spinus tristis",
                             COMMONNAME == "Barred Owl" ~ 'Strix varia',
                             COMMONNAME == "Greater Prairie-chicken" ~ 'Tympanuchus cupido',
                             COMMONNAME == "Black-and-White Warbler" ~ "Mniotilta varia",
                             COMMONNAME == "Canada Goose" ~ "Branta canadensis",
                             TRUE ~ SCINAME)) %>%
  filter(SPECNAME != "NONE") %>%
  #and remove birds of prey, shorebirds, and gamebirds as a test
  filter(!AOUCODE %in% c(#birds of prey and shorebirds
    "AMKE", "BEKI", "CANG",
    "CONI", "COHA", "GOEA",
    "GBHE", "GHOW", "GRHE",
    "KILL", "NOHA", "RTHA",
    "RLHA", "SEOW", "SWHA",
    "UPSA", "WODU", "TUVU",
    #gamebirds
    "COPO", "GRPC", "NOBO", 
    "RNEP", "WITU",
    #woodpeckers
    "DOWO", "HAWO", "RBWO",
    "RHWO")) %>%
  distinct(COMMONNAME, SPEC, SCINAME) 
#78 species


# Manupulate data to abundance structure ----------------------------------

#Data are individuals seen at different distances from a line transect

#need to summarise data as count data for each recording
#period in each year for each species

birds3 <- birds2 %>%
  #and remove birds of prey, shorebirds, and gamebirds as a test
  filter(!AOUCODE %in% c(#birds of prey and shorebirds
    "AMKE", "BEKI", "CANG",
    "CONI", "COHA", "GOEA",
    "GBHE", "GHOW", "GRHE",
    "KILL", "NOHA", "RTHA",
    "RLHA", "SEOW", "SWHA",
    "UPSA", "WODU", "TUVU",
    #gamebirds
    "COPO", "GRPC", "NOBO", 
    "RNEP", "WITU",
    #woodpeckers
    "DOWO", "HAWO", "RBWO",
    "RHWO")) %>%
  group_by(RECYEAR, RECMONTH, RECDAY, TRANSNUM,
           WATERSHED, COMMONNAME, AOUCODE) %>%
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
  group_by(COMMONNAME) %>%
  fill(AOUCODE, .direction = "updown") %>%
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
#transect 18 started in 1982


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

hist(sizes$Mass)
#get info from bird dataset on ids of different surveys to blend to
#survey dataset
survey_meta <- birds4 %>%
  distinct(RECYEAR, RECMONTH, RECDAY, TRANSNUM, WATERSHED,
           TransID, yrID)

#survey length
survey <- birds %>%
  #filter(WATERSHED %in% c("001D", "004B", "020B")) %>%
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
  y[spec[i], site[i], yr[i], rep[i]] <- as.numeric(birds5[i,8])
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



# Random effects variables ------------------------------------------------

n.sites <- n.transects

Site.ID <- birds5 %>%
  distinct(TransID) %>%
  as_vector()


Year.ID <- years$yearnum

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
             n.sites = n.sites,
             Year.ID = Year.ID,
             Site.ID = Site.ID,
             #initials
             ymax = ymax)


#export that for using with the model
saveRDS(data, here('examples',
                   "data_output",
                   'bird_fundiv',
                   "model_inputs",
                   "bird_msam_data_input_list.RDS"))

# Export metadata ---------------------------------------------------------

write.csv(birds5, here('examples',
                       "data_output",
                       'bird_fundiv',
                       'tidy_data',
                       'bird_msam_tidy_data.csv'))

sites <- birds5 %>%
  distinct(RECYEAR, TRANSNUM,
           WATERSHED,
           TransID, yrID)

write.csv(sites, here('examples',
                      "data_output",
                      'bird_fundiv',
                      'tidy_data',
               'bird_msam_site_year_IDs.csv'))


# Summaries ---------------------------------------------------------------

birds5 %>%
  distinct(RECYEAR, RECMONTH, TRANSNUM, WATERSHED, DURATION, REP) %>%
  summarise(dur = mean(DURATION, na.rm = T),
            min = min(DURATION, na.rm = T),
            max = max(DURATION, na.rm = T))

sizes %>%
  summarise(mean = mean(Mass),
            min = min(Mass),
            max = max(Mass))


# Transect selection based on number of birds observed --------------------

birds5 %>%
  ungroup() %>%
  group_by(RECYEAR, TRANSNUM, WATERSHED) %>%
  summarise(mean = mean(NOBS, na.rm =T)) %>%
  mutate(TRANSNUM = as.factor(TRANSNUM)) %>%
  ggplot(aes(x = TRANSNUM, y = mean)) +
  geom_boxplot()

birds5 %>%
  distinct(WATERSHED, TRANSNUM) %>%
  arrange(TRANSNUM)

#from this - N04D, 004B, N01B-10, 001D, N20B, 020B
  
