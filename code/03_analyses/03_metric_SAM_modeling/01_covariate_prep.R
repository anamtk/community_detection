
#Prepping environmental data for stability SAM
#Ana Miller-ter Kuile
#June 26, 2023

#this is a script that can prep environmental
#data for the SAM model

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 

package.list <- c("here", "tidyverse",
                  "data.table")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

## And loading them

for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

#giant kelp biomass (yearly)

#bottom temperature (every 15 minutes - maybe 
#summarise monthly for now? and later think about
#more biologically-relevant "seasons"?)

biomass <- read.csv(here("data_raw",
                         "environmental",
                         "Annual_All_Species_Biomass_at_transect_20230201.csv"))

bottemp <- read.csv(here("data_raw",
                         "environmental",
                         "Bottom_temp_all_years_20220729.csv"))

#to get the sites and transects we need for the other 
#two datasets
stability <- read.csv(here("data_outputs",
                           "community_stability",
                           "corrected_stability_metrics.csv"))
# Biomass by Year by Site -------------------------------------------------

colnames(biomass)
#i'll filter only giant kelp "MAPY"
#and the dry biomass value "DRY_GM2"

biomass <- biomass %>%
  dplyr::select(YEAR, MONTH, SITE, 
                TRANSECT, SP_CODE, DRY_GM2) %>%
  filter(SP_CODE == "MAPY")


# Summarise bottom temps --------------------------------------------------

colnames(bottemp)

bottemp2 <- bottemp %>%
  separate(DATE_LOCAL, into = c("YEAR", "MONTH", "DAY"),
           sep = "-")

#summarise by site, year, and month, and get SD in temp
#maybe want to do that by year - we'll see...
#or add in the variability as a lag later? who knows
bottemp3 <- bottemp2 %>%
  group_by(SITE, YEAR, MONTH) %>%
  summarise(sd_TEMP = sd(TEMP_C, na.rm = T),
            TEMP_C = mean(TEMP_C, na.rm = T)) %>%
  ungroup()
  
ggplot(bottemp3, aes(x = MONTH, y = TEMP_C)) +
  geom_point()

bottemp3 %>%
  group_by(MONTH) %>%
  summarise(mean = mean(TEMP_C))
#looks like splitting up months <15C on average
#and >=15C might be a good call?
# Get site and transect IDs -----------------------------------------------

sites <- stability %>%
  distinct(SITE_TRANS, siteID) %>%
  separate(SITE_TRANS, into= c("SITE", "TRANSECT"),
           sep = "_",
           remove = F) %>%
  mutate(siteID2 = as.numeric(as.factor(SITE)))

site <- sites$SITE
sitetrans <- sites$SITE_TRANS

# Filter environmental to sites and transects in dataset ------------------

bottemp4 <- bottemp3 %>%
  filter(SITE %in% site)

bottemp4 %>%
  group_by(MONTH) %>%
  summarise(mean = mean(TEMP_C))

biomass2 <- biomass %>%
  unite(c(SITE, TRANSECT),
        col = "SITE_TRANS",
        sep = "_",
        remove = F) %>%
  filter(SITE_TRANS %in% sitetrans) %>%
  mutate(DRY_GM2 = case_when(DRY_GM2 == -99999 ~ NA_real_,
                             TRUE ~ DRY_GM2))
  
# Make Lags ---------------------------------------------------------------

#biomass2:
#yearly lags - maybe start with 5-6 years back?
#site, transect, site_Transect, year
#Bio_1, lag back 5 years

bio_lags <- biomass2 %>%
  group_by(SITE_TRANS) %>%
  arrange(SITE_TRANS, YEAR) %>%
  #this creates a column for every lag 1:6 years ago
  do(data.frame(., setNames(shift(.$DRY_GM2, 1:5), c("DRY_GM2_l1",
                                                 "DRY_GM2_l2", "DRY_GM2_l3",
                                                 "DRY_GM2_l4", "DRY_GM2_l5")))) %>%
  ungroup() %>%
  dplyr::select(YEAR, SITE_TRANS, SITE,
                TRANSECT, DRY_GM2:DRY_GM2_l5)

#temperatuer lags - monthly
#going back 12 months right now
#site, year, month, lags

temp_lags <- bottemp4 %>%
  group_by(SITE) %>%
  arrange(SITE, YEAR, MONTH) %>%
  #this creates a column for a lag 1:12 months ago
  do(data.frame(., setNames(shift(.$TEMP_C, 1:18), c("TEMP_C_l1",
                                                     "TEMP_C_l2", "TEMP_C_l3",
                                                     "TEMP_C_l4", "TEMP_C_l5",
                                                     "TEMP_C_l6", "TEMP_C_l7",
                                                     "TEMP_C_l8", "TEMP_C_l9",
                                                     "TEMP_C_l10", "TEMP_C_l11",
                                                     "TEMP_C_l12", "TEMP_C_l13",
                                                     "TEMP_C_l14", "TEMP_C_l15",
                                                     "TEMP_C_l16", "TEMP_C_l17",
                                                     "TEMP_C_l18")))) %>%
  ungroup() %>%
  dplyr::select(SITE, YEAR, MONTH,
                TEMP_C:TEMP_C_l18)

#set first month of temp lags to be 07 - since
#this is the most common first month of survey
#for all the sites - and to make it consistent
#might make these more "seasonal" at some point
temp_lags2 <- temp_lags %>%
  filter(MONTH == '07') %>%
  dplyr::select(-MONTH) %>%
  mutate(YEAR = as.integer(YEAR))

# Combine all data --------------------------------------------------------

all_data <- stability %>%
  left_join(bio_lags, by = c("YEAR",
                             "SITE_TRANS")) %>%
  left_join(temp_lags2, by = c("SITE", "YEAR"))
  

# Prep data for jags ------------------------------------------------------

n.data <- nrow(all_data)

turn <- as.vector(all_data$tot_turnover)

n.transects <- length(unique(all_data$SITE_TRANS))

Transect.ID <- all_data$siteID

n.years <- length(unique(all_data$YEAR))

Year.ID <- all_data$yearID

n.sites <- length(unique(all_data$SITE))

Site.ID <- all_data %>%
  distinct(siteID, SITE) %>%
  arrange(siteID) %>%
  mutate(SITE = as.numeric(as.factor(SITE))) %>%
  dplyr::select(SITE) %>%
  as_vector()

n.kelplag <- all_data %>%
  dplyr::select(DRY_GM2:DRY_GM2_l5) %>%
  ncol()

Kelp <- all_data %>%
  dplyr::select(SITE_TRANS, YEAR, DRY_GM2:DRY_GM2_l5) %>%
  pivot_longer(DRY_GM2:DRY_GM2_l5,
               names_to = 'lag',
               values_to = 'kelp') %>%
  mutate(kelp = scale(kelp)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "kelp") %>%
  dplyr::select(DRY_GM2:DRY_GM2_l5) %>%
  as.matrix()

sum(is.na(Kelp))
sum(!is.na(Kelp))

#~10% missing data

n.templag <- all_data %>%
  dplyr::select(TEMP_C:TEMP_C_l18) %>%
  ncol()

Temp <- all_data %>%
  dplyr::select(SITE_TRANS, YEAR, TEMP_C:TEMP_C_l18) %>%
  pivot_longer(TEMP_C:TEMP_C_l18,
               names_to = 'lag',
               values_to = 'temp') %>%
  mutate(temp = scale(temp)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "temp") %>%
  dplyr::select(TEMP_C:TEMP_C_l18) %>%
  as.matrix()

sum(is.na(Temp))/(sum(is.na(Temp)) + sum(!is.na(Temp)))
#~12% missing data

data <- list(n.data = n.data,
             n.transects = n.transects,
             n.sites = n.sites,
             n.years = n.years,
             Transect.ID = Transect.ID,
             Year.ID = Year.ID,
             Site.ID = Site.ID,
             turn = turn,
             n.kelplag = n.kelplag,
             Kelp = Kelp,
             n.templag = n.templag,
             Temp = Temp)

saveRDS(data, here("data_outputs",
                   "model_inputs",
                   "turnover_SAM",
                   "turnover_SAM_input_data.RDS"))


write.csv(all_data, here("data_outputs",
                        "community_stability",
                        "stability_metrics_with_covariates.csv"))


