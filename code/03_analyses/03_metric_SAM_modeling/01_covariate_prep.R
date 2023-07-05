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
  geom_point(position = position_jitter(width = 0.25))

#standardize:
#whatever months have + are one season, whatever months
#have - values are the other season

bottemp3 %>%
  mutate(TEMP_C = scale(TEMP_C)) %>%
  group_by(MONTH) %>%
  summarise(mean = mean(TEMP_C))
  
bottemp3 %>%
  group_by(MONTH) %>%
  summarise(mean = mean(TEMP_C))

write.csv(bottemp3, here("data_outputs",
                         "community_stability",
                         "monthly_bottom_temps.csv"))


#looks like splitting up months <15C on average
#and >=15C might be a good call?

bottemp4 <- bottemp2 %>%
  mutate(SEASON = case_when(MONTH %in% c("12", "01", "02",
                                         "03", "04", "05") ~ "COLD",
                            MONTH %in% c("06", "07", "08", "09",
                                         "10", "11") ~ "WARM")) %>%
  group_by(SITE, YEAR, SEASON) %>%
  summarise(sd_TEMP = sd(TEMP_C, na.rm = T),
            TEMP_C = mean(TEMP_C, na.rm = T)) %>%
  ungroup()

ggplot(bottemp4, aes(x = SEASON, y = TEMP_C)) +
  geom_point(position = position_jitter(width = 0.25))

bottemp5 <- bottemp2 %>%
  group_by(SITE, YEAR) %>%
  summarise(sd_TEMP = sd(TEMP_C, na.rm = T),
            TEMP_C = mean(TEMP_C, na.rm = T)) %>%
  ungroup()

ggplot(bottemp5, aes(x = YEAR, y = TEMP_C, color = SITE)) +
  geom_point() +
  geom_line(aes(group = SITE))

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

bottemp6 <- bottemp4 %>%
  filter(SITE %in% site)

bottemp6 %>%
  group_by(SEASON) %>%
  summarise(mean = mean(TEMP_C))

write.csv(bottemp6, here("data_outputs",
                         "community_stability",
                         "seasonal_bottom_temps.csv"))

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

#temperatuer lags - seasonally - going back 
#how many seasons now??? HMMM....

temp_lags <- bottemp6 %>%
  group_by(SITE) %>%
  arrange(SITE, YEAR, SEASON) %>%
  #this creates a column for a lag 1:12 months ago
  do(data.frame(., setNames(shift(.$TEMP_C, 1:5), c("TEMP_C_l1",
                                                     "TEMP_C_l2", "TEMP_C_l3",
                                                     "TEMP_C_l4", "TEMP_C_l5")))) %>%
  ungroup() %>%
  dplyr::select(SITE, YEAR, SEASON,
                TEMP_C:TEMP_C_l5)

#set first season of temp lags to be 
#"WARM" since this is when all the surveys took place
temp_lags2 <- temp_lags %>%
  filter(SEASON == "WARM") %>%
  dplyr::select(-SEASON) %>%
  mutate(YEAR = as.integer(YEAR))

# Combine all data --------------------------------------------------------

all_data <- stability %>%
  left_join(bio_lags, by = c("YEAR",
                             "SITE_TRANS")) %>%
  left_join(temp_lags2, by = c("SITE", "YEAR"))
  

write.csv(all_data, here("data_outputs",
                        "community_stability",
                        "stability_metrics_with_covariates.csv"))


