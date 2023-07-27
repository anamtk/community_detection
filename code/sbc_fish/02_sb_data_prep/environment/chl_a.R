#Extracting monthly Chlrophyll-A data
#Ana Miller-ter Kuile
#July 6, 2023
#based on code by Kristen Michaud
#https://doi.org/10.5281/zenodo.7117750


# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse",
                  'ncdf4', 'raster',
                  'rgdal','janitor')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

#functions to pull in all the rasters of interest
source(here("code",
            "00_functions",
            "tidy_functions.R"))

# Load data ---------------------------------------------------------------

#coordinates for each site to extract ChlA from
#doing site-level since the resolution is larger
#than the area that encompasses all the transects
#at a site
coords <- read.csv(here("data_raw", 
                        "environmental", 
                        "Benthic_Transect_Depths_coors.csv"))

#list of all the rasters in the folder with
#the rasters in it
folder <- here('data_raw', "environmental", "chl_a")

all_files <- list.files(folder, pattern = "\\.nc$", full.names = TRUE)


# Prep site-level coordinates ---------------------------------------------

coordinates_a <- coords %>% 
  janitor::clean_names() %>% 
  dplyr::select(site, latitude, longitude) %>% 
  group_by(site) %>% 
  summarise(lat = mean(latitude), lon = mean(longitude))

# Modify the reef locations by 2km offshore to reduce
#likelihood of bin overlap with land 
#(Chl files are in 4km bins)
coordinates <- coordinates_a %>% 
  mutate(lat_off = case_when(
    site == "ABUR" ~ (lat - 0.018),
    site == "AHND" ~ (lat - 0.018),
    site == "AQUE" ~ (lat - 0.018),
    site == "BULL" ~ (lat - 0.018),
    site == "CARP" ~ (lat - 0.018),
    site == "GOLB" ~ (lat - 0.018),
    site == "IVEE" ~ (lat - 0.018),
    site == "MOHK" ~ (lat - 0.018),
    site == "NAPL" ~ (lat - 0.018),
    site == "SCDI" ~ (lat + 0.018),
    site == "SCTW" ~ (lat + 0.018)))

lat.pts <- coordinates$lat_off
lon.pts <- coordinates$lon
extract.pts <- cbind(lon.pts, lat.pts)


# Set ROI -----------------------------------------------------------------

#region of interest to clip rasters to
ROI <- extent(-121,-119, 33.9, 34.6)
ROI


# Clip rasters ------------------------------------------------------------

chl1 <- extract_chl(all_files)

chl2 <- lapply(chl1, function(x) chla_df(x))

chl3 <- dplyr::bind_rows(chl2, .id = "date") %>%
  mutate(YEAR = str_sub(date, start = 1, end = 4)) %>%
  mutate(MONTH = str_sub(date, start = 5, end = 6)) %>%
  dplyr::select(-date)

ggplot(chl3, aes(x = YEAR, y = chla)) +
  geom_violin()

ggplot(chl3, aes(x = MONTH, y = chla)) +
  geom_point() +
  facet_wrap(~YEAR)

write.csv(chl3, here("data_outputs",
                     "community_stability",
                     "monthly_chla.csv"))
