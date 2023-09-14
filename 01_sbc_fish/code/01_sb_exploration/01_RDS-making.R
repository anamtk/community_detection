
# An Bui
# May 11, 2023
# reading in raw data and creating RDS


# raw data ----------------------------------------------------------------

#fish survey data
fish <- read.csv(here("data_raw",
                      "Monthly_Fish_All_Years_20221018.csv"))

# fish from LTE
fish_LTE <- read.csv(here("data_raw",
                          "LTE_All_Species_Biomass_at_transect_20230323.csv"))


# RDS making --------------------------------------------------------------


# create a fish survey data RDS to push
saveRDS(fish, here("data_outputs",
                   "data_RDS",
                   "fish_monthly.RDS"))

saveRDS(fish_LTE, here("data_outputs",
                       "data_RDS",
                       "fish_LTE.RDS"))