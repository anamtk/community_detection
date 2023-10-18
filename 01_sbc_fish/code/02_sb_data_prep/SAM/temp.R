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

bottemp <- read.csv(here('01_sbc_fish',
                         "data_raw",
                         "environmental",
                         "Bottom_temp_all_years_20230724.csv"))


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

write.csv(bottemp3, here('01_sbc_fish',
                         "data_outputs",
                         'SAM',
                         "data_prep",
                         "monthly_bottom_temps.csv"))


bottemp4 <- bottemp2 %>%
  mutate(SEASON = case_when(MONTH %in% c("12", "01", "02",
                                         "03", "04", "05") ~ "COLD",
                            MONTH %in% c("06", "07", "08", "09",
                                         "10", "11") ~ "WARM")) %>%
  group_by(SITE, YEAR, SEASON) %>%
  summarise(sd_TEMP = sd(TEMP_C, na.rm = T),
            TEMP_C = mean(TEMP_C, na.rm = T)) %>%
  ungroup()

write.csv(bottemp4, here('01_sbc_fish',
                         "data_outputs",
                         'SAM',
                         'data_prep',
                         "seasonal_bottom_temps.csv"))
