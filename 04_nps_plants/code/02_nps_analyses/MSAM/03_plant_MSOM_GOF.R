#Plant MSOM GOF
#Ana Miller-ter Kuile
#November 16, 2023


#this script looks at observed vs. predicted y and calculates
#some z stuff for plants

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 

package.list <- c("here", "tidyverse",
                  "data.table")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

## And loading them

for(i in package.list){library(i, character.only = T)}

theme_set(theme_bw())


# Load data ---------------------------------------------------------------

observed <- read.csv(here('04_nps_plants',
                          'data_outputs',
                          'MSAM',
                          'pfnp_tidy_data_for_model.csv'))

modeled <- readRDS(here('04_nps_plants',
                        'monsoon',
                        'nps_MSAM',
                        'outputs_yrsite',
                        'plant_MSOM_GOF_summary.RDS'))



# Prep observed -----------------------------------------------------------

obs2 <- observed %>%
  dplyr::select(EventYear, Plot, Transect, Quadrat, quadnum,
                yrID, REP, SpecID, presence)
  
# Pull out yrep -----------------------------------------------------------

#y.rep is indexed species, years, quadrats, reps

yrep <- as.data.frame(modeled$statistics) %>%
  rownames_to_column(var = 'parm') %>%
  filter(str_detect(parm, "y.rep")) %>%
  separate(parm,
           into = c("SpecID", 'yrID', "quadnum", "REP"),
           sep = ",") %>%
  mutate(SpecID = str_sub(SpecID, 7, nchar(SpecID))) %>%
  mutate(REP = str_sub(REP, 1, (nchar(REP)-1))) %>%
  mutate(SpecID = as.numeric(SpecID),
         quadnum = as.numeric(quadnum),
         yrID = as.numeric(yrID),
         REP = as.numeric(REP)) %>%
  left_join(obs2, by = c("SpecID", "quadnum", "yrID", "REP"))


# Plot --------------------------------------------------------------------

ggplot(yrep, aes(x = presence, y = Mean)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point()

yrep %>%
  mutate(presence = as.factor(presence)) %>%
  ggplot(aes(x = presence, y = Mean)) +
  geom_boxplot()

m1 <- lm(Mean ~ presence,
         data = yrep)

summary(m1)



# Pull out z --------------------------------------------------------------

z <- as.data.frame(modeled$statistics) %>%
  rownames_to_column(var = 'parm') %>%
  filter(str_detect(parm, "z"))%>%
  separate(parm,
           into = c("SpecID", "yrID", "quadnum"),
           sep = ",") %>%
  mutate(SpecID = str_sub(SpecID, 3, nchar(SpecID))) %>%
  mutate(quadnum = str_sub(quadnum, 1, (nchar(quadnum)-1))) %>%
  mutate(SpecID = as.numeric(SpecID),
         quadnum = as.numeric(quadnum),
         yrID = as.numeric(yrID)) 

z_tot <- z %>%
  group_by(SpecID) %>%
  summarise(total = sum(Mean, na.rm = T)) %>%
  ungroup()

obs <- observed %>%
  group_by(SpecID) %>%
  summarise(total = sum(presence),
            mean = mean(presence))

ggplot() +
  geom_density(data = obs, aes(x = total), fill = "black", alpha = 0.4) +
  geom_density(data = z_tot, aes(x = total), fill = 'blue', alpha = 0.4)

