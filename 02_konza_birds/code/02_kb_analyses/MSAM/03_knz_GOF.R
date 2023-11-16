#Bird MSAM GOF
#Ana Miller-ter Kuile
#November 16, 2023


#this script looks at observed vs. predicted y and calculates
#some total N stuff for the bird MSAM

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

observed <- read.csv(here('02_konza_birds',
                          'data_outputs',
                          'MSAM',
                          'knz_tidy_data_for_model.csv'))

modeled <- readRDS(here('02_konza_birds',
                        'monsoon',
                        'MSAM',
                        'outputs',
                        'bird_MSAM_GOF_summary.RDS'))


# Pull out yrep -----------------------------------------------------------

#y.rep is indexed species, transects, years, replicates

yrep <- as.data.frame(modeled$statistics) %>%
  rownames_to_column(var = 'parm') %>%
  filter(str_detect(parm, "y.rep")) %>%
  separate(parm,
           into = c("SpecID", "TransID", "yrID", "REP"),
           sep = ",") %>%
  mutate(SpecID = str_sub(SpecID, 7, nchar(SpecID))) %>%
  mutate(REP = str_sub(REP, 1, (nchar(REP)-1))) %>%
  mutate(SpecID = as.numeric(SpecID),
         TransID = as.numeric(TransID),
         yrID = as.numeric(yrID),
         REP = as.numeric(REP)) %>%
  left_join(observed, by = c("SpecID", "TransID", "yrID", "REP"))


# Plot --------------------------------------------------------------------

ggplot(yrep, aes(x = NOBS, y = Mean)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point()

m1 <- lm(Mean ~ NOBS,
         data = yrep)

summary(m1)



# Pull out N --------------------------------------------------------------

N <- as.data.frame(modeled$statistics) %>%
  rownames_to_column(var = 'parm') %>%
  filter(str_detect(parm, "N"))%>%
  separate(parm,
           into = c("SpecID", "TransID", "yrID"),
           sep = ",") %>%
  mutate(SpecID = str_sub(SpecID, 3, nchar(SpecID))) %>%
  mutate(yrID = str_sub(yrID, 1, (nchar(yrID)-1))) %>%
  mutate(SpecID = as.numeric(SpecID),
         TransID = as.numeric(TransID),
         yrID = as.numeric(yrID)) 

N_tot <- N %>%
  group_by(SpecID) %>%
  summarise(total = sum(Mean, na.rm = T)) %>%
  ungroup()

obs <- observed %>%
  group_by(AOUCODE) %>%
  summarise(total = sum(NOBS),
            mean = mean(NOBS))

ggplot(obs, aes(x = total)) +
  geom_histogram()

