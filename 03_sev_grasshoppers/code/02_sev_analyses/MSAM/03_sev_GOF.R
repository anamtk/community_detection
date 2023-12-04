#SEV MSAM GOF
#Ana Miller-ter Kuile
#November 16, 2023


#this script looks at observed vs. predicted y and calculates
#some total N stuff for the sev MSAM

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

observed <- read.csv(here('03_sev_grasshoppers',
                          'data_outputs',
                          'MSAM',
                          'sev_tidy_data_for_model.csv'))

modeled <- readRDS(here('03_sev_grasshoppers',
                        'monsoon',
                        'MSAM',
                        'outputs',
                        'sev_MSAM_GOF_summary.RDS'))


# Pull out yrep -----------------------------------------------------------

#y.rep is indexed species, transects, years, replicates

yrep <- as.data.frame(modeled$statistics) %>%
  rownames_to_column(var = 'parm') %>%
  filter(str_detect(parm, "y.rep")) %>%
  separate(parm,
           into = c("speciesID", "siteID", "yrID", "rep"),
           sep = ",") %>%
  mutate(speciesID = str_sub(speciesID, 7, nchar(speciesID))) %>%
  mutate(rep = str_sub(rep, 1, (nchar(rep)-1))) %>%
  mutate(speciesID = as.numeric(speciesID),
         siteID = as.numeric(siteID),
         yrID = as.numeric(yrID),
         rep = as.numeric(rep)) %>%
  left_join(observed, by = c("speciesID", "siteID", "yrID", "rep"))


# Plot --------------------------------------------------------------------

m1 <- lm(Mean ~ CNT,
         data = yrep)

summary(m1)

lb1 <- paste("R^2 == 0.45")

ggplot(yrep, aes(x = CNT, y = Mean)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point() +
  annotate(geom = 'text', x = 20, y = 2.5, label = lb1, parse = T) +
  labs(x = "observed", y = "predicted", title = "SEV LTER grasshopper MSAM GOF")

ggsave(plot = last_plot(),
       filename = here("pictures",
                       "detection_models",
                       "SEV_GOF_graph.jpg"),
       height = 4,
       width = 6,
       units = "in")





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

ggplot() +
  geom_density(data = obs, aes(x = total), fill = "black", alpha = 0.4) +
  geom_density(data = N_tot, aes(x = total), fill = 'blue', alpha = 0.4) +
  scale_x_log10()

ggplot() +
  geom_density(data = obs, aes(x = mean), fill = 'black', alpha = 0.4) +
  geom_density(data = N, aes(x = Mean), fill = "blue", alpha = 0.4) +
  scale_x_log10() +
  labs(x = "Mean detection")
