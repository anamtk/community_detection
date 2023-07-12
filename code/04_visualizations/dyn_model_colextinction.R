#Visuals of the dynamic model persistence and colonization
#Ana Miller-ter Kuile
#July 5, 2023

#persistence and colonizaiton posteriors

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 

package.list <- c("here", "tidyverse", "patchwork")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

## And loading them

for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------


summaries <- readRDS(here("monsoon",
                          "MSOM",
                          "outputs",
                          "fish_specieslevel_colext_summaries.RDS"))

ids <- read.csv(here("data_outputs",
                     "metadata",
                     "species_IDs.csv"))


# Pull out medians --------------------------------------------------------

meds <- as.data.frame(summaries$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(parm != "deviance") %>%
  separate(parm,
           into = c("parm", "species_year"),
           sep = "\\[") %>%
  separate(species_year,
           into = c("speciesID", "yearID"),
           sep = ",") %>%
  mutate(yearID = str_sub(yearID, 1, end = -2L)) %>%
  rename('median' = '50%',
         "lci" = '2.5%',
         'uci' = '97.5%') %>%
  dplyr::select(parm, speciesID, yearID, median, lci, uci)  %>%
  mutate(speciesID = as.numeric(speciesID)) %>%
  left_join(ids, by = c("speciesID" = "specID"))

meds %>%
  filter(parm == "phi") %>%
  summarise(mean_phi = mean(median))

a <- meds %>%
  filter(parm == "phi") %>%
  group_by(speciesID, SP_CODE) %>%
  summarise(mean_phi = mean(median),
            mean_lci = mean(lci),
            mean_uci = mean(uci)) %>%
  ggplot(aes(x = SP_CODE, y = mean_phi)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_lci, ymax = mean_uci),
                width = 0.2) +
  geom_hline(yintercept = 0.8141435) +
  labs(x = "Species", y = "Persistance probability") +
  coord_flip() 


meds %>%
  filter(parm == "gamma") %>%
  summarise(mean_phi = mean(median))

b <- meds %>%
  filter(parm == "gamma") %>%
  group_by(speciesID, SP_CODE) %>%
  summarise(mean_gamma = mean(median),
            mean_lci = mean(lci),
            mean_uci = mean(uci)) %>%
  ggplot(aes(x = SP_CODE, y = mean_gamma)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_lci, ymax = mean_uci),
                width = 0.2) +
  geom_hline(yintercept = 0.0434824) +
  labs(x = "Species", y = "Colonization probability") +
  coord_flip()

a + b


# Relationships bw probs? -------------------------------------------------

phi <- meds %>%
  filter(parm == "phi") %>%
  rename(median_phi = median) %>%
  dplyr::select(speciesID, SP_CODE, yearID, median_phi)

gamma <- meds %>%
  filter(parm == "gamma") %>%
  rename(median_gamma = median) %>%
  dplyr::select(speciesID, SP_CODE, yearID, median_gamma)

both <- phi %>%
  left_join(gamma, by = c("speciesID", "SP_CODE", "yearID"))

ggplot(both, aes(x = median_phi, y = median_gamma)) + 
  geom_point() +
  labs(x = "Persistence", y = "Colonization") 

# Another visual ----------------------------------------------------------

m1 <- meds %>%
  filter(parm == "phi") %>%
  ggplot(aes(x = SP_CODE, y = median)) +
  geom_hline(yintercept = 0.8141435) +
  geom_violin() +
  coord_flip() 

m2 <- meds %>%
  filter(parm == "gamma") %>%
  ggplot(aes(x = SP_CODE, y = median)) +
  geom_hline(yintercept = 0.0434824) +
  geom_violin() +
  coord_flip() 

m1 + m2

fishiD <- readRDS(here("data_outputs",
                       "metadata",
                       "species_metadata.RDS"))  

fishiD <- fishiD %>%
  distinct(sp_code, common_name, scientific_name)

meds %>%
  filter(parm == "phi") %>%
  group_by(speciesID, SP_CODE) %>%
  filter(all(median >= 0.8141435)) %>%
  ungroup() %>%
  left_join(fishiD, by = c("SP_CODE" = "sp_code")) %>%
  ggplot(aes(x = common_name, y = median)) +
  geom_violin() +
  labs(x = "Fish species", y = "Median yearly persistence",
       title = "Fish with above-average persistence") +
  coord_flip()

meds %>%
  filter(parm == "gamma") %>%
  group_by(speciesID, SP_CODE) %>%
  filter(all(median <= 0.0434824)) %>%
  ungroup() %>%
  left_join(fishiD, by = c("SP_CODE" = "sp_code")) %>%
  ggplot(aes(x = common_name, y = median)) +
  geom_violin() +
  labs(x = "Fish species", y = "Median yearly colonization",
       title = "Fish with below-average colonization") +
  coord_flip()



