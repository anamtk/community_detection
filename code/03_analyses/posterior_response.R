# look at posteriors for MSOM
# Ana Miller-ter Kuile
# January 23, 2023

# this script looks at the preliminary results of the MSOM

# for benthic communities at the SBC LTER

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 

package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

## And loading them

for(i in package.list){library(i, character.only = T)}

# Load posterior data -----------------------------------------------------

summary_df <- readRDS(here("data_outputs",
                           "monsoon_outputs",
                           "benthic_MSOM_1_20_summary.RDS"))

species_ID <- read.csv(here("data_outputs",
                            "raw_community",
                            "single_site_MSOM_matrix.csv"))

# Look at distributions of posteriors -------------------------------------

df_sum <- as.data.frame(summary_df) %>%
rownames_to_column(var = "parameter") 

# Plot visibility by species effect ---------------------------------------

#get species ID for the names of the variables

species <- colnames(species_ID)


species <- species[2:length(species)]

df_vis <- df_sum %>%
filter(str_detect(parameter, "a1.Vis")) %>%
cbind(species)

ggplot(df_vis, aes(x = reorder(species, `50%`), y = `50%`)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`)) +
  coord_flip()

df_psi <- df_sum %>%
  filter(str_detect(parameter, "psi")) %>%
  filter(parameter != "sd.lpsi") %>%
  cbind(species)


ggplot(df_psi, aes(x = reorder(species, `50%`), y = `50%`)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`)) +
  coord_flip() +
  labs(y = "median occupancy probability",
       x = "species")

# Relationship between visibility effect and occupancy probability? -------

psi_vis_df <- df_vis %>%
left_join(df_psi, by = "species")

ggplot(psi_vis_df, aes(x = `50%.x`, y = `50%.y`)) +
  geom_point() +
  labs(x = "median visibility effect",
       y = "median occupancy probability")

#doesnt seem to be actually
m <- lm(`50%.x` ~ `50%.y`,
        data = psi_vis_df)

summary(m)