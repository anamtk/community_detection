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

summary_list <- readRDS(here("data_outputs",
                             "monsoon_outputs",
                             "fish_MSOM_1_26_stats.RDS"))

# Look at distributions of posteriors -------------------------------------

df_sum <- as.data.frame(summary_list$summary) %>%
  rownames_to_column(var = "parameter") 


# Plot species-level parameters -------------------------------------------

df_vis <- df_sum %>%
filter(str_detect(parameter, "a1.Vis")) 

ggplot(df_vis, aes(x = reorder(parameter, `50%`), y = `50%`)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`)) +
  coord_flip()

df_size <- df_sum %>%
  filter(str_detect(parameter, "a2.Size"))

ggplot(df_size, aes(x = reorder(parameter, `50%`), y = `50%`)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`)) +
  coord_flip()


# Plot community-wide effects ---------------------------------------------

df_com <- as.data.frame(summary_list$community_sum) %>%
  rownames_to_column(var = "parameter")

df_com %>%
  filter(parameter %in% c("mu.vis", "mu.size")) %>%
  ggplot(aes(x = parameter, y = `50%`)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`)) 

df_psi <- as.data.frame(summary_list$psi_sum) %>%
  rownames_to_column(var = "parameter") %>%
  filter(parameter != "deviance")


ggplot(df_psi, aes(x = reorder(parameter, `50%`), y = `50%`)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`)) +
  coord_flip() +
  labs(y = "median occupancy probability",
       x = "species")


# Get z matrices out ------------------------------------------------------

zs <- as.data.frame(summary_list$z_sum) %>%
  rownames_to_column(var = "parameter") %>%
  filter(parameter != "deviance")
