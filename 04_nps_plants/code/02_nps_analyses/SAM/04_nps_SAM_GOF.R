#Fish SAM GOF
#Ana Miller-ter Kuile
#November 16, 2023


#this script looks at observed vs. predicted bray for the
#fish SAM model

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


observed <- readRDS(here('04_nps_plants',
                   "data_outputs",
                   'SAM',
                   "model_inputs",
                   "nps_diss_SAM_input_data.RDS"))

modeled <- readRDS(here('04_nps_plants',
                        'monsoon',
                        'nps_SAM',
                        'outputs',
                        'nps_SAM_GOF_summary.RDS'))


# Pull out y and yrep -----------------------------------------------------

diss <- observed$diss

diss.rep <- as.data.frame(modeled$statistics) %>%
  rownames_to_column(var = 'parm') %>%
  filter(str_detect(parm, "diss.rep")) %>%
  rename(diss.rep.Mean = Mean,
         diss.rep.SD = SD) %>%
  cbind(diss)

m1 <- lm(diss.rep.Mean ~ diss,
         data = diss.rep)

summary(m1)

lb1 <- paste("R^2 == 0.77")

ggplot(diss.rep, aes(x = diss, y = diss.rep.Mean)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  geom_errorbar(aes(ymin = diss.rep.Mean - diss.rep.SD,
                    ymax = diss.rep.Mean + diss.rep.SD)) +
  annotate(geom = 'text', x = 0.8, y = 0.3, label = lb1, parse = T) +
  labs(x = 'observed', y = 'modeled', title = "PFNP plant SAM GOF")

ggsave(plot = last_plot(),
       filename = here("pictures",
                       "supplementary",
                       'SAM',
                       "PFNP_SAM_GOF_graph.jpg"),
       height = 4,
       width = 6,
       units = "in")

