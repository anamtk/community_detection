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


observed <- readRDS(here('02_konza_birds',
                         "data_outputs",
                         'SAM',
                         "model_inputs",
                         "knz_bray_SAM_input_data.RDS"))

modeled <- readRDS(here('02_konza_birds',
                        'data_outputs',
                        'SAM',
                        'model_outputs',
                        'knz_SAM_GOF_summary.RDS'))


# Pull out y and yrep -----------------------------------------------------

bray <- observed$bray

bray.rep <- as.data.frame(modeled$statistics) %>%
  rownames_to_column(var = 'parm') %>%
  filter(str_detect(parm, "bray.rep")) %>%
  rename(bray.rep.Mean = Mean,
         bray.rep.SD = SD) %>%
  cbind(bray)

m1 <- lm(bray.rep.Mean ~ bray,
         data = bray.rep)

summary(m1)

lb1 <- ("R^2 == 0.38")
ggplot(bray.rep, aes(x = bray, y = bray.rep.Mean)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  geom_errorbar(aes(ymin = bray.rep.Mean - bray.rep.SD,
                    ymax = bray.rep.Mean + bray.rep.SD)) +
  annotate(geom = "text", x = 0.46, y = 0.68, label = lb1, parse = T) +
  labs(x = "observed", y = "modeled", title = "KNZ LTER bird SAM GOF")

ggsave(plot = last_plot(),
       filename = here("pictures",
                       "supplementary",
                       'SAM',
                       'KNZ_SAM_GOF_graph.jpg'),
       height = 4,
       width = 6,
       units = "in")

