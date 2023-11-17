#Grasshopper SAM GOF
#Ana Miller-ter Kuile
#November 16, 2023


#this script looks at observed vs. predicted bray for the
#grasshopper SAM model

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


observed <- readRDS(here('03_sev_grasshoppers',
                         "data_outputs",
                         'SAM',
                         "model_inputs",
                         "sev_bray_SAM_input_data.RDS"))

modeled <- readRDS(here('03_sev_grasshoppers',
                        "monsoon",
                        'SAM',
                        "outputs",
                        "hopper_SAM_GOF_summary.RDS"))


# Pull out y and yrep -----------------------------------------------------

bray <- observed$bray

bray.rep <- as.data.frame(modeled$statistics) %>%
  rownames_to_column(var = 'parm') %>%
  filter(str_detect(parm, "bray.rep")) %>%
  rename(bray.rep.Mean = Mean,
         bray.rep.SD = SD) %>%
  cbind(bray)

ggplot(bray.rep, aes(x = bray, y = bray.rep.Mean)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  geom_errorbar(aes(ymin = bray.rep.Mean - bray.rep.SD,
                    ymax = bray.rep.Mean + bray.rep.SD))

m1 <- lm(bray.rep.Mean ~ bray,
         data = bray.rep)

summary(m1)
