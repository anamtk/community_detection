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


observed <- readRDS(here('01_sbc_fish',
                         "data_outputs",
                         'SAM',
                         "model_inputs",
                         "bray_SAM_input_data_long.RDS"))

modeled <- readRDS(here('01_sbc_fish',
                        'data_outputs',
                        'SAM',
                        'model_outputs',
                        'fish_SAM_GOF_summary.RDS'))


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

lb1 <- paste("R^2 == 0.19")

ggplot(bray.rep, aes(x = bray, y = bray.rep.Mean)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  geom_errorbar(aes(ymin = bray.rep.Mean - bray.rep.SD,
                    ymax = bray.rep.Mean + bray.rep.SD)) +
  annotate(geom = 'text', x = 0.8, y = 0.3, label = lb1, parse = T) +
  labs(x = 'observed', y = 'modeled', title = "SBC LTER fish SAM GOF")

ggsave(plot = last_plot(),
       filename = here("pictures",
                       "supplementary",
                       'SAM',
                       "SBC_SAM_GOF_graph.jpg"),
       height = 4,
       width = 6,
       units = "in")


# Check for autocorrelation -----------------------------------------------

auto_check <- as.data.frame(cbind(observed$Transect.ID, 
                                  observed$Year.ID)) %>%
  rename("Transect.ID" = "V1",
         "Year.ID" = "V2")

resid <- as.data.frame(modeled$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "resid")) %>%
  mutate(id = str_sub(parm, 7, (nchar(parm)-1))) %>%
  dplyr::select(id, `50%`)

auto_check <- auto_check %>%
  bind_cols(resid) %>%
  arrange(Transect.ID, Year.ID)

auto_function <- function(site){
  
  df <- auto_check %>%
    filter(Transect.ID == site)
  
  auto <- acf(df$`50%`, type = "correlation")

  return(auto)
}

sites <- unique(auto_check$Transect.ID)
site_list <- lapply(sites, auto_function)
#low autocorrelation for a few sites (~7/43)


