
# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  "jagsUI",
                  'rjags',
                  'mcmcplots',
                  "coda",
                  'patchwork') #mcmc output


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

data_list <- readRDS(here('02_konza_birds',
                          "data_outputs",
                          'SAM',
                          "model_inputs",
                          "knz_bray_SAM_input_data.RDS"))


# Parameters to save ------------------------------------------------------

params <- c('b',
            'b0',
            'wA',
            'wB',
            'wC',
            'var.process')



# JAGS model --------------------------------------------------------------

model <- here('02_konza_birds',
              "code", 
              "02_kb_analyses",
              'SAM',
              "jags",
              "bird_SAM.R")

Sys.time()
mod <- jagsUI::jags(data = data_list,
                    inits = NULL,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.burnin = 1000,
                    n.iter = 8000,
                    DIC = TRUE)

Sys.time()

# Check convergence -------------------------------------------------------

mcmcplot(mod$samples)

gelman.diag(mod$samples, multivariate = F)


# Look at some plots ------------------------------------------------------

sum <- summary(mod$samples)

med <- as.data.frame(sum$quantiles) %>%
  rownames_to_column(var = "parameter") %>%
  dplyr::select(parameter, `2.5%`, `50%`, `97.5%`) %>%
  filter(parameter != "deviance")
# 
# write.csv(med, here("01_sbc_fish",
#                     "data_outputs",
#                     "SAM",
#                     "model_outputs",
#                     "fish_SAM_summary.csv"))


# Quick plots -------------------------------------------------------------

med2 <- med %>%
  filter(parameter %in% c("b[1]", 'b[2]', 'b[3]'))

ggplot(med2, aes(x = parameter, y= `50%`)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
  labs(y = "Median (and 95% BCI)", x = "Covariate") +
  scale_x_discrete(labels = c("Temperature", "Precipitation", "NPP")) +
  coord_flip()

weights <- med %>%
  filter(str_detect(parameter, "wB|wA")) %>%
  mutate(type = case_when(str_detect(parameter, "wA") ~ "wA",
                          str_detect(parameter, "wB") ~ "wB")) %>%
  mutate(parameter = factor(parameter, levels = c("wA[1]", "wA[2]", "wA[3]",
                                                  "wA[4]", "wA[5]", "wA[6]",
                                                  "wA[7]", "wA[8]", "wA[9]",
                                                  "wA[10]", "wA[11]", "wA[12]",
                                                  "wB[1]", "wB[2]", "wB[3]",
                                                  "wB[4]", "wB[5]", "wB[6]",
                                                  "wB[7]", "wB[8]", "wB[9]",
                                                  "wB[10]", "wB[11]", "wB[12]")))

1/12
ggplot(weights, aes(x = parameter, y = `50%`)) +
  geom_hline(yintercept = 1/12, linetype = 2) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
  theme(axis.text.x = element_text(angle = 45, hjust =1)) +
  facet_wrap(~type, scales = "free")
