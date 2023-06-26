
# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  "jagsUI",
                  'rjags',
                  'mcmcplots',
                  "coda") #mcmc output


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

data_list <- readRDS(here("data_outputs",
                          "model_inputs",
                          "turnover_SAM",
                          "turnover_SAM_input_data.RDS"))


# Parameters to save ------------------------------------------------------

params <- c("b0",
            'b',
            'wA',
            'wB')

# JAGS model --------------------------------------------------------------

model <- here("code", 
              "03_analyses",
              '03_metric_SAM_modeling',
              "jags",
              "turnover_SAM.R")

Sys.time()
mod <- jagsUI::jags(data = data_list,
                    inits = NULL,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 4000,
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

theme_set(theme_bw())
med %>%
  filter(parameter %in% c("b[1]", "b[2]", "b[3]")) %>%
  mutate(beta = case_when(parameter == "b[1]" ~ "Year",
                          parameter == "b[2]" ~ "Kelp Biomass",
                          parameter == "b[3]" ~ "Temperature")) %>%
  ggplot(aes(x = beta, y = `50%`)) +
  geom_point() + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
  labs(x = "Variable", y = "Median and 95% BCI") +
  geom_hline(yintercept = 0, linetype = 2) +
  coord_flip()

#higher turnover with increased kelp biomass; 
#higher turnover with higher bottom temperatures

med %>%
  filter(str_detect(parameter, "wA")) %>%
  mutate(lag = str_sub(parameter, start = -2, end = -2)) %>%
  mutate(lag = as.numeric(lag)) %>%
  mutate(lag = lag - 1) %>%
  ggplot(aes(x = lag, y = `50%`)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
  geom_hline(yintercept = 1/6, linetype = 2) +
  labs(x = "Years in past", y = "Importance weight \n (median and 95% BCI)",
       title = "Kelp biomass lags")

med %>%
  filter(str_detect(parameter, "wB")) %>%
  mutate(lag = str_sub(parameter, start = 4, end = (nchar(parameter)-1))) %>%
  mutate(lag = as.numeric(lag)) %>%
  mutate(lag = lag - 1) %>%
  mutate(lag = as.factor(lag)) %>%
  ggplot(aes(x = lag, y = `50%`)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
  geom_hline(yintercept = 1/13, linetype = 2) +
  labs(x = "Month", 
       y = "Importance weight \n (median and 95% BCI)",
       title = "Temperature lags") 

1/6
1/13
