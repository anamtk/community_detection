#TUrnover model visuals
#Ana Miller-ter Kuile
#July 5, 2023

#results of the SAM models

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  'patchwork')


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

med <- read.csv(here("data_outputs",
                    "SAM_outputs",
                    "loss_SAM_summary.csv"))

medg <- read.csv(here("data_outputs",
                     "SAM_outputs",
                     "gain_SAM_summary.csv"))


# Loss model --------------------------------------------------------------

theme_set(theme_bw())
al <- med %>%
  filter(parameter %in% c("b[1]", "b[2]")) %>%
  mutate(beta = case_when(parameter == "b[1]" ~ "Kelp Biomass",
                          parameter == "b[2]" ~ "Temperature")) %>%
  ggplot(aes(x = beta, y = `50%`)) +
  geom_point() + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
  labs(x = "Variable", y = "Median and 95% BCI",
       title = "A. Variable coefficients") +
  geom_hline(yintercept = 0, linetype = 2) +
  coord_flip()

#higher loss with higher bottom temperatures

lagID <- c("WARM", "COLD", "WARM-1", "COLD-1", "WARM-2", "COLD-2")
bl <- med %>%
  filter(str_detect(parameter, "wB")) %>%
  mutate(lag = str_sub(parameter, start = 4, end = (nchar(parameter)-1))) %>%
  mutate(lag = as.numeric(lag)) %>%
  mutate(lag = lag - 1) %>%
  mutate(lag = as.factor(lag)) %>%
  ggplot(aes(x = lag, y = `50%`)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
  geom_hline(yintercept = 1/13, linetype = 2) +
  scale_x_discrete(labels = lagID) +
  labs(x = "Season", 
       y = "Importance weight \n (median and 95% BCI)",
       title = "B. Temperature lags")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

al + bl +
  plot_annotation('Loss model')


# Gain model --------------------------------------------------------------

ag <- medg %>%
  filter(parameter %in% c("b[1]", "b[2]")) %>%
  mutate(beta = case_when(parameter == "b[1]" ~ "Kelp Biomass",
                          parameter == "b[2]" ~ "Temperature")) %>%
  ggplot(aes(x = beta, y = `50%`)) +
  geom_point() + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
  labs(x = "Variable", y = "Median and 95% BCI",
       title = "A. Variable coefficients") +
  geom_hline(yintercept = 0, linetype = 2) +
  coord_flip()

#higher loss with higher bottom temperatures

bg <- medg %>%
  filter(str_detect(parameter, "wA")) %>%
  mutate(lag = str_sub(parameter, start = -2, end = -2)) %>%
  mutate(lag = as.numeric(lag)) %>%
  mutate(lag = lag - 1) %>%
  ggplot(aes(x = lag, y = `50%`)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
  geom_hline(yintercept = 1/6, linetype = 2) +
  labs(x = "Years in past", y = "Importance weight \n (median and 95% BCI)",
       title = "B. Kelp biomass lags")

cg <- medg %>%
  filter(str_detect(parameter, "wB")) %>%
  mutate(lag = str_sub(parameter, start = 4, end = (nchar(parameter)-1))) %>%
  mutate(lag = as.numeric(lag)) %>%
  mutate(lag = lag - 1) %>%
  mutate(lag = as.factor(lag)) %>%
  ggplot(aes(x = lag, y = `50%`)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
  geom_hline(yintercept = 1/13, linetype = 2) +
  scale_x_discrete(labels = lagID) +
  labs(x = "Season", 
       y = "Importance weight \n (median and 95% BCI)",
       title = "C. Temperature lags")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ag / (bg + cg) +
  plot_annotation('Gain model')
