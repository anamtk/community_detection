# Running the MSAM for fish
# Ana Miller-ter Kuile
# July 25, 2023

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


# Load Data ---------------------------------------------------------------

data <- readRDS(here("data_outputs",
                     "model_inputs",
                     "fish_msam_dynmultisite.RDS"))

# Parameters to save ------------------------------------------------------

params <- c(
            #COMMUNITY parameters
            'a1.Vis',
            'a2.Size',
            'lambda.mean',
            'sig.lambda',
            'a0.mean',
            'sig.a0')


#we found ymax to set initials, since otherwise the model will hate us
inits <- function() list(N = data$ymax)

# JAGS model --------------------------------------------------------------

model <- here("code", 
              "03_analyses",
              '01b_abundance_modeling',
              "models",
              "MSAM_multisite.R")

Sys.time()
mod <- jagsUI::jags(data = data,
                         inits = inits,
                         model.file = model,
                         parameters.to.save = params,
                         parallel = TRUE,
                         n.chains = 3,
                         n.iter = 1,
                         DIC = TRUE)

Sys.time()

# Check convergence -------------------------------------------------------

mcmcplot(mod$samples)

gelman.diag(mod$samples, multivariate = F)

# Raftery -----------------------------------------------------------------

raf <- raftery.diag(mod$samples)

names <- rownames(raf[[1]]$resmatrix)
ch1 <- raf[[1]]$resmatrix[,2]
ch2 <- raf[[2]]$resmatrix[,2]
ch3 <- raf[[3]]$resmatrix[,2]

raf_all <- as.data.frame(cbind(names, 
                               ch1, ch2, ch3)) %>%
  mutate(ch1 = as.numeric(ch1),
         ch2 = as.numeric(ch2),
         ch3 = as.numeric(ch3)) %>%
  filter(!str_detect(names, "z")) %>%
  pivot_longer(ch1:ch3,
               names_to = "chain",
               values_to = 'iterations') 

ggplot(raf_all, aes(x = iterations/3)) +
  geom_histogram() 

raf_all %>%
  summarise(iterations_90 = quantile(iterations, 
                                     probs = 0.9, 
                                     na.rm = T)/3,
            iterations_95 = quantile(iterations,
                                     probs = 0.95,
                                     na.rm = T)/3,
            max = max(iterations, 
                      na.rm = T)/3)
# A tibble: 1 × 3
# iterations_90 iterations_95   max
# <dbl>         <dbl> <dbl>
#   1        20717.        29211. 86112

bu1 <- raf[[1]]$resmatrix[,1]
bu2 <- raf[[2]]$resmatrix[,1]
bu3 <- raf[[3]]$resmatrix[,1]

burn <- as.data.frame(cbind(names, bu1, bu2, bu3)) %>%
  mutate(bu1 = as.numeric(bu1),
         bu2 = as.numeric(bu2),
         bu3 = as.numeric(bu3)) %>%
  filter(!str_detect(names, "z")) %>%
  pivot_longer(bu1:bu3,
               names_to = "chain",
               values_to = 'iterations') 

burn %>%
  summarise(max(iterations, na.rm = T))
#792


# Update model to track turnover components -------------------------------

theme_set(theme_bw())

parms2 <- c("tot_turnover",
            "gain",
            "loss",
            "jaccard")

mod2 <- update(mod,
               n.iter = 335,
               parameters.to.save = parms2)

community_fun <- function(mod, variable){
  
  q50 <- mod$q50
  q2.5 <- mod$q2.5
  q97.5 <- mod$q97.5
  
  med <- as.data.frame(get(variable, q50))  %>%
    rownames_to_column(var = "site") %>%
    pivot_longer(2:22,
                 names_to = "year",
                 values_to = "median")

  low <- as.data.frame(get(variable, q2.5))  %>%
    rownames_to_column(var = "site") %>%
    pivot_longer(2:22,
                 names_to = "year",
                 values_to = "lci")

  high <- as.data.frame(get(variable, q97.5))  %>%
    rownames_to_column(var = "site") %>%
    pivot_longer(2:22,
                 names_to = "year",
                 values_to = "uci")

  df <- med %>%
    left_join(low, by = c("site", "year")) %>%
    left_join(high, by = c("site", "year")) %>%
    mutate(year = str_sub(year, start = 2, end = length(year))) %>%
    mutate(year = as.numeric(year))

  plot <- ggplot(df, aes(x = year, y = median, color = site)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = lci, ymax = uci))
  
  return(plot)
  
}

turnp <- community_fun(mod = mod2, variable = 'tot_turnover')
turnp <- turnp + 
  labs(x = "Year",
       y = "Total turnover")

gainp <- community_fun(mod = mod2, variable = "gain")
gainp <- gainp + 
  labs(x = "Year",
       y = "Gains")
  
lossp <- community_fun(mod = mod2, variable = "loss")
lossp <- lossp + 
  labs(x = "Year",
       y = "Losses")

jaccp <- community_fun(mod = mod2, variable = "jaccard")
jaccp <- jaccp+ 
  labs(x = "Year",
       y = "Beta diversity")

library(patchwork)

turnp | (gainp / lossp) 

turnp + jaccp


turnover_samps <- mod2$sims.list$tot_turnover
gain_samps <- mod2$sims.list$gain
loss_samps <- mod2$sims.list$loss
beta_samps <- mod2$sims.list$jaccard

n.iter <- dim(turnover_samps)[1]
n.year <- dim(turnover_samps)[3]

iterations <- as.data.frame(1:n.iter) %>%
  rename("iteration" = "1:n.iter")
years <- as.data.frame(1:n.year) %>%
  rename('year' = '1:n.year')

sites <- as.data.frame(1:5) %>%
  rename(site = '1:5')

df <- iterations %>%
  cross_join(years) %>%
  cross_join(sites)

iter <- df$iteration
yr <- df$year
st <- df$site

df$turnover <- rep(NA, nrow(df))

for(i in 1:nrow(df)){
  df$turnover[i] <- turnover_samps[iter[i], st[i], yr[i]]
}

df$gain <- rep(NA, nrow(df))

for(i in 1:nrow(df)){
  df$gain[i] <- gain_samps[iter[i], st[i], yr[i]]
}

df$loss <- rep(NA, nrow(df))

for(i in 1:nrow(df)){
  df$loss[i] <- loss_samps[iter[i], st[i], yr[i]]
}

df$jaccard <- rep(NA, nrow(df))

for(i in 1:nrow(df)){
  df$jaccard[i] <- beta_samps[iter[i], st[i], yr[i]]
}

df %>%
  mutate(year = as.factor(year),
         site = as.factor(site)) %>%
  ggplot(aes(x = year, y= turnover, fill = site)) +
  geom_violin() +
  facet_wrap(~site)
