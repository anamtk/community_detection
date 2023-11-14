
# Load packages ---------------------------------------------------------------
Sys.time()

# Load packages,
#package.list <- c('jagsUI',"coda",'dplyr', 
#                  'stringr','magrittr',
#                  'tidyr','ggplot2',
#                  'tibble', 'reshape2',
#                 'data.table') 

library(jagsUI)
library(coda)
library(dplyr)
library(purrr)
library(stringr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(tibble)
library(reshape2)
library(data.table)
## Installing them if they aren't already on the computer
#new.packages <- package.list[!(package.list %in% 
#                                 installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)

## And loading them
#for(i in package.list){library(i, character.only = T)}



# Load model --------------------------------------------------------------

model <- readRDS(file ="/scratch/atm234/nps_plants/outputs/nps_JAGS_model.RDS")

print('model loaded')
# Detection summary -------------------------------------------------------

mod_sum <- summary(model$samples) 

saveRDS(mod_sum, file = "/scratch/atm234/nps_plants/outputs/nps_detection_summary.RDS")

print('model summarized')

# update for p0 -----------------------------------------------------------

parms <- c('p0')

p0_update <- update(model,
                    parameters.to.save = parms,
                    n.iter = 1335)

print("detection updated")

p0_sum <- summary(p0_update$samples)

print("detection summarized")

saveRDS(p0_sum, file = "/scratch/atm234/nps_plants/outputs/nps_p0_summary.RDS")

print("detection saved")

# update for z ------------------------------------------------------------

parms2 <- c("z")

z_update <- update(model, 
                   parameters.to.save = parms2,
                   n.iter = 1335)

# Jaccard calculations ----------------------------------------------------

quadnums <- read.csv(file = "/scratch/atm234/nps_plants/inputs/site_year_IDs.csv")

samples <- z_update$sims.list[[1]]

test <- reshape2::melt(samples) %>%
  rename("iter" = "Var1",
         'SpecID' = "Var2",
         'yrID' = "Var3",
         'quadnum' = "Var4")

t2 <- test %>%
  full_join(quadnums, by = c("yrID", "quadnum")) %>%
  unite(c("Plot", "Transect", "Quadrat"),
        col = "plot_trans_quad",
        sep = "_",
        remove = F) %>%
  filter(plot_trans_quad != "NA_NA_NA")

years <- t2 %>%
  distinct(plot_trans_quad, EventYear) %>%
  group_by(plot_trans_quad) %>%
  arrange(EventYear) %>%
  mutate(yrNum = 1:n())

t3 <- t2 %>%
  left_join(years, by = c("plot_trans_quad", "EventYear"))

# Get jaccard for each iteration ------------------------------------------

sites <- unique(t2$plot_trans_quad)
iterations <- unique(t2$iter)

g <- expand.grid(sites, iterations)

jacc_fun <- function(x, data = g){
  
  site <- as.character(g[x,1])
  iteration <- as.numeric(g[x,2])
  
  matrix <- t3 %>%
    filter(plot_trans_quad == site) %>%
    filter(iter == iteration) %>%
    dplyr::select(SpecID, yrNum, value) %>%
    arrange(yrNum) %>%
    pivot_wider(names_from = yrNum,
                values_from = value,
                values_fn = list) %>%
    column_to_rownames(var = 'SpecID')
  
  a <- matrix(NA, nrow = nrow(matrix),
              ncol = ncol(matrix))
  
  b <- matrix(NA, nrow = nrow(matrix),
              ncol = ncol(matrix))
  
  c <- matrix(NA, nrow = nrow(matrix),
              ncol = ncol(matrix))
  
  for(r in 1:nrow(matrix)){
    for(t in 2:ncol(matrix)){
      #is species k shared in site i between t and t+1
      #if shared, value of a will be 1
      a[r, t] <- as.integer(matrix[r, t-1]==1 & matrix[r, t]==1)
      #is species k gained in site i between t and t+1
      #if gained, value of b will be 1
      b[r,t] <- as.integer(matrix[r,t-1] == 1 & matrix[r,t] == 0)
      #is species k lost in site i between t and t+1?
      #if lost, value of c will be 1
      c[r,t] <- as.integer(matrix[r,t-1]==0 & matrix[r,t]==1)
    }
  }
  
  A <- colSums(a)
  B <- colSums(b)
  C <- colSums(c)
  
  turnover <- (B + C)/(A + B + C)
  
  gain <- C/(A+B+C)
  
  loss <- B/(A+B+C)
  
  jacc_df <- t3 %>%
    filter(plot_trans_quad == site) %>%
    filter(iter == iteration) %>%
    distinct(yrNum, plot_trans_quad) %>%
    arrange(yrNum) %>%
    cbind(turnover) 
  
  return(jacc_df)
}


results_list <- lapply(1:nrow(g), FUN = jacc_fun)
results_df <- rbindlist(results_list, idcol = "groups")
# jaccard summary ---------------------------------------------------------


#Get summary of all of those iterations
results_sum <- results_df %>%
  left_join(years, by = c('plot_trans_quad', 'yrNum')) %>%
  dplyr::select(-yrNum) %>%
  group_by(EventYear, plot_trans_quad) %>%
  summarise(mean = mean(turnover, na.rm = T),
            sd = sd(turnover, na.rm = T))

saveRDS(file = "/scratch/atm234/nps_plants/outputs/nps_Jaccard_summary.RDS")

# jaccard samples ---------------------------------------------------------

ids <- results_df %>%
  filter(plot_trans_quad == "S02_B_3") %>%
  distinct(groups) %>%
  as_vector()

iter_sample <- sample(ids,100,replace=FALSE)

one_sitedf <- results_df %>%
  filter(plot_trans_quad == "S02_B_3") %>%
  filter(groups %in% iter_sample)

saveRDS(file = "/scratch/atm234/nps_plants/outputs/nps_Jaccard_samples_S02_B_3.RDS")

onesite_summary <- results_sum %>%
  filter(plot_trans_quad == "S02_B_3")

saveRDS(file = "/scratch/atm234/nps_plants/outputs/nps_Jaccard_summary_S02_B_3.RDS")
