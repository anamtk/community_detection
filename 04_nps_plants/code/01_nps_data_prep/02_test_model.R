# Running the MSOM for birds
# Ana Miller-ter Kuile
# September 8, 2023

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  "jagsUI",
                  'rjags',
                  'mcmcplots',
                  "coda", 'reshape2') #mcmc output


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load Data ---------------------------------------------------------------

data <- readRDS(here('04_nps_plants',
                     'data_outputs',
                     'MSAM',
                     'model_inputs',
                     'nps_msam_multisite_subset.RDS'))

# data <- list(n.species = data$n.species,
#              n.quads = data$n.quads,
#              n.yr = data$n.yr,
#              n.rep = data$n.rep,
#              cover = data$cover,
#              lifegroup = data$lifegroup,
#              n.groups = data$n.groups,
#              y = data$y,
#              z = data$z)

# Set Initials ------------------------------------------------------------

inits <- list(list(z = data$z),
              list(z = data$z),
              list(z = data$z))

#ideally you would use the omega.init object to set inits, but it 
#seems to break something with this dataset in particular,
#let me know if this is happening to you and we can troubleshoot
#inits <- function() list(omega = omega.init)

# Parameters to save ------------------------------------------------------

params <- c(
  'b0.star',
  'eps.site.star',
  'eps.year.star',
  'mu.b0species',
  'sig.b0species',
  'sig.eps.site',
  'sig.eps.year',
  'mu.a0',
  'sig.a0',
  'a1.Cover',
  'a2.LifeGroup'
)


# JAGS model --------------------------------------------------------------

model <- here("04_nps_plants",
              'code',
              '02_nps_analyses',
              'jags',
              #"nps_MSOM_simple_SL.R")
              "nps_MSOM_simple_yrsite.R")

Sys.time()
start<-proc.time()
mod <- jagsUI::jags(data = data,
                    #inits = inits,
                    inits = NULL,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 1,
                    DIC = TRUE)

Sys.time()
end<-proc.time()
end-start

mcmcplot(mod$samples)
gelman.diag(mod$samples, multivariate = F)

#get 3-d array into a dataframe, dimensions need to be named
#https://community.rstudio.com/t/fastest-way-to-convert-3d-array-to-matrix-or-data-frame/38398/5
parms2 <- c("z")

mod2 <- update(mod,
               parameters.to.save = parms2,
               n.iter = 1)

samples <- mod2$sims.list[[1]]

test <- melt(samples) %>%
  rename("iter" = "Var1",
         'SpecID' = "Var2",
         'yrID' = "Var3",
         'quadnum' = "Var4")

t2 <- test %>%
  left_join(quadnums, by = c("yrID", "quadnum")) %>%
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

sites <- unique(t3$plot_trans_quad)
iterations <- unique(t3$iter)

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
results_df <- do.call(rbind, results_list)

results_sum <- results_df %>%
  left_join(years, by = c('plot_trans_quad', 'yrNum')) %>%
  dplyr::select(-yrNum) %>%
  group_by(EventYear, plot_trans_quad) %>%
  summarise(mean = mean(turnover, na.rm = T),
            sd = sd(turnover, na.rm = T))
