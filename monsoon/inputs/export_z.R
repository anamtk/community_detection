#update converged community detection model to output
# the z matrix

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "dplyr",
                  "tidyr", "ggplot2", 
                  'mcmcplots', "stringr",
                  "coda", "htmltools",
                  "jagsUI") #mcmc output


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Import model ------------------------------------------------------------

jags <- readRDS(file = "/scratch/atm234/sbc_benthic/outputs/benthic_MSOM_1_20.RDS")



# Update model tracking just the z matrix ---------------------------------


parms <- c("z")

z.model <- update(jags,
                parameters.to.save = parms,
                n.iter = 500)


# Export the three z-matrices ---------------------------------------------

z.med <- z.model$q50$z

z.lower <- z.model$q2.5$z

z.upper <- z.model$q97.5$z

all_z <- list(z.med = z.med,
              z.lower = z.lower,
              z.upper = z.upper)

saveRDS(all_z, file = "/scratch/atm234/sbc_benthic/outputs/benthic_MSOM_1_20_zmatrix.RDS")

saveRDS(z.model$samples, "/scratch/atm234/sbc_benthic/outputs/benthic_MSOM_1_20_zmatrices.RDS")
# Export parameters of interest -------------------------------------------

vis_posterior <- jags$sims.list$a1.Vis


# Export summaries --------------------------------------------------------

jags_sum <- jags$summary
saveRDS(jags_sum,
        file = "/scratch/atm234/sbc_benthic/outputs/benthic_MSOM_1_20_summary.RDS")


