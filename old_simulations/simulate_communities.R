#Simulate communities through time
# Ana Miller-ter Kuile
# September 13, 2022

library(tidyverse)
#this script simulates communities through time 

comm <- c(1,1,1,1,1,1,1,1)
nsteps <- 10
thresh <- c(rep(0.2, 4), rep(0.4, 4))

matrix <- matrix(data = NA, nrow = nsteps, ncol = length(comm))

 for(i in 1:nsteps){
   randdraw <- rnorm(length(comm))
   trans <- ifelse(randdraw > thresh, 1, 0)
   comm <- ifelse(trans == 1, abs(comm-1), comm)
   matrix[i,] <- comm
 }

as.data.frame(matrix) %>%
  rownames_to_column(var = "timestep") %>%
  pivot_longer(V1:V8,
               names_to = "species",
               values_to = "occupancy") %>%
  ggplot(aes(x = timestep, y = occupancy)) +
  geom_point() +
  facet_grid(~species)
