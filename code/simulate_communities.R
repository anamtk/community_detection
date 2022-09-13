#Simulate communities through time
# Ana Miller-ter Kuile
# September 13, 2022

#this script simulates communities through time 

comm <- c(0,0,0,1,0,1,1,1)
nsteps <- 10
thresh <- 0.2

matrix <- matrix(data = NA, nrow = nsteps, ncol = length(comm))

 for(i in 1:nsteps){
   randdraw <- runif(length(comm))
   trans <- ifelse(randdraw > thresh, 1, 0)
   comm <- ifelse(trans == 1, abs(comm-1), comm)
   matrix[i,] <- comm
 }


#THIS is not quite working yet- keep working on it
comm <- c(1,1,1,1,0,1,1,1)
nsteps <- 10
thresh <- 0.2

matrix <- matrix(data = NA, nrow = nsteps, ncol = length(comm))

for(i in 1:nsteps){
  #random variation in occupancy
  randdraw <- runif(length(comm))
  #occupancy baseline to simulate "common" vs "rare"
  occupancy <- c(rep(2, 4), rep(1, 4))
  #multiply the occupancy by the random draw - which will
  #always be greater for the more common species 
  occchange <- occupancy*randdraw
  #if the change < threshhold value set (smaller = more variable),
  #then transition likelihood is 1, otherwise 0
  trans <- ifelse(occchange < thresh, 1, 0)
  #set 1-0 of community based on previous community based on
  #whether transition is = 1 or 0
  comm <- ifelse(trans == 1, abs(comm-1), comm)
  #populate the ith row in the matrix with the new community
  matrix[i,] <- comm
}

matrix
# for(i in 1:nsteps){
#   comm <- rbinom(length(comm), size =3 , prob = c(rep(0.8, 4), rep(0.3,4)))
#   matrix[i,] <- comm
# }
