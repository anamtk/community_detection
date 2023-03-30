#DYNAMIC MSOM FOR KELP FISH
#Ana Miller-ter Kuile
#March 30, 2023

# This script is a multi-site, dynamic (multi-year) multi-species model 
# for fish in the kelp forest in the SB LTER. 
# This script considers first year occupancy and then any following year
# occupancy as dependent on the previous year occupancy state and 
# persistence and colonization probabilities.

# Right now, I'm modeling transect-level data, but I would really like to 
# be able to add a random hierarchy of transects within a reef site - not
# sure where to put this random effect into the model - might be a good
# convo for a future meeting with Kiona

# Right now this model has these aspects:
## 1. Multi-site, multi-species, multi-year
## 2. Persistence and colonization are constant across years
## 3. No covariates for initial occupancy, persistence, or colonization
## 4. Covariates on detection including average size of the fish species 
### and visibility during the survey period

model{
  for(k in 1:n.species){
    for(i in 1:n.transects){
      #year 1 true occupancy
      #year one occupancy probability is
      #based on initial occupancy probability, psi1 for species k in
      # site i
      #true occupancy is bernoulli around first-year
      #occupancy probability
      z[k,i,1] ~ dbern(psi1[k,i])
      
      #year 2+ true occupancy
      for(t in 2:n.years){
        #occupancy dependent on persistence, phi
        # and colonization, gamma
        #persistence only applies if previous year z = 1
        # colonization only applies if previous year z = 0
        z[k,i,t] ~ dbern(z[k,i, t-1]*phi[k,i] +
                           (1- z[k,i, t-1])*gamma[k,i])
      }#Year 2+ loop
      for(t in 1:n.years){
        # Observation model
        logit(p[i,k,t]) <- a0[k] + #species-level intercept
          a1.Vis[k]*vis[i,t] + #visibility effect
          a2.Size[k]*size[k] #size effect
        #might also consider a "average group size covariate in here
        # down the road?
        
        #detection is bionmiaml based on detection probability times
        # true occupancy, with the number of trials being the 
        # nubmer of visits to that site that year
        y[i, k] ~ dbin(p[i,k,t] * z[i, k,t], n.reps)
        
      } #years detection model loop
    } #transects likelihood loop
    
    #Species-level PRIORS
    # Species level priors for occupancy, persistence, colonization,
    # and detection are centered around community-level priors for each
    # of these variables
    
    #occupancy
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi) 
    psi[k] <- ilogit(lpsi[k])
    
    #persistence
    lphi[k] ~ dnorm(mu.lphi, tau.lphi)
    phi[k] <- ilogit(lphi[k])
    
    #colonization
    lgamma[k] ~ dnorm(mu.lgamma, tau.lgamma)
    gamma[k] <- ilogit(lgamma[k])

    #detection intercept and slopes
    a0[k] ~ dnorm(mu.a0, tau.a0)
    a1.Vis[k] ~ dnorm(mu.a1, tau.a1)
    a2.Size[k] ~ dnorm(mu.a2, tau.a2)
    
  } #species loop
  
  #missing data 
  #SOme data for visibility are missing, so we're imputing them
  for(i in 1:n.years){
    for(j in 1:n.reps){
      #missing data in the visibility column
      vis[i,j] ~ dnorm(mu.missingvis, tau.missingvis)
    }
  }
  
  #Community-level hyperpriors
  #All species-level priors are centered around hyperpriors for 
  # the community for that variaable
  
  #initial occupancy
  psi.mean ~ dbeta(1,1)
  mu.lpsi <- logit(psi.mean)
  sd.lpsi ~ dunif(0, 10)
  tau.lpsi <- pow(sd.lpsi, -2)
  
  #persistence
  phi.mean ~ dbeta(1,1)
  mu.lphi <- logit(phi.mean)
  sd.lphi ~ dunif(0, 10)
  tau.lphi <- pow(sd.lphi, -2)
  
  #colonization
  gamma.mean ~ dbeta(1,1)
  mu.lgamma <- logit(gamma.mean)
  sd.lgamma ~ dunif(0, 10)
  tau.lgamma <- pow(sd.lgamma, -2)
  
  #Detection intercept
  a0.mean ~ dbeta(1,1)
  mu.a0 <- logit(a0.mean)
  tau.a0 <- pow(sig.a0, -2)
  sig.a0 ~ dunif(0, 10)
  
  #covariate means
  mu.a1 ~ dunif(-5,5)
  tau.a1 <- pow(sig.a1, -2)
  sig.a1 ~ dunif(0, 10)
  
  mu.a2 ~ dunif(-5,5)
  tau.a2 <- pow(sig.a2, -2)
  sig.a2 ~ dunif(0, 10)
  
  #PRIORS FOR IMPUTING MISSING DATA
  #Priors for mean and tau of missing covariates in the model
  mu.missingvis ~ dunif(-10, 10)
  sig.missingvis ~ dunif(0, 20)
  tau.missingvis <- pow(sig.missingvis, -2)
  
} #whole model