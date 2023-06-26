# Multi-species occupancy model for SBC-LTER benthic surveys
# Ana Miller-ter Kuile
# January 17, 2022

# this is a multi-species occupancy model for the SBC-LTER benthic surveys
# that we developed to try and understand how detection error might shape
# expectations about community dynamics through time (and potentially space)
# WRT species richness, species composition, and potentially functional
# composition, which might be the most important piece here

# Attributes of this model include:
# 1. Biological process model:
## - first-year occupancy is bernoulli distributed
## - no covariates currently
## - true occupancy, z, is shaped by persistence (phi) and colonization (gamma)
### dynamics, which both vary by year
# 2. Detection process:
## - given lack of intra-annual re-surveys of transects, spatial replicates of 
### four quadrats on a given transect, with transect occupancy being the unit
### we are estimating
## - covariate of visibility on that dive day for whether the species is detected
## - mean detection is driven by detection probability and true occupancy
## - observed data follow bernoulli distribution around mean detection probability
# 3. Multi-species components:
## - species-level priors for first year occupancy, population persistence, 
### probability of coloniazation, intercept and visibility covariates
## - community-level hyperpriors for population persistence, probability of 
### colonization, and covariate priors

######
model{
  for(i in 1:n.species){
      #first-year occupancy
      z[i, 1] ~ dbern(psi1[i]) #initial occupancy for that species
      for(t in 2:n.years){ #years indexed by j right now bc i think transects have diff. # survey years
        #biological process
        #B) no covariates on biological process:
        #phi is the probability population persists if z[t-1] = 1
        #and gamma is the probability of colonization if z[t-1] = 0
        #true occupancy state is based on this occupancy
        # probability that is shaped by persistence and colonization
        z[i,t] ~ dbern(z[i,t-1]*phi[i,t-1] +
          (1 - z[i,t-1]*gamma[i,t-1]))
        #right now persistence and colonization rates are 
        #dependent on species and year

         } #biological process year
      #Detection process
      for(t in 1:n.years){
        for(k in 1:n.reps){ #replicates are spatial in this model
          # detection probability is shaped by covariate of visibility
          logit(p[i,t,k]) <- a0[i] + aVis[i]*vis[t]
          #mean detection probability is a product of this detection
          # probability and the true occupancy state
          mu.p[i,t,k] <- p[i,t,k]*z[i,t]
          #observed data follow a bernoulli distribution around mean
          # detection probability
          y[i,t,k] ~ dbin(mu.p[i,t,k],1)
        } #detection quadrat loop
      } #detection Year loop
    
    ##PRIORS##
    
    ## Species-level priors
    
    #persistence and colonization are currently
    #dependent on species and year
    for(t in 1:(n.years-1)){
      phi[i,t] ~ dnorm(mu.phi, tau.phi)
      gamma[i,t] ~ dnorm(mu.gamma, tau.gamma)
    } #persistence/colonization year loop
    
    #initial occupancy for species i
    # is centered around the community mean 
    # and SD for occupancy
    lpsi1[i] ~ dnorm(mu.lpsi, tau.lpsi)
    psi1[i] <- ilogit(lpsi1[i])
    #species-level detection covariate priors: 
    a0[i] ~ dnorm(mu.a0, tau.a0)
    aVis[i] ~ dnorm(mu.Vis, tau.Vis)
    
  } #species
  
  #Community-level hyperpriors
  lpsi.mean ~ dbeta(1,1)
  mu.lpsi <- logit(lpsi.mean)
  sd.lpsi ~ dunif(0, 5)
  tau.lpsi <- 1/pow(sd.lpsi,2)
  
  #for detection intercept
  a0.mean ~ dbeta(1,1)
  mu.a0 <- logit(a0.mean)
  sd.a0 ~ dunif(0, 5)
  tau.a0 <- 1/pow(sd.a0,2)
  
  #for visibility covariate
  mu.Vis ~ dunif(-5,5)
  sd.Vis ~ dunif(0,5)
  tau.Vis <- 1/pow(sd.Vis,2)
  
  #for persistence
  phi.mean ~ dbeta(1, 1)
  mu.phi <- logit(phi.mean)
  sd.phi ~ dunif(0, 5)
  tau.phi <- 1/pow(sd.phi,2)
  
  #for colonization
  gamma.mean ~ dbeta(1,1)
  mu.gamma <- logit(gamma.mean)
  sd.gamma ~ dunif(0, 5)
  tau.gamma <- 1/pow(sd.phi,2)
  
}