model{
  for(k in 1:n.species){
    for(i in 1:n.transects){
      #year 1 true occupancy
      #year one occupancy probability is
      #based on a site-level random effects intercept
      logit(psi1[k,i]) <- b0[k,SiteID[i]]#site-level random effects intercept
      
      #true occupancy is bernoulli around first-year
      #occupancy probability
      z[k,i,1] ~ dbern(psi1[k,i])
      
      #year 2+ true occupancy
      for(t in 2:n.years){
        #occupancy dependent on persistence, phi
        # and colonization, gamma
        #persistence only applies if previous year z = 1
        # colonization only applies if previous year z = 0
        logit(phi[k,i]) <- c0[k,SiteID[i]] #site-level random effects intercept
        logit(gamma[k,i]) <- d0[k,SiteID[i]] #site-level random effects intercept
        
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
      
      #occupancy intercept
      b0[k,SiteID[i]] ~ dnorm(mu.b0, tau.b0)
      #persistence intercept
      c0[k,SiteID[i]] ~ dnorm(mu.c0, tau.c0)
      #colonization intercept
      d0[k,SiteID[i]] ~ dnorm(mu.d0, tau.d0)

    # Zero-centered hierarchical prior for random effects, just
    # as used in the original model (do not monitor or report these):
    #initial occupancy
    eps.b0[k, SiteID[i]] ~ dnorm(0,tau.eps.b0)
    # Compute identifiable random effects (monitor and report these,
    # if desired):
    eps.b0.star[k, SiteID[i]] <- eps.b0[k, SiteID[i]] - mean.eps.b0

    #persistence
    eps.c0[k, SiteID[i]] ~ dnorm(0,tau.eps.c0)
    # Compute identifiable random effects (monitor and report these,
    # if desired):
    eps.c0.star[k, SiteID[i]] <- eps.c0[k, SiteID[i]] - mean.eps.c0
    
    #colonization
    eps.d0[k, SiteID[i]] ~ dnorm(0,tau.eps.d0)
    # Compute identifiable random effects (monitor and report these,
    # if desired):
    eps.d0.star[k, SiteID[i]] <- eps.d0[k, SiteID[i]] - mean.eps.d0
    
    } #transects likelihood loop
    
    
    #Species-level PRIORS
    #detection intercept and slopes
    a0[k] ~ dnorm(mu.a0, tau.a0)
    a1.Vis[k] ~ dnorm(mu.a1, tau.a1)
    a2.Size[k] ~ dnorm(mu.a2, tau.a2)
  } #species loop
  
  #missing data 
  for(i in 1:n.years){
    for(j in 1:n.reps){
      #missing data in the visibility column
      vis[i,j] ~ dnorm(mu.missingvis, tau.missingvis)
    }
  }
  
  #Community-level hyperpriors
  
  # Prior for non-identifiable intercept (don't monitor or report)
  #initial occupancy
  mu.b0 ~ dnorm(0,1E-6)
  # Compute identifiable intercept (monitor and report this):
  b0.star <- mu.b0 + mean.eps.b0
  # Mean of non-identifiable random effects (required)
  mean.eps.b0 <- mean(eps.b0[])
  
  #persistence
  mu.c0 ~ dnorm(0,1E-6)
  # Compute identifiable intercept (monitor and report this):
  c0.star <- mu.c0 + mean.eps.c0
  # Mean of non-identifiable random effects (required)
  mean.eps.c0 <- mean(eps.c0[])
  
  #colonization
  mu.d0 ~ dnorm(0,1E-6)
  # Compute identifiable intercept (monitor and report this):
  d0.star <- mu.d0 + mean.eps.d0
  # Mean of non-identifiable random effects (required)
  mean.eps.d0 <- mean(eps.d0[])
  

  #Detection intercept and slope community priors
  a0.mean ~ dbeta(1,1)
  mu.a0 <- logit(a0.mean)
  tau.a0 <- pow(sig.a0, -2)
  sig.a0 ~ dunif(0, 10)
  
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