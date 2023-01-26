model{
  for(k in 1:n.species){  # Loop through species
    # Likelihood
    for(i in 1:n.years) {
      # Ecological model
      z[k,i] ~ dbern(psi[k])
      # Observation model
      for(j in 1:n.reps){
      logit(p[k,i,j]) <- eta[k] + 
                        a1.Vis[k]*vis[i,j] + 
                        a2.Size[k]*size[k]
      
      y[k,i, j] ~ dbin(p[k,i,j] * z[k,i], n.reps)
      

      } #replicate visits
    } #years
    
    # Priors (species level)
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    psi[k] <- ilogit(lpsi[k])
    
    #links detection and occupancy probabilities to each other
    #(more abundant = more observed)
    mu.eta[k] <- mu.lp + rho * sd.lp/sd.lpsi * 
      (lpsi[k] - mu.lpsi)
    
    eta[k] ~ dnorm(mu.eta[k], tau.eta)
    #lp[k] ~ dnorm(mu.eta[k], tau.eta)
    #p[k] <- ilogit(lp[k])
    a1.Vis[k] ~ dnorm(mu.vis, tau.vis)
    a2.Size[k] ~ dnorm(mu.size, tau.size)
  }
  
  for(i in 1:n.years){
    for(j in 1:n.reps){
  #missing data in the visibility column
  vis[i,j] ~ dnorm(mu.missingvis, tau.missingvis)
    }
    }
  
  # Hyperpriors (community level)
  psi.mean ~ dbeta(1, 1)
  mu.lpsi <- logit(psi.mean)
  sd.lpsi ~ dunif(0, 5)
  tau.lpsi <- pow(sd.lpsi, -2)
  
  #does this part need to be deterministic
  #linked somehow to the p above related
  #to the regression of visibility and body size?
  p.mean ~ dbeta(1, 1)
  mu.lp <- logit(p.mean)
  sd.lp ~ dunif(0, 5)
  tau.lp <- pow(sd.lp, -2)
  
  mu.vis ~ dunif(-5,5)
  sd.vis ~ dunif(0,5)
  tau.vis <- pow(sd.vis, -2)
  
  mu.size ~ dunif(-5,5)
  sd.size ~ dunif(0, 5)
  tau.size <- pow(sd.size, -2)
  
  rho ~ dunif(-1, 1)
  tau.eta <- tau.lp/(1 - rho^2)
  
  #PRIORS FOR IMPUTING MISSING DATA
  #Priors for mean and tau of missing covariates in the model
  mu.missingvis ~ dunif(-10, 10)
  sig.missingvis ~ dunif(0, 20)
  tau.missingvis <- pow(sig.missingvis, -2)
}