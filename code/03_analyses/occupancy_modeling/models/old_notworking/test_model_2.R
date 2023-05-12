model{
  for(k in 1:n.species){  # Loop through species
    # Likelihood
    for(i in 1:n.years) {
      # Ecological model
      z[k,i] ~ dbern(psi[k])
      # Observation model
      for(j in 1:n.reps){
      logit(p[k,i,j]) <- a0[k] + 
                        a1.Vis[k]*vis[i,j] + 
                        a2.Size[k]*size[k]
      
      y[k,i, j] ~ dbin(p[k,i,j] * z[k,i], n.reps)
      

      } #replicate visits
    } #years
    
    # Priors (species level)
    #occupancy probability
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    psi[k] <- ilogit(lpsi[k])
    
    #detection intercept
    a0[k] ~ dnorm(mu.a0, tau.a0)
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

  a0.mean ~ dbeta(1, 1)
  mu.a0 <- logit(a0.mean)
  sd.a0 ~ dunif(0, 5)
  tau.a0 <- pow(sd.a0, -2)
  
  mu.vis ~ dunif(-5,5)
  sd.vis ~ dunif(0,5)
  tau.vis <- pow(sd.vis, -2)
  
  mu.size ~ dunif(-5,5)
  sd.size ~ dunif(0, 5)
  tau.size <- pow(sd.size, -2)
  
  #PRIORS FOR IMPUTING MISSING DATA
  #Priors for mean and tau of missing covariates in the model
  mu.missingvis ~ dunif(-10, 10)
  sig.missingvis ~ dunif(0, 20)
  tau.missingvis <- pow(sig.missingvis, -2)
}