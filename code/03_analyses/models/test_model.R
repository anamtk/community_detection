model{
  for(k in 1:n.species){  # Loop through species
    # Likelihood
    for(i in 1:n.years) {
      # Ecological model
      z[i, k] ~ dbern(psi[k])
      # Observation model
      logit(p[i,k]) <- eta[k]*a1.Vis[k]*vis[i]
      y[i, k] ~ dbin(p[i,k] * z[i, k], n.reps)
    }
    
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
  }
  
  # Hyperpriors (community level)
  psi.mean ~ dbeta(1, 1)
  mu.lpsi <- logit(psi.mean)
  sd.lpsi ~ dunif(0, 5)
  tau.lpsi <- 1/sd.lpsi^2
  
  p.mean ~ dbeta(1, 1)
  mu.lp <- logit(p.mean)
  sd.lp ~ dunif(0, 5)
  tau.lp <- 1/sd.lp^2
  
  mu.vis ~ dunif(-5,5)
  sd.vis ~ dunif(0,5)
  tau.vis <- pow(sd.vis, -2)
  
  rho ~ dunif(-1, 1)
  tau.eta <- tau.lp/(1 - rho^2)
}