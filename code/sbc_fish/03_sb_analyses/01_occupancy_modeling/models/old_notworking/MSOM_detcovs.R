model{
  for(k in 1:n.species){  # Loop through species
    # Likelihood
    for(i in 1:n.years) {
      # Ecological model
      z[k,i] ~ dbern(psi[k,i])
      # Observation model
      for(j in 1:n.reps){
        logit(p[k,i,j]) <- a0[k] + aVis[k]*vis[i]
        
      y[k, i, j] ~ dbern(p[k,i,j] * z[k,i])
      }
    }
    
    #logit psi for species k is normal dist around
    # community mean
    lpsi[k] ~ dnorm(mu.eta[k], tau.eta)
    #mu.eta links psi to p 
    mu.eta[k] <- mu.a0 + rho * sd.a0/sd.lpsi * 
      (lpsi[k] - mu.lpsi)
    
    psi[k] <- ilogit(lpsi[k])
    
    #covariates on detection
    a0[k] ~ dnorm(mu.a0, tau.a0)
    aVis[k] ~ dnorm(mu.vis, tau.vis)
    
  }
  
  # Hyperpriors (community level)
  a0.mean ~ dbeta(1, 1)
  mu.a0 <- logit(a0.mean)
  sd.a0 ~ dunif(0, 5)
  tau.a0 <- 1/sd.a0^2
  
  mu.vis ~ dunif(-5, 5)
  sd.vis ~ dunif(0, 5)
  tau.vis <- 1/sd.vis^2

  psi.mean ~ dbeta(1, 1)
  mu.lpsi <- logit(psi.mean)
  sd.lpsi ~ dunif(0, 5)
  tau.lpsi <- 1/sd.lpsi^2
  
  rho ~ dunif(-1, 1)
  tau.eta <- tau.lpsi/(1 - rho^2)
}