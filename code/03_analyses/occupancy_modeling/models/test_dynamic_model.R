model{
  for(k in 1:n.species){  # Loop through species
    #initial occupancy
    z0[k] ~ dbern(gamma0[k])
    logit(psi[1,k]) <- lpsi[k] + gamma[k]*z0[k]
    z[1,k] ~ dbern(psi[1,k])
    # Likelihood
    for(i in 2:n.years) {
      # Ecological model
      logit(psi[i,k]) <- lpsi[k]+ gamma[k]*z[i-1,k]
      z[i, k] ~ dbern(psi[i,k])
      # Observation model
      logit(p[i,k]) <- eta[k] + a1.Vis[k]*vis[i]
      y[i, k] ~ dbin(p[i,k] * z[i, k], n.reps)
    }
    

    # Priors (species level)
    # for(t in 1:n.years){
    # psi[t,k] <- ilogit(lpsi[k])
    # }
    #links detection and occupancy probabilities to each other
    #(more abundant = more observed)
    
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    mu.eta[k] <- mu.lp + rho * sd.lp/sd.lpsi * 
      (lpsi[k] - mu.lpsi)
    
    eta[k] ~ dnorm(mu.eta[k], tau.eta)
    a1.Vis[k] ~ dnorm(mu.vis, tau.vis)
    
    # species specific random effects
    gamma0[k] ~ dbeta(1, 1)
    gamma[k] ~ dnorm(mugamma, taugamma)

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
  
  p_gamma ~ dbeta(1, 1)
  mugamma <- log(p_gamma / (1 - p_gamma))
  sigmagamma~dunif(0,10)
  taugamma<-1/(sigmagamma*sigmagamma)
}