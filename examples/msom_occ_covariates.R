# from: https://bcss.org.my/tut/community-multi-species-models/multi-species-occupancy-with-covariates/
# MSOM model with covariates
# File name "MSOM_cov.jags"

model{
  for(k in 1:nSobs){  # Loop through species
    # Likelihood
    for(i in 1:nSites) {
      # Ecological model
      logit(psi[i, k]) <- b0[k] + bArea[k] * logArea[i] +
        bDist[k] * MLdist[i]
      z[i, k] ~ dbern(psi[i, k])
      # Observation model
      y[i, k] ~ dbin(p[k] * z[i, k], nOcc)
    }
    
    # Priors (species level)
    b0[k] ~ dnorm(mu.b0, tau.b0)
    mu.eta[k] <- mu.lp + rho * sd.lp/sd.b0 *
      (b0[k] - mu.b0)
    lp[k] ~ dnorm(mu.eta[k], tau.eta)
    p[k] <- ilogit(lp[k])
    
    bArea[k] ~ dnorm(mu.area, tau.area)
    bDist[k] ~ dnorm(mu.dist, tau.dist)
    
  }
  
  # Hyperpriors (community level)
  b0.mean ~ dbeta(1, 1)
  mu.b0 <- logit(b0.mean)
  sd.b0 ~ dunif(0, 5)
  tau.b0 <- 1/sd.b0^2
  
  mu.area ~ dunif(-5, 5)
  sd.area ~ dunif(0, 5)
  tau.area <- 1/sd.area^2
  mu.dist ~ dunif(-5, 5)
  sd.dist ~ dunif(0, 5)
  tau.dist <- 1/sd.dist^2
  
  p.mean ~ dbeta(1, 1)
  mu.lp <- logit(p.mean)
  sd.lp ~ dunif(0, 5)
  tau.lp <- 1/sd.lp^2
  
  rho ~ dunif(-1, 1)
  tau.eta <- tau.lp/(1 - rho^2)
}