#simplified from: https://github.com/valentinlauret/IntegratedSingleVisitOccupancy/blob/master/dolphins_codes/Aerial_SV.R
model {
  for(s in 1:n.species){
  #biological process
    for (i in 1:nsite){
      #true state follows bernoulli of occupancy probability
    z[s, i] ~ dbern(psi[s])
      #the logit link of that psi can have covariates
      # in either one or two rows of code as follows:
    logit(psi[s]) <- b0[s] + b1[s]*biological_covariate[i] 
    
    #Detection process
    #mean occupancy probability is true state times
    #occupancy probability
    mu.p[s] <- z[s,i] * p[s] 
    #logit link of occupancy probability is 
    #set of intercept plus effort covariate
    logit(p[s]) <-  a0[s] + a1[s]*detection_covariate[i]
    #observed data follow bernoilli with mean
    #detection probability
    y[s,i] ~ dbern(mu.p[s])
    } #i
  
  # priors
  #species-level priors
  #priors for intercepts
  b0[s] ~ dnorm(mu.b0, tau.b0)
  a0[s] ~ dnorm(mu.a0, tau.a0)
  ## WORK FROM HERE ON UPDATING PRIORS FOR
  # COMMUNITY PRIORS##
  mu.eta[s] <- mu.lp + rho * sd.lp/sd.b0 *
    (b0[s] - mu.b0)
  mu.eta2[s] <- mu.lp2 + rho2 *sd.lp2/sd.a0 *
    (a0[s] - mu.a0)
  
  lp[s] ~ dnorm(mu.eta[s], tau.eta)
  lp2[s] ~ dnorm(mu.eta2[s], tau.eta)
  
  #priors for covariates
  b1[s] ~ dnorm(mu.b1, tau.b1)
  a1[s] ~ dnorm(mu.a1, tau.a1)
  } # species loop
  
  #biological process covariates
  b0 ~ dnorm(0, 0.001)
  for(bb in 1:n.beta){
    b[bb] ~ dnorm(0, 0.001)
  }
  #detection process covariates
  a0 ~ dnorm(0, 0.001)
  for(aa in 1:n.alpha){
    a[aa] ~ dnorm(0, 0.001)
  }

  
  # Priors (species level)
  lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
  psi[k] <- ilogit(lpsi[k])
  
  lp[k] ~ dnorm(mu.lp, tau.lp)
  p[k] <- ilogit(lp[k])
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
}