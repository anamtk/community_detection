model{
  for(k in 1:n.species){
    for(i in 1:n.transects){
      #BIOLOGICAL MODEL
      z[k,i,1] ~ dbern(psi1[k])
      
      #year 2+ true occupancy
      for(t in 2:n.years){
        z[k,i,t] ~ dbern(z[k,i, t-1]*phi[k] +
                           (1- z[k,i, t-1])*gamma[k])
      }#Year 2+ loop
      
      #OBSERVATION MODEL
      
      for(t in 1:n.years){
        for(r in 1:n.rep[i,t]){
          y[k,i,t,r] ~ dbin(p[k] * z[k,i,t], n.rep[i,t])
          
        } #reps detection loop
      } #years detection model loop
    }#transects likelihood loop
    
    #Species-level PRIORS
    # Species level priors for occupancy, persistence, colonization,
    # and detection are centered around community-level priors for each
    # of these variables
    
    #occupancy
    
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi) 
    psi1[k] <- ilogit(lpsi[k])
    
    #persistence
    lphi[k] ~ dnorm(mu.lphi, tau.lphi)
    phi[k] <- ilogit(lphi[k])
      
    lgamma[k] ~ dnorm(mu.lgamma, tau.lgamma)
    gamma[k] <- ilogit(lgamma[k])
    
    lp[k] ~ dnorm(mu.lp, tau.lp)
    p[k] <- ilogit(lp[k])
    
  } #species loop
  
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
  
  p.mean ~ dbeta(1,1)
  mu.lp <- logit(p.mean)
  sd.lp ~ dunif(0, 10)
  tau.lp <- pow(sd.lp, -2)
  
} 