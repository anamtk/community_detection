#species s on transect t in year y
#detection: species s on quadrat q on transect t in year y

model{
  for(s in 1:n.species){
    for(i in 1:n.transects){
      for(y in 1:n.years[i]){ #years indexed by t right now bc i think transects have diff. survey years
        #biological process
        #A) with covariates:
        #logit(psi[s,i,y]) <- b0[s] + b1*cov1[s,i,y]
        #true occupancy
        #z[s,i,y] ~ dbern(psi[s,i,y])
        #B) without covariates:
        z[s,i,y] ~ dbern(psi[s,i,y])
        #detection process
        for(q in 1:n.quads[i]){
          logit(p[s,i,q,y]) <- a0[s] + aVis[s]*vis[i,y]
          #y[s,i,q,y] ~ dbin(p[s,i,q,y] * z[s,i,y], nOcc[i])
          mu.p[s, i, q,y] <- p[s,i,q,y]*z[s,i,y]
          y[s,i,q,y] ~ dbern(mu.p[s,i,q,y])
        } #replicate loop
        
        #Species-level priors
        lpsi[s,i,y] ~ dnorm(mu.lpsi, tau.lpsi)
        psi[s,i,y] <- ilogit(lpsi[s,i,y])
      } #year loop
      
    } #transect loop
    #Species-level priors
    a0[s] ~ dnorm(mu.a0, tau.a0)
    aVis[s] ~ dnorm(mu.Vis, tau.Vis)
    
  } #species loop
  
  #Community-level priors
  a0.mean ~ dbeta(1,1)
  mu.a0 <- logit(a0.mean)
  sd.a0 ~ dunif(0, 5)
  tau.a0 <- 1/sd.a0^2
  
  mu.Vis ~ dunif(-5,5)
  sd.Vis ~ dunif(0,5)
  tau.Vis <- 1/sd.Vis^2
  

} #model loop
