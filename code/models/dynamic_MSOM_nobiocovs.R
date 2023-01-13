#species s on transect t in year y
#detection: species s on quadrat q on transect t in year y

######
model{
  for(j in 1:n.transects){
    for(i in 1:n.species){
      #first-year occupancy
      z[j,i, 1] ~ dbern(psi1[j,i]) #occupancy for that species
      for(t in 2:n.years[j]){ #years indexed by j right now bc i think transects have diff. # survey years
        #biological process
        #B) no covariates on biological process:
        #phi is the probability population persists if z[t-1] = 1
        z[j,i,t] ~ dbern(z[j,i,t-1]*phi[j,i,t-1] +
                           #gamma is the probability of colonization
                           # if z[t-1] = 0
                           (1 - z[j,i,t-1]*gamma[j,i,t-1]))

         } #biological process year
      #Detection process
      for(t in 1:n.years[j]){
        for(k in 1:n.quads[j]){
          logit(p[j,i,t,k]) <- a0[i] + aVis[i]*vis[j,t]
          
          mu.p[j,i,t,k] <- p[j,i,t,k]*z[j,i,t]
          y[j,i,t,k] ~ dbern(mu.p[j,i,t,k])
        } #detection quadrat loop
      } #detection Year loop
      
      #Species-level priors
      #initial occupancey for species j on site i
      # is centered around the community mean 
      # and SD for occupancy
      lpsi1[j,i] ~ dnorm(mu.lpsi, tau.lpsi)
      psi1[j,i] <- ilogit(lpsi[j,i])
      #species-level detection covariate priors: 
      a0[i] ~ dnorm(mu.a0, tau.a0)
      aVis[i] ~ dnorm(mu.Vis, tau.Vis)
      
    } #species

  } #transect
  
  #Community-level hyperpriors
  a0.mean ~ dbeta(1,1)
  mu.a0 <- logit(a0.mean)
  sd.a0 ~ dunif(0, 5)
  tau.a0 <- 1/sd.a0^2
  
  mu.Vis ~ dunif(-5,5)
  sd.Vis ~ dunif(0,5)
  tau.Vis <- 1/sd.Vis^2
  
} #model