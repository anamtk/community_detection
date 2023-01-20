model {
  
  #data
  for(i in 1:(n.species+n.zeros)){
    for(k in 1:n.years){
      for(j in 1:n.reps){
        Ytemp[i,k,j]
      }
    }
  }
  
  for (i in 1:(n.species+n.zeros)) {
    for (k in 1:n.years) {
      #logit detection is linked to phi, which ties occ. to det. probability
      logit(psi[i,k]) <- phi[i]
      #mean occupancy for species i in site k is related to occupancy times 
      # community inclusion for species i
      mu.psi[i,k] <- psi[i,k]*w[i]
      #true occupancy for species follows a bernoilli distribution
      #with occupancy probabiliy mu.psi
      Z[i,k] ~ dbern(mu.psi[i,k]) #occupancy
      
      for (j in 1:n.reps){
          logit(p[i,k,j]) <- eta[i]+a1.Vis[i]*vis[k]	#detection
          mu.p[i,k,j] <- p[i,k,j]*Z[i,k]

        #observed data
        Ytemp[i,k,j] ~ dbern(mu.p[i,k,j]) 
      }		
    }
    
    #PRIORS
    #SPecies-level
    #original membership in the community
    w[i] ~ dbern(omega)				#membership in community
    #covariate for occupancy is phi - which links occupancy and detection
    phi[i] ~ dnorm(beta, tau.u)
    a1.Vis[i] ~ dnorm(mu.vis, tau.vis)
    
    #mu.eta and eta links occupancy to detection
    mu.eta[i] <- alpha + (rho*sigma.v/sigma.u)*(phi[i] - beta)
    #links detection to occupancy
    eta[i] ~ dnorm(mu.eta[i], var.eta)
  }
  
  #Community-level hyperpriors
  omega ~ dunif(0,1)
  
  psi.mean ~ dunif(0.0001,1)
  beta <- log(psi.mean) - log(1-psi.mean)
  sigma.u ~ dunif(0,5)
  tau.u <- pow(sigma.u,-2)
  
  p.mean ~ dunif(0.0001,1)
  alpha <- log(p.mean) - log(1-p.mean)
  sigma.v ~ dunif(0,5)
  tau.v <- pow(sigma.v,-2)
  
  mu.vis ~ dunif(-5,5)
  sd.vis ~ dunif(0,5)
  tau.vis <- pow(sd.vis, -2)
  
  rho ~ dunif(-1,1)
  var.eta <- tau.v/(1.-pow(rho,2))
}
