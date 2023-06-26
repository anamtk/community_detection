model {

  for (k in 1:(n.species+n.zeros)) {
    for (i in 1:n.years) {
      #logit detection is linked to phi, which ties occ. to det. probability
      logit(psi[i,k]) <- phi[k]
      #mean occupancy for species i in site k is related to occupancy times 
      # community inclusion for species i
      mu.psi[i,k] <- psi[i,k]*w[k]
      #true occupancy for species follows a bernoilli distribution
      #with occupancy probabiliy mu.psi
      z[i,k] ~ dbern(mu.psi[i,k]) #occupancy
      
      #DETECTION MODEL
      logit(p[i,k]) <- eta[k]+a1.Vis[k]*vis[i]	#detection
      mu.p[i,k] <- p[i,k]*z[i,k]
      
      #observed data
      y[i,k] ~ dbin(mu.p[i,k], n.reps)

    }
    
    #original membership in the community
    w[k] ~ dbern(omega)				#membership in community
    #covariate for occupancy is phi - which links occupancy and detection
    phi[k] ~ dnorm(beta, tau.u)
    a1.Vis[k] ~ dnorm(mu.vis, tau.vis)
    
    #mu.eta and eta links occupancy to detection
    mu.eta[k] <- alpha + (rho*sigma.v/sigma.u)*(phi[k] - beta)
    #links detection to occupancy
    eta[k] ~ dnorm(mu.eta[k], var.eta)
    
  }
  
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

