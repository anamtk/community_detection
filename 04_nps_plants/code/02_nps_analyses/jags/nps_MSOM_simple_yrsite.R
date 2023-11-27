model{
  
  for(k in 1:n.species){ #loop through all species in the community
    for(t in 1:n.years){
        for(i in 1:quads[t]){
          
          z[k,t,i] ~ dbern(psi[k,t, i])
          
          logit(psi[k,t,i]) <- b0.species[k] +
            eps.site[Site.ID[t,i], k] +
            eps.year[Year.ID[t], k]
          
          for(r in 1:n.rep[t,i]){
            
            logit(p[k,t,i,r]) <- a0[k] +
              a1.Cover*cover[k] + #average cover for a species across all plots, proxy for "abundance/size"
              a2.LifeGroup[lifegroup[k]]
            
            y[k,t,i,r] ~ dbern(p[k,t,i,r]*z[k,t,i])
            
            y.rep[k,t,i,r] ~ dbern(p[k,t,i,r]*z[k,t,i])
            
          }
        }

    }
    
    #if we add in detection covariates
    # #Detection intercept
    a0[k] ~ dnorm(mu.a0, tau.a0)
    
    #"baseline" detection at all covariates == 0
    p0[k] <- ilogit(a0[k])
    
    #abundance model (biological process)
    #non-identifiable:
    b0.species[k] ~ dnorm(mu.b0species, tau.b0species)
    #identifiable species-level intercept:
    #track this for convergence
    b0.star[k] <- b0.species[k] + ave.eps.site[k] + ave.eps.year[k]
    
  }
  
  #SITE W/IN SPECIES RANDOM EFFECTS
  #sites nested within species (sites sum to zero w/in each species)
  for(s in 1:n.species){
    for(p in 1:n.sites){
      #non-identifiable random effect
      eps.site[p,s] ~ dnorm(0, tau.eps.site)
      #identifiable site random effect (monitor this)
      eps.site.star[p,s] <- eps.site[p,s] - ave.eps.site[s]
    }
    #mean site level random effects within each species
    ave.eps.site[s] <- mean(eps.site[,s])
  }
  
  #YEARS W/IN SPECIES RANDOM EFFECTS
  #sites nested within species (sites sum to zero w/in each species)
  for(s in 1:n.species){
    for(p in 1:n.years){
      #non-identifiable random effect
      eps.year[p,s] ~ dnorm(0, tau.eps.year)
      #identifiable year random effect (monitor this)
      eps.year.star[p,s] <- eps.year[p,s] - ave.eps.year[s]
    }
    #mean year level random effects within each species
    ave.eps.year[s] <- mean(eps.year[,s])
  }
  

  #Community-level hyperpriors
  #All species-level priors are centered around hyperpriors for 
  # the community for that variaable
  
  #species-level abundance
  mu.b0species ~ dnorm(0, 0.001)
  tau.b0species <- pow(sig.b0species, -2)
  sig.b0species ~ dunif(0,50)
  
  #site and year variances
  sig.eps.site ~ dunif(0, 10)
  tau.eps.site <- pow(sig.eps.site, -2)
  sig.eps.year ~ dunif(0, 10)
  tau.eps.year <- pow(sig.eps.year, -2)

  #If we add detection covariates:
  # #Detection intercept
  mu.a0 ~ dnorm(0, 0.001)
  tau.a0 <- pow(sig.a0, -2)
  sig.a0 ~ dunif(0, 50)

  #covariate means
  a1.Cover ~ dnorm(0, 1E-3)

  #categorical covariate of lifegroup
  #this is cell-referenced - so each level other than the first
  #gets a prior that looks like any other a prior, but the 
  #first level is set to 0, thus, all other level values are
  #in comparison to this baseline value. 
  #we set the baseline to be the group with the most observations because
  #this helps the model statistically 
  for(g in 2:n.groups){ #number of life groups
    a2.LifeGroup[g] ~ dnorm(0, 1E-3)
    }

  # SL: error with running when this is set to 0
  # added condition in regression for logit(p) so this is 0 when lifegroup = 1
  # but will have a value other than 0 in the posterior results
  a2.LifeGroup[1] <- 0

  #PRIORS FOR IMPUTING MISSING DATA
  #Priors for mean and tau of missing covariates in the model
  #mu.missingcover ~ dunif(-10, 10)
  #sig.missingcover ~ dunif(0, 20)
  #tau.missingcover <- pow(sig.missingcover, -2)



}