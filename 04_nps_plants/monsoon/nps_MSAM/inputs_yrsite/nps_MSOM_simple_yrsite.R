model{
  
  for(k in 1:n.species){ #loop through all species in the community
    for(t in 1:n.years){
        for(i in 1:quads[t]){
          
          z[k,t,i] ~ dbern(psi[k,t])
          
          for(r in 1:n.rep[t,i]){
            
            logit(p[k,t,i,r]) <- a0[k] +
              a1.Cover*cover[k] + #average cover for a species across all plots, proxy for "abundance/size"
              a2.LifeGroup[lifegroup[k]]
            
            y[k,t,i,r] ~ dbern(p[k,t,i,r]*z[k,t,i])
            
          }
        }
      #yearly values for each species for psi
      lpsi[k,t] ~ dnorm(mu.lpsi, tau.lpsi)
      psi[k,t] <- ilogit(lpsi[k,t])
    }
    
    #if we add in detection covariates
    # #Detection intercept
    a0[k] ~ dnorm(mu.a0, tau.a0)
    
    #"baseline" detection at all covariates == 0
    p0[k] <- ilogit(a0[k])
  }

  #Community-level hyperpriors
  #All species-level priors are centered around hyperpriors for 
  # the community for that variaable
  
  #occupancy probability
  #psi.mean ~ dbeta(1,1)
  #mu.lpsi <- logit(psi.mean)
  sig.lpsi ~ dunif(0, 100)
  tau.lpsi <- pow(sig.lpsi, -2)

  mu.lpsi ~ dnorm(0, 0.0001)

  #Detection community means
  # p.mean ~ dbeta(1, 1)
  # mu.lp <- logit(p.mean)
  # sig.lp ~ dunif(0, 5)
  # tau.lp <- pow(sig.lp, -2)

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