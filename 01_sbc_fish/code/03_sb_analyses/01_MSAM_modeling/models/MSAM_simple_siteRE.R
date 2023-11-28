model{

  for(k in 1:n.species){
    for(i in 1:n.transects){
      for(t in n.start[i]:n.end[i]){
        
        #biological process model
        N[k,i,t] ~ dpois(lambda[k,i,t])
        
        log(lambda[k,i,t]) <- b0.species[k] + 
          eps.site[Site.ID[i], k] +
          eps.year[Year.ID[t], k]
        
        for(r in 1:n.rep[i,t]){ #for the number of surveys on each transect in each year
          # Observation model
          logit(p[k,i,t,r]) <- a0[k] + #species-level intercept
            a1.Vis*vis[i,t,r] + #visibility effect, not species dependent
            a2.Size*size[k] #size effect, not species dependent
            
            #abundance is binomial based on detection probability
            #and total true abundance at the site
            y[k,i,t,r] ~ dbin(p[k,i,t,r], N[k,i,t])
            
            y.rep[k,i,t,r] ~ dbin(p[k,i,t,r], N[k,i,t])
            
        }
        
      }
      
    }
    
    #SPECIES-LEVEL PRIORS:
    #Detection intercept and slopes
    a0[k] ~ dnorm(mu.a0,tau.a0)
    
    #"baseline" detection at vis = 0 and size = 0 
    #on standardized scale
    p0[k] <- ilogit(a0[k])
    
    #species-level intercept - see note below about 
    #how to incorprate the site and year effects into
    #this or the community level for this
    b0.species[k] ~ dnorm(mu.b0species, tau.b0species)

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
  
  #Species random effects:
  for(s in 1:n.species){
    #non-identifiable random effect
    gamma.site[s] ~ dnorm(0, tau.gamma.site)
    #Identifiable random effect (monitor this):
    gamma.site.star[s] <- gamma.site[s] + 
      ave.eps.site[s] - 
      ave.gamma.site -
      ave.ave.eps.site
  }
  
  #Mean species random effects:
  ave.gamma.site <- mean(gamma.site[])
  #mean overall site random effect:
  ave.ave.eps.site <- mean(ave.eps.site[])
  
  
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
  
  #Species random effects:
  for(s in 1:n.species){
    #non-identifiable random effect
    gamma.year[s] ~ dnorm(0, tau.gamma.year)
    #Identifiable random effect (monitor this):
    gamma.year.star[s] <- gamma.year[s] + 
      ave.eps.year[s] - 
      ave.gamma.year -
      ave.ave.eps.year
  }
  
  #Mean species random effects:
  ave.gamma.year <- mean(gamma.year[])
  #mean overall year random effect:
  ave.ave.eps.year <- mean(ave.eps.year[])
  
  
  #OUTSTANDING QUESTION
  #what to do about overall b0 when it is indexed by species
  #and then hierarchically centered around the community value?
  # also how to incorporate both year and site effects into it?
  
  #code from lecture has b0 as:
  #b0 ~ dnorm(0, 1E-6)
  #identifiable intercept
  #b0.star <- b0 + ave.gamma + ave.ave.eps
  
  #initial occupancy
  #Detection intercept
  mu.a0 ~ dnorm(0, 0.001)
  tau.a0 <- pow(sig.a0, -2)
  sig.a0 ~ dunif(0, 50)
  
  #species-level abundance
  mu.b0species ~ dnorm(0, 0.001)
  tau.b0species <- pow(sig.b0species, -2)
  sig.b0species ~ dunif(0,50)
  
  #covariate effect priors
  a1.Vis ~ dnorm(0, 0.001)
  a2.Size ~ dnorm(0, 0.001)
  
  #missing data 
  #SOme data for visibility are missing, so we're imputing them
  for(i in 1:n.transects){
    for(t in n.start[i]:n.end[i]){
      for(r in 1:n.rep[i,t]){
        #missing data in the visibility column
        vis[i,t,r] ~ dnorm(mu.missingvis, tau.missingvis)
      }
    }
  }
  
  #PRIORS FOR IMPUTING MISSING DATA
  #Priors for mean and tau of missing covariates in the model
  mu.missingvis ~ dunif(-10, 10)
  sig.missingvis ~ dunif(0, 20)
  tau.missingvis <- pow(sig.missingvis, -2)
 
  
  #BRAY CURTIS DERIVED QUANTIIES
  #lots of ways to calculate this, but I did this way
  #IF WE WANT TO partition components of Bray:
  for(i in 1:n.transects){
    for(t in (n.start[i]+1):n.end[i]){
      for(k in 1:n.species){
        # num individuals in both time periods per species
        a[k,i,t] <- min(N[k,i,t-1], N[k,i,t])
        # num individuals only in first time point
        b[k,i,t] <- N[k,i,t-1] - a[k,i,t]
        # num individuals only in second time point
        c[k,i,t] <- N[k,i,t] - a[k,i,t]
      }
      #for all years 2 onward:
      #total number of shared individuals across time periods
      A[i,t] <- sum(a[,i,t])
      #total number of individuals in only first time period
      B[i,t] <- sum(b[,i,t])
      #total number of individuals in only second time period
      C[i,t] <- sum(c[,i,t])
      
      #total bray-curtis (B+C)/(2A+B+C)
      num[i,t] <- B[i,t] + C[i,t]
      denom1[i,t] <- 2*A[i,t]+B[i,t]+C[i,t]
      #if all values are zero - this just keeps the eqn. from
      #dividing by zero
      denom[i,t] <- ifelse(denom1[i,t]==0,1, denom1[i,t])
      bray[i,t] <- num[i,t]/denom[i,t]
      
      #how much is dissimilarity shaped by
      # individuals of one species being replaced by individuals
      #of another species?
      num.bal[i,t] <- min(B[i,t], C[i,t])
      denom.bal1[i,t] <- A[i,t] + num.bal[i,t]
      denom.bal[i,t] <- ifelse(denom.bal1[i,t] == 0,1, denom.bal1[i,t])
      bray_balanced[i,t] <- num.bal[i,t]/denom.bal[i,t]
      
      #how much is dissimilarity shaped by
      # individuals that are lost without substitution?
      bray_gradient[i,t] <- bray[i,t] - bray_balanced[i,t]
    }
  }
  
  

  
}