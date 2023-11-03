model{

  for(k in 1:n.species){
    for(i in 1:n.transects){
      for(t in n.start[i]:n.end[i]){
        
        #biological process model
        # z[k,i,t] ~ dbern(omega[k,t])
        # N[k,i,t] ~ dpois(lambda[k,t]*z[k,i,t])
        # 
        N[k,i,t] ~ dpois(lambda[k,t])
        
        for(r in 1:n.rep[i,t]){ #for the number of surveys on each transect in each year
          # Observation model
          logit(p[k,i,t,r]) <- a0[k] + #species-level intercept
            a1.Effort*effort[i,t,r] + #visibility effect, not species dependent
            a2.Size*size[k] #size effect, not species dependent
            
            #abundance is binomial based on detection probability
            #and total true abundance at the site
            y[k,i,t,r] ~ dbin(p[k,i,t,r], N[k,i,t])
            
            #y.rep[k,i,t,r] ~ dbin(p[k,i,t,r], N[k,i,t])
            
        }
        
      }
      
    }
    
    #SPECIES-LEVEL PRIORS:
    #Detection intercept and slopes
    a0[k] ~ dnorm(mu.a0, tau.a0)
    #"baseline" detection at covariates = 0
    p0[k] <- ilogit(a0[k])
    
    for(t in 1:n.years){
    llambda[k,t] ~ dnorm(mu.llambda, tau.llambda) #centered around community mean
    lambda[k,t] <- exp(llambda[k,t])
    }

  }
  
  #Community-level hyperpriors
  #initial abundance
  mu.llambda ~ dnorm(0, 0.00001)
  sig.llambda ~ dunif(0, 10)
  tau.llambda <- pow(sig.llambda, -2)
  
  #initial occupancy
  #Detection intercept
  mu.a0 ~ dnorm(0, 0.001)
  tau.a0 <- pow(sig.a0, -2)
  sig.a0 ~ dunif(0, 50)
  
  #covariate means
  a1.Effort ~ dnorm(0, 0.001)
  a2.Size ~ dnorm(0, 0.001)
  
  #missing data 
  #SOme data for effort are missing, so we're imputing them
  for(i in 1:n.transects){
    for(t in n.start[i]:n.end[i]){
      for(r in 1:n.rep[i,t]){
        #missing data in the effort column
        effort[i,t,r] ~ dnorm(mu.missingeffort, tau.missingeffort)
      }
    }
  }
  
  #PRIORS FOR IMPUTING MISSING DATA
  #Priors for mean and tau of missing covariates in the model
  mu.missingeffort ~ dunif(-10, 10)
  sig.missingeffort ~ dunif(0, 20)
  tau.missingeffort <- pow(sig.missingeffort, -2)

  
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