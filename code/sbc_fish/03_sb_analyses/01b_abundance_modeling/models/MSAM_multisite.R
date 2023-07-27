model{
  
  
  
  for(k in 1:n.species){
    for(i in 1:n.transects){
      for(t in n.start[i]:n.end[i]){
      #latent abundance
      N[k,i,t] ~ dpois(lambda[k,t])
      
      for(r in 1:n.rep[i,t]){
        # Observation model
        logit(p[k,i,t,r]) <- a0[k] + #species-level intercept
          a1.Vis*vis[i,t,r] + #visibility effect, not species dependent
          a2.Size*size[k] #size effect, not species dependent
        
        #abundance is binomial based on detection probability
        #and total true abundance at the site
        y[k,i,t,r] ~ dbin(p[k,i,t,r], N[k,i,t])
      } #reps
    }# years
    }#transects
    
    #SPECIES-LEVEL PRIORS:
    for(t in 1:n.years){
      llambda[k,t] ~ dnorm(mu.llambda, tau.llambda)
      lambda[k,t] <- ilogit(llambda[k,t])
    }
    
    #Detection intercept and slopes
    la0[k] ~ dnorm(mu.a0, tau.a0)
    a0[k] <- ilogit(la0[k])
    
    
  } #species
  
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
  
  #Community-level hyperpriors
  #All species-level priors are centered around hyperpriors for 
  # the community for that variaable
  
  #initial abundance
  lambda.mean ~ dbeta(1,1)
  mu.llambda <- logit(lambda.mean)
  sig.llambda ~ dunif(0, 10)
  tau.llambda <- pow(sig.llambda, -2)
  
  #Detection intercept
  a0.mean ~ dbeta(1,1)
  mu.a0 <- logit(a0.mean)
  tau.a0 <- pow(sig.a0, -2)
  sig.a0 ~ dunif(0, 50)
  
  #covariate means
  a1.Vis ~ dnorm(0, 1E-2)
  a2.Size ~ dnorm(0, 1E-2)
  
  #PRIORS FOR IMPUTING MISSING DATA
  #Priors for mean and tau of missing covariates in the model
  mu.missingvis ~ dunif(-10, 10)
  sig.missingvis ~ dunif(0, 20)
  tau.missingvis <- pow(sig.missingvis, -2)
  
  
  # #DERIVED PARAMETERS##
  
  #BRAY-CURTIS DISSIMILARITY
  for(i in 1:n.transects){
    for(t in (n.start[i]+1):n.end[i]){
      for(k in 1:n.species){
        # num individuals in both time periods per species
        a[k,i,t] <- min(N[k,i,t-1], N[k,i,t])
        # num individuals only in first time point
        b[k,i,t] <- N[k,i,t-1] - (min(N[k,i,t-1], N[k,i,t]))
        # num individuals only in second time point
        c[k,i,t] <- N[k,i,t] - (min(N[k,i,t-1], N[k,i,t]))
      }
      #for all years 2 onward:
      #total number of shared individuals across time periods
      A[i,t] <- sum(a[,i,t])
      #total number of individuals in only first time period
      B[i,t] <- sum(b[,i,t])
      #total number of individuals in only second time period
      C[i,t] <- sum(c[,i,t])
      
      #total bray-curtis (B+C)/(2A+B+C)
      #0 means the two sites have the same composition 
      #(that is they share all the same num of individuals), 
      #and 1 means the two sites do not share same
      #number of individuals.
      bray[i,t] <- (B[i,t] + C[i,t])/(2*A[i,t]+B[i,t]+C[i,t])
      # 
      # balanced variation in abundance:
      # #how much is dissimilarity shaped by
      # # individuals of one species being replaced by individuals
      # #of another species?
      min[i,t] <- min(B[i,t], C[i,t])
      denom[i,t] <- A[i,t] + min[i,t]
      bray_b[i,t] <- min[i,t]/denom[i,t]
      #bray_balanced[i,t] <- min[i,t]/(A[i,t] + min[i,t])
      #bray_balanced[i,t] <- min(B[i,t],C[i,t])/(A[i,t] + min(B[i,t],C[i,t]))
      # 
      # abundance gradient:
      # #how much is dissimilarity shaped by
      # # individuals that are lost without substitution?
      # bray_gradient[i,t] <- bray[i,t] - bray_balanced[i,t]
    }
  }
  
  
  
}