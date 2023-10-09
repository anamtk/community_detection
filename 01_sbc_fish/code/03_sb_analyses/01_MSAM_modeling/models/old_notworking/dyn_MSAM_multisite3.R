model{
 
  
  for(k in 1:n.species){
    for(i in 1:n.transects){
      for(t in (n.start[i]:n.end[i])){
        
        N[k,i,t] ~ dpois(lambda[k,t])
        
        for(r in 1:n.rep[i,t]){ #for the number of surveys on each transect in each year
          # Observation model
          logit(p[k,i,t,r]) <- a0[k] + #species-level intercept
            a1.Vis*vis[i,t,r] + #visibility effect, not species dependent
            a2.Size*size[k] #size effect, not species dependent
          
          #abundance is binomial based on detection probability
          #and total true abundance at the site
          y[k,i,t,r] ~ dbin(p[k,i,t,r], N[k,i,t])
          #y[k,i,t,r] ~ dbin(p[k,i,t,r], ymax2[k,i,t])
        } #reps
        

      }
    }
    
    #prior for lambda
    for(t in 1:22){
    lambda[k,t] ~ dnorm(mu.llambda, tau.llambda)
    }
    
    #priors for species-level values
    la0[k] ~ dnorm(mu.a0, tau.a0)
    a0[k] <- ilogit(la0[k])
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
  sig.a0 ~ dunif(0, 10)
  
  #covariate means
  a1.Vis ~ dnorm(0, 1E-2)
  a2.Size ~ dnorm(0, 1E-2)
  
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
  
  
}