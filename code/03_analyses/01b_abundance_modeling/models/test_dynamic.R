model{
  for(i in 1:n.sites){
    #BIOLOGICAL MODEL
    #First year abundance is based on year
    #1 lambda
    log(lambda[1]) <- llambda
    N[i,1] ~ dpois(lambda[1])
    
    llambda ~ dnorm(0, 1E-2)
    
    #in subsequent years
    for(t in 2:n.years){
      #lambda of that year is a combo of 
      #survival and coloniziation and lambda
      #from previous year 
      #look at "WHY" of this fanciness
      lambda[t] <- lambda[t-1]*phi[t] +
        lambda[t-1]*gamma[t]
      
      #surviving individuals based on 
      #survival prob phi, and N from previous
      #year
      S[i,t] ~ dbin(phi[t], N[i,t-1])
    
      #gains dependent on lambda of previous year?
      #times colonization/recruitment prob
      G[i,t] ~ dpois(lambda[t-1]*gamma[t])
      
      #number of total individuals is 
      #then these added together
      N[i,t] <- S[i,t] + G[i,t]
    } #2+years
    
    #OBSERVATION MODEL
    for(t in 1:n.years){
      for(r in 1:n.reps){
        #detection prob is dependnent on
        #some covariate at the replicate level
        logit(p[i,t,r]) <- b0 + b1*cov[i,t,r]
        
        #Data are "trials" of observed number
        #of individuals, based on detection error
        #and the total number of real individual
        #in the population
        y[i,t,r] ~ dbin(p[i,t,r], N[i,t])
        } #reps
    } #years
  } #sites
  
  #PRIORS:
  for(t in 1:n.years){
    lphi[t] ~ dnorm(mu.lphi, tau.lphi)
    phi[t] <- ilogit(lphi[t])
    
    lgamma[t] ~ dnorm(mu.lgamma, tau.lgamma)
    gamma[t] <- ilogit(lgamma[t])
  }

  b0 ~ dnorm(0, 1E-2)
  b1 ~ dnorm(0, 1E-2)
  
  #Hyperpriors
  #survival
  phi.mean ~ dbeta(1,1)
  mu.lphi <- logit(phi.mean)
  sig.lphi ~ dunif(0, 10)
  tau.lphi <- pow(sig.lphi, -2)
  
  #recruitment
  gamma.mean ~ dbeta(1,1)
  mu.lgamma <- logit(gamma.mean)
  sig.lgamma ~ dunif(0, 10)
  tau.lgamma <- pow(sig.lgamma, -2)
}