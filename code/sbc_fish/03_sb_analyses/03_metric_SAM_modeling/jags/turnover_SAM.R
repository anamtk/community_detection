model{
  
  for(i in 1:n.data){
    
    #-------------------------------------## 
    # Likelihood ###
    #-------------------------------------##
    
      turn[i] ~ dbeta(alpha[i], beta[i])
      
      alpha[i] <- mu[i] * phi
      beta[i]  <- (1-mu[i]) * phi
    #add in: transect within site crossed with yearly
    #random effects
      logit(mu[i]) <- b0.transect[Transect.ID[i]] +
        #b0.year[Year.ID[i]] +
        b[1]*AntKelp[i] +
        b[2]*AntTemp[i]
      

      
      #-------------------------------------## 
      # SAM summing ###
      #-------------------------------------##
      
      #summing the antecedent values
      AntKelp[i] <- sum(KelpTemp[i,]) #summing across the total number of antecedent years
      AntTemp[i] <- sum(TempTemp[i,]) #summing across the total num of antecedent months
      
      
      #Generating each year's weight to sum above
      for(t in 1:n.kelplag){ #number of time steps we're going back in the past
        KelpTemp[i,t] <- Kelp[i,t]*wA[t] 
      
        #missing data
        Kelp[i,t] ~ dnorm(mu.kelp, tau.kelp)
      }
        
      #generating each month's weight to sum above
      for(t in 1:n.templag){ #number of time steps we're going back in the past
        TempTemp[i,t] <- Temp[i,t]*wB[t] 
        
        #missing data
        Temp[i,t] ~ dnorm(mu.temp, tau.temp)
      }
      
      #-------------------------------------## 
      # Goodness of fit parameters ###
      #-------------------------------------##
      
      #replicated data
      turn.rep[i] ~ dbeta(alpha[i], beta[i])
      
      #residuals - is this still right?
      resid[i] <- turn[i] - mu[i]
 
  }
  
  #-------------------------------------## 
  # Priors ###
  #-------------------------------------##
  
  # ANTECEDENT CLIMATE PRIORS
  #Sum of the weights for kelp lag
  sumA <- sum(deltaA[]) #all the kelp weights
  #Employing "delta trick" to give vector of weights dirichlet priors
  #this is doing the dirichlet in two steps 
  #see Ogle et al. 2015 SAM model paper in Ecology Letters
  for(t in 1:n.kelplag){ #for the total number of lags
    #the weights for kelp - getting the weights to sum to 1
    wA[t] <- deltaA[t]/sumA
    #and follow a relatively uninformative gamma prior
    deltaA[t] ~ dgamma(1,1)
  }
  
  #Sum of the weights for temp lag
  sumB <- sum(deltaB[]) #all the temp weights
  #Employing "delta trick" to give vector of weights dirichlet priors
  #this is doing the dirichlet in two steps 
  #see Ogle et al. 2015 SAM model paper in Ecology Letters
  for(t in 1:n.templag){ #for the total number of lags
    #the weights for kelp - getting the weights to sum to 1
    wB[t] <- deltaB[t]/sumB
    #and follow a relatively uninformative gamma prior
    deltaB[t] ~ dgamma(1,1)
  }
  
  #BETA PRIORS
  #HIERARCHICAL STRUCTURE PRIORS
  #hierarchical centering of transects on sites on b0
   for(t in 1:n.transects){
     b0.transect[t] ~ dnorm(b0.site[Site.ID[t]], tau.transect)
   }
   
   for(s in 1:n.sites){
     b0.site[s] ~ dnorm(b0, tau.site)
   }
  
  #if not using site RE
  # for(t in 1:n.transects){
  #    b0.transect[t] ~ dnorm(b0, tau.transect)
  # }
  
  b0 ~ dnorm(0, 1E-2)
  
  #SUM TO ZERO for years
  #for every year but the last one:
  # for(y in 2:(n.years)){
  #   b0.year[y] ~ dnorm( 0, tau.year)
  # }
  #set the last year to be the -sum of all other years so the 
  # overall fo all year levels == 0
  #b0.year[n.years+1] <- -sum(b0.year[2:(n.years)])
  
  #for low # of levels, from Gellman paper - define sigma
  # as uniform and then precision in relation to this sigma
  sig.transect ~ dunif(0, 10)
  sig.site ~ dunif(0, 10)
  #sig.year ~ dunif(0, 10)
  
  tau.transect <- 1/pow(sig.transect,2)
  tau.site <- 1/pow(sig.site,2)
  #tau.year <- 1/pow(sig.year, 2)
  
  for(i in 1:2){
    b[i] ~ dnorm(0, 1E-2)
  }
  
  #overall phi value
  phi ~ dgamma(.1,.1)
  
  #MISSING DATA PRIORS
  mu.kelp ~ dunif(-10, 10)
  sig.kelp ~ dunif(0, 20)
  tau.kelp <- pow(sig.kelp, -2)
  mu.temp ~ dunif(-10, 10)
  sig.temp ~ dunif(0, 20)
  tau.temp <- pow(sig.temp, -2)
  
}