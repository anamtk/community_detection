model{
  
  for(i in 1:n.data){
    
    #-------------------------------------## 
    # Likelihood ###
    #-------------------------------------##
    
    #bray is proportional, so beta distribution works here
    bray[i] ~ dbeta(alpha[i], beta[i])
    
    #var.process is scalar but could be made dependent on site/other variables
    #phi incorporates mu (mean estimate), var.estimate (which is from
    #"data" on standard deviation (squared) from the original detection 
    #correction model) 
    #and var.process is something we're tryign to estimate,
    #basically, the rest of the variation not accounted for
    phi[i] <- (((1-mu[i])*mu[i])/(var.estimate[i] + var.process))-1
    
    #alpha and beta are based on mu and phi values
    #sometimes these values send alpha and beta outside
    #the domain, so we have extra code below to get them to
    #stay where they belong
    alphaX[i] <- mu[i] * phi[i]
    betaX[i] <- (1 - mu[i]) * phi[i]
    
    #here is where we get alpha and beta to stay in their
    #domain
    alpha[i] <- max(0.01, alphaX[i])
    beta[i] <- max(0.01, betaX[i])
    
    #to get a good estimate of a prior for var.process, we
    #track the difference between these two values for each
    #data point
    diff[i] <- (1-mu[i])*mu[i] - var.estimate[i]
    
    #Regression of mu, which is dependent on antecedent
    #temperature, precipitation and npp
    logit(mu[i]) <- b0 +
      b[1]*AntTemp[i] +
      b[2]*AntPPT[i] +
      b[3]*AntNPP[i] 
    
    #-------------------------------------## 
    # SAM summing ###
    #-------------------------------------##
    
    #summing the antecedent values
    AntTemp[i] <- sum(TempTemp[i,]) #summing across the total number of antecedent months
    AntPPT[i] <- sum(TempPPT[i,]) #summing across the total num of antecedent months
    AntNPP[i] <- sum(TempNPP[i,]) #summing across total num of antecedent years
    
    #generating each month's weight to sum above
    for(t in 1:n.templag){ #number of time steps we're going back in the past
      TempTemp[i,t] <- Temp[i,t]*wA[t] 
      TempPPT[i,t] <- PPT[i,t]*wB[t]
      
      #missing data
      Temp[i,t] ~ dnorm(mu.temp, tau.temp)
      PPT[i,t] ~ dnorm(mu.ppt, tau.ppt)
    }
    
    #Generating each lag's weight to sum above
    for(t in 1:n.npplag){ #number of time steps we're going back in the past
      TempNPP[i,t] <- NPP[i,t]*wC[t] 
      
      #missing data
      NPP[i,t] ~ dnorm(mu.npp, tau.npp)
    }
    
    #-------------------------------------## 
    # Goodness of fit parameters ###
    #-------------------------------------##
    # 
    # #replicated data
    # beta.rep[i] ~ dbeta(alpha[i], beta[i])
    # 
    # #residuals - is this still right?
    # resid[i] <- bray[i] - mu[i]
    
  }
  
  #-------------------------------------## 
  # Priors ###
  #-------------------------------------##
  
  # ANTECEDENT CLIMATE PRIORS
  #Sum of the weights for temp and ppt lags
  sumA <- sum(deltaA[]) #all the temp weights
  
  sumB <- sum(deltaB[]) #all the ppt weights
  #Employing "delta trick" to give vector of weights dirichlet priors
  #this is doing the dirichlet in two steps 
  #see Ogle et al. 2015 SAM model paper in Ecology Letters
  for(t in 1:n.templag){ #for the total number of lags
    #the weights for temp - getting the weights to sum to 1
    wA[t] <- deltaA[t]/sumA
    #and follow a relatively uninformative gamma prior
    deltaA[t] ~ dgamma(1,1)
    
    #the weights for ppt - getting the weights to sum to 1
    wB[t] <- deltaB[t]/sumB
    #and follow a relatively uninformative gamma prior
    deltaB[t] ~ dgamma(1,1)
  }
  
  #sum of weights for the npp lag
  sumC <- sum(deltaC[])
  #Employing "delta trick" to give vector of weights dirichlet priors
  #this is doing the dirichlet in two steps 
  #see Ogle et al. 2015 SAM model paper in Ecology Letters
  for(t in 1:n.npplag){ #for the total number of lags

    #the weights for npp
    wC[t] <- deltaC[t]/sumC
    deltaC[t] ~ dgamma(1,1)
  }
  
  #BETA PRIORS
  #only three sites in this dataset, so not adding in 
  #any site-level hierarchies for now - but could??
  b0 ~ dnorm(0, 1E-2)
  
  for(i in 1:3){
    b[i] ~ dnorm(0, 1E-2)
  }
  
  #PRior for overall process error
  # sig.process ~ dunif(0, 10)
  # var.process <- pow(sig.process, 2)
  
  var.process ~ dunif(0, min(diff[]))
  
  #MISSING DATA PRIORS
  mu.temp ~ dunif(-10, 10)
  sig.temp ~ dunif(0, 20)
  tau.temp <- pow(sig.temp, -2)
  mu.ppt ~ dunif(-10, 10)
  sig.ppt ~ dunif(0, 20)
  tau.ppt <- pow(sig.ppt, -2)
  mu.npp ~ dunif(-10, 10)
  sig.npp ~ dunif(0, 20)
  tau.npp <- pow(sig.npp, -2)
  
}