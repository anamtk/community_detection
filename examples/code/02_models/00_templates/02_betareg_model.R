model{
  
  #this is general model code for a simple beta regression
  #with antecedent terms
  
  #the code as written is for modeled data, we have provided
  #the code for using empirical data (without uncertainty 
  #estimates modeled) in commented code where approrpiate in this
  #model
  
  for(i in 1:n.data){
    
    #-------------------------------------## 
    # Likelihood ###
    #-------------------------------------##
    
    #diss is proportional, so beta distribution works here
    #diss can be any dissimilarlity measure that ranges
    #in domain [0,1]
    diss[i] ~ dbeta(alpha[i], beta[i])
      
    #var.process is scalar but could be made dependent on site/other variables
    #phi incorporates mu (mean estimate), var.estimate (which is from
    #"data" on standard deviation (squared) from the original detection 
    #correction model) 
    #and var.process is something we're tryign to estimate,
    #basically, the rest of the variation not accounted for
    phi[i] <- (((1-mu[i])*mu[i])/(var.estimate[i] + var.process))-1

    #for empirical estiamtes of dissimilarity, phi is simplified:
    #phi[i] <- (((1-mu[i])*mu[i])/var.process)-1
    
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
    
    #for empirical data, you will not need to compute diff[i]

    #Regression of mu, which is dependent on antecedent
    #variable AntVar
      logit(mu[i]) <- b0.site[Site.ID[i]] +
        b[1]*AntVar[i] 
      
      #-------------------------------------## 
      # SAM summing ###
      #-------------------------------------##
      
      #summing the antecedent values
      AntVar[i] <- sum(VarTemp[i,]) #summing across the total number of antecedent steps

      #Generating each year's weight to sum above
      for(t in 1:n.lag){ #number of time steps we're going back in the past
        VarTemp[i,t] <- Var[i,t]*wA[t] 
      
        #missing data
        Var[i,t] ~ dnorm(mu.var, tau.var)
      }
        
      #-------------------------------------## 
      # Goodness of fit parameters ###
      #-------------------------------------##
      # 
      # #replicated data
      diss.rep[i] ~ dbeta(alpha[i], beta[i])
      # 
      # #residuals
      resid[i] <- diss[i] - mu[i]
 
  }
  
  #-------------------------------------## 
  # Priors ###
  #-------------------------------------##
  
  # ANTECEDENT CLIMATE PRIORS
  #Sum of the weights for var lag
  sumA <- sum(deltaA[]) #all the  weights
  #Employing "delta trick" to give vector of weights dirichlet priors
  #this is doing the dirichlet in two steps 
  #see Ogle et al. 2015 SAM model paper in Ecology Letters
  for(t in 1:n.lag){ #for the total number of lags
    #the weights for var - getting the weights to sum to 1
    wA[t] <- deltaA[t]/sumA
    #and follow a relatively uninformative gamma prior
    deltaA[t] ~ dgamma(1,1)
    
    #get the cumulative weights through time
    cumm.wt[t] <- sum(wA[1:t])
  }
  
  #BETA PRIORS
  #HIERARCHICAL STRUCTURE PRIORS
  #hierarchical centering of sites on b0

   for(s in 1:n.sites){
     b0.site[s] ~ dnorm(b0, tau.site)
   }
  

  b0 ~ dnorm(0, 1E-2)
  
  #for low # of levels, from Gellman paper - define sigma
  # as uniform and then precision in relation to this sigma
  sig.site ~ dunif(0, 10)
  tau.site <- 1/pow(sig.site,2)
  
  for(i in 1:1){
    b[i] ~ dnorm(0, 1E-2)
  }
  
  #PRior for overall process error
  var.process ~ dunif(0, min(diff[]))
  
  #for empirical estimates, var.process prior is simply:
  #var.process ~ dunif(0, 100)

  #MISSING DATA PRIORS
  mu.var ~ dunif(-10, 10)
  sig.var ~ dunif(0, 20)
  tau.var <- pow(sig.var, -2)
  
}