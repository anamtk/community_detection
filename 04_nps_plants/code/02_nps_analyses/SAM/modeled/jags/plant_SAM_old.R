model{
  
  #this is a SAM model for the PFNP plant dataset
  
  #attributes of the model:
  # - beta distribution for Jaccard dissimilarlity input data
  # - incorporate uncertainty from the detection model
  # - SAM components to look at effects of VPD and PPT
  # - hierarchical spatial effects of quadrat within transect within plot
  # - various other options to code this RE structure, depending
  ### on the distribution of quadrats across plots and adding more data
  
  for(i in 1:n.data){
    
    #-------------------------------------## 
    # Likelihood ###
    #-------------------------------------##
    
    #diss is proportional, so beta distribution works here
    diss[i] ~ dbeta(alpha[i], beta[i])
    
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
    #VPD and PPT
    logit(mu[i]) <- b0.transect[Transect.ID[i]] +
      b[1]*AntVPD[i] +
      b[2]*AntPPT[i]
    #maybe also check the interaction - but I've been finding
    #that it overfits the data when I do this. 
    
    #-------------------------------------## 
    # SAM summing ###
    #-------------------------------------##
    
    #summing the antecedent values
    AntVPD[i] <- sum(VPDTemp[i,]) #summing across the total number of antecedent seasons
    AntPPT[i] <- sum(PPTTemp[i,]) #summing across the total num of antecedent seasons
    #Generating each year's weight to sum above
    for(t in 1:n.lag){ #number of time steps we're going back in the past
      VPDTemp[i,t] <- VPD[i,t]*wA[t] 
      PPTTemp[i,t] <- PPT[i,t]*wB[t]
      #missing data - you may not have any so you can check
      VPD[i,t] ~ dnorm(mu.vpd, tau.vpd)
      PPT[i,t] ~ dnorm(mu.ppt, tau.ppt)
    }
    
    #-------------------------------------## 
    # Goodness of fit parameters ###
    #-------------------------------------##
    # 
    # #replicated data
    diss.rep[i] ~ dbeta(alpha[i], beta[i])
    # 
    # #residuals - is this still right?
    resid[i] <- diss[i] - mu[i]
    
  }
  
  #-------------------------------------## 
  # Priors ###
  #-------------------------------------##
  
  # ANTECEDENT CLIMATE PRIORS
  #Sum of the weights for VPD lag
  sumA <- sum(deltaA[]) #all the VPD weights
  #Sum of the weights for PPT
  sumB <- sum(deltaB[]) #all the PPT weights
  
  #Employing "delta trick" to give vector of weights dirichlet priors
  #this is doing the dirichlet in two steps 
  #see Ogle et al. 2015 SAM model paper in Ecology Letters
  for(t in 1:n.lag){ #for the total number of lags
    #the weights for vpd - getting the weights to sum to 1
    wA[t] <- deltaA[t]/sumA
    #and follow a relatively uninformative gamma prior
    deltaA[t] ~ dgamma(1,1)
    
    #the weights for ppt - getting the weights to sum to 1
    wB[t] <- deltaB[t]/sumB
    #and follow a relatively uninformative gamma prior
    deltaB[t] ~ dgamma(1,1)
  }

  #BETA PRIORS
  #HIERARCHICAL STRUCTURE PRIORS
  
  #double check how many quadrats you have per plot - 
  # i would make "plot" your higher hierarchy in this if there
  #are ~5+ quadrats in each plot
  #if there aren't that many per plot (I don't remember how many
  #there are), probably having just the quadrat random effect is
  #fine. If we add more data (e.g. all quads on a transect),
  #then we will need to add in a transect hierarcy
  
  #OPTION 1:
  #hierarchical centering of quads on plots on b0
  # for(q in 1:n.quads){
  #   b0.quad[q] ~ dnorm(b0.plot[Plot.ID[q]], tau.quad)
  # }
  # 
  # for(p in 1:n.plots){
  #   b0.plot[p] ~ dnorm(b0, tau.plot)
  # }
  
  #OPTION 2: (if we add more quads - i added all quads on
  #transects with repeat measurements on at least one quadrat)
  # #hierarchical centering of quads on transects on plots on b0
  # for(q in 1:n.quads){
  #   b0.quad[q] ~ dnorm(b0.transect[Transect.ID[q]], tau.quad)
  # }

  for(t in 1:n.transects){
    b0.transect[t] ~ dnorm(b0, tau.transect)
  }
  
  # for(t in 1:n.transects){
  #   b0.transect[t] ~ dnorm(b0.plot[Plot.ID[t]], tau.transect)
  # }

  #check how many plots are in the dataset - if it isn't that many,
  #then we can remove this hierarchy (~ <5)
  # for(s in 1:n.plots){
  #   b0.plot[s] ~ dnorm(b0, tau.plot)
  # }

  #OPTION 3: (if low sample size/plot)
  #if not using plot RE
  # for(t in 1:n.quads){
  #    b0.quad[t] ~ dnorm(b0, tau.quad)
  # }
  
  b0 ~ dnorm(0, 1E-2)
  
  #for low # of levels, from Gellman paper - define sigma
  # as uniform and then precision in relation to this sigma
  #sig.quad ~ dunif(0, 10)
  sig.transect ~ dunif(0, 10)
  #sig.plot ~ dunif(0, 10)
  
  #tau.quad <- 1/pow(sig.quad,2)
  tau.transect <- 1/pow(sig.transect,2)
  #tau.plot <- 1/pow(sig.plot,2)
  
  for(i in 1:2){
    b[i] ~ dnorm(0, 1E-2)
  }
  
  #PRior for overall process error
  # sig.process ~ dunif(0, 10)
  # var.process <- pow(sig.process, 2)
  
  var.process ~ dunif(0, min(diff[]))
  
  #MISSING DATA PRIORS
  mu.vpd ~ dunif(-10, 10)
  sig.vpd ~ dunif(0, 20)
  tau.vpd <- pow(sig.vpd, -2)
  mu.ppt ~ dunif(-10, 10)
  sig.ppt ~ dunif(0, 20)
  tau.ppt <- pow(sig.ppt, -2)

  
}