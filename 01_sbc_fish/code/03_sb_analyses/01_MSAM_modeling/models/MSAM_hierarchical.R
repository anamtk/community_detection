model{
  #MSAM for Kelp Fish
  #Ana Miller-ter Kuile
  #August 31, 2023
  
  # This model is a multi-site, dynamic (multi-year) multi-species model 
  # for fish in the kelp forest in the SB LTER. 
  # This script considers first year abundance and then any following year
  # abundance as dependent on the previous year abundance state 
  
  # Right now this model has these aspects:
  ## 1. Multi-site, multi-species, multi-year
  ## 2. Species-level abundance parameter (lambda) varies by year
  ## 3. No covariates for the biological process (abundance)
  ## 4. Covariates on detection including average size of the fish species 
  ### and visibility during the survey period
  
  
  for(k in 1:n.species){ #loop through all species in the community
    for(i in 1:n.transects){ #loop through each transect
      for(t in n.start[i]:n.end[i]){ #loop through the start to end year for each transect
        
        #Biological process model
        #latent abundance
        N[k,i,t] ~ dpois(lambda[k,t])
        #abundance rate parameter, lambda, is dependent on species and year
        
        for(r in 1:n.rep[i,t]){ #for the number of surveys on each transect in each year
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
    
    #Detection intercept and slopes
    # KO: If you include additional hierarchies for species (e.g., species
    # nested in functional group or genus), I'd also include those additional 
    # hierarchies here.
    ## KO: Important: a0 is the intercept for in the logit-scale linear model,
    ## and thus should not constrain a0 to (0,1), with the following ilogit 
    ## coding is doing.
    #la0[k] ~ dnorm(mu.a0, tau.a0) ## Not appropriate
    #a0[k] <- ilogit(la0[k]) ## Not appropriate
    ## KO: This would be more appropriate:
    a0[k] ~ dnorm(mu.a0, tau.a0)
    ## KO: If you want to "look" at the predicted "base-line" probability of detection
    ## e.g., when vis = 0 and size = 0, then you could compute:
    p0[k] <- ilogit(a0[k])
    
    #to make recursive portion (below), we need to know a mean species-level
    #lambda to populate the first year of species-year lambdas
    sp.llambda[k] ~ dnorm(mu.llambda, tau.llambda) #centered around community mean
    ## KO: Important, lambda is the expected abundance associated with the Poisson
    ## model for latent abundance; we do not want to use ilogit here as this will 
    ## constrain sp.lambda to (0,1), which doesn't make sense for an abundance
    ## quantity.
    #sp.lambda[k] <- ilogit(sp.llambda[k]) ## Not appropriate
    ## KO: This is more appropriate:
    #sp.lambda[k] <- log(sp.llambda[k])
    sp.lambda[k] <- exp(sp.llambda[k]) #is this correct, actually?
    
    
  } #species
  
  #Species-year level priors
  #Create recursive (year to year dependence)
  #by centering all yearly lambdas by species around a mean value for each species
  #or around the lambda of the year before.
  
  #DOUBLE CHECK WITH KIONA - THAT these should be on lambda (ilogit) scale, not the llambda scale
  #also - what should be precision on these??? are they hierarchical??
  #year 1 lambda for each species
  for(k in 1:n.species){
    ## KO: we might "anchor" the initial lambda at sp.lambda, so that:
    ## lambda[k,1] = sp.lambda[k]
    ## KO: also, yes, I would model on the LOG scale and compute lambda; 
    ## If we give a normal prior to lambda, this could allow for lambda < 0.
    #lambda[k,1] ~ dnorm(sp.lambda[k], sp.tau.lambda[k]) ## Not appropriate?
    ## KO: More appropriate:
    # llambda[k,1] ~ dnorm(sp.llambda[k], sp.tau.lambda[k])
    # lambda[k,1] <- exp(llambda[k,1])
    ## KO: This is the version that would "anchor" lambda at t = 1:
    lambda[k,1] <- sp.lambda[k]
    llambda[k,1] <- sp.llambda[k]
    #years 2+ lambda for each species
    for(t in 2:n.years){
      ## KO: Again, would model lambda on log scale to obey it's domain:
      #lambda[k,t] ~ dnorm(lambda[k, t-1], sp.tau.lambda[k]) ## Not appropriate
      ## KO: More appropriate:
      llambda[k,t] ~ dnorm(llambda[k, t-1], sp.tau.lambda[k])
      lambda[k,t] <- exp(llambda[k,t])
    }
  }
  
  #Species level precision priors
  #hierarhical prior for sig.lambda:
  # folded t distribution with 2 degrees of freedom for standard deviation
  # for(k in 1:n.species){
  #   t.obs[k] ~ dt(0,D,2) # folded t distribution with 2 degrees of freedom for standard deviation
  #   #(sig:add a small value re: Doing Bayesian Analysis Book)
  #   ## KO: I don't see why we need to add the "small" value; I've never had to do this, and
  #   ## the "original" Gelman paper about modeling hierarchical variances doesn't do this.
  #   sp.sig.lambda[k] <- abs(t.obs[k]) + 0.001  # abs value, hierarchical folded t priors for sd
  #   #compute tau from sig
  #   sp.tau.lambda[k] <- pow(sp.sig.lambda[k],-2)# compute precision based on folded-t sd
  # }
  #from Kiona's meta-analysis papers
  ## KO: We want to use the above folded-t as a hierarchical prior for the 
  ## species-level standard deviations, so that the "population-level" parameter
  ## (D) is treated as unknown and given a prior (not fixed):
  ## D ~ dunif(0,)
  #D <- 1/(1*1)
  
  ## KO: Swap above with hierarchical folded-t prior below:
  for(k in 1:n.species){
    # folded t prior for species-level standard deviations, with 2 degrees of 
    # freedom and precision-type parameter, D:
    t.lambda[k] ~ dt(0,D,2) 
    # abs value to produce "folded" t:
    sp.sig.lambda[k] <- abs(t.lambda[k]) + 0.001
    #compute tau from sig
    sp.tau.lambda[k] <- pow(sp.sig.lambda[k],-2)
  }
  # Give prior to population-level scale parameter (A) and compute precision-
  # type parameter, D:
  ## KO: want to make sure U(0,10) is wide enough; check posterior for A,
  ## to see if it bumps up against upper limit (e.g., 10)
  E ~ dunif(0,10)
  D <- 1/(E*E)
  
  
  #Community-level hyperpriors
  #All species-level priors are centered around hyperpriors for 
  # the community for that variable
  
  #initial abundance
  ## KO: Why is lambda.mean being constrained to (0,1) with the beta prior?
  ## KO: mu.llambda is on a log-scale (e.g., log expected abundance), not logit
  ## KO: scale, right?
  # lambda.mean ~ dbeta(1,1)
  # mu.llambda <- logit(lambda.mean)
  sig.llambda ~ dunif(0, 10)
  tau.llambda <- pow(sig.llambda, -2)
  ## KO: The following would be more appropriate
  mu.llambda ~ dnorm(0,0.00001)
  lambda.mean <- exp(mu.llambda)
  
  
  #Detection intercept
  ## KO: Again, beta and logit transform lot appropriate here.
  # a0.mean ~ dbeta(1,1)
  # mu.a0 <- logit(a0.mean)
  ## KO: More appropriate:
  mu.a0 ~ dnorm(0,0.001)
  # Can compute population-level "baseline" probability of detection:
  mu.p0 <- ilogit(mu.a0)
  
  # KO: The priors below are appropriate, just need to check that posterior
  # for sig.a0 isn't bumping up against prior upper limit (e.g., 50)
  tau.a0 <- pow(sig.a0, -2)
  sig.a0 ~ dunif(0, 50)
  
  #covariate means
  ## KO: These priors for the covariate EFFECTS are fine, but might want to 
  ## specify a slightly smaller precision (e.g., 0.001)
  a1.Vis ~ dnorm(0, 0.001)
  a2.Size ~ dnorm(0, 0.001)
  
  #missing data 
  #SOme data for visibility are missing, so we're imputing them
  for(i in 1:n.transects){
    for(t in n.start[i]:n.end[i]){
      for(r in 1:n.rep[i,t]){
        #missing data in the visibility column
        ## KO: Seems fine, assuming vis can be positive or negative (which would
        ## be the case if continuous and standardized)
        vis[i,t,r] ~ dnorm(mu.missingvis, tau.missingvis)
      }
    }
  }
  
  #PRIORS FOR IMPUTING MISSING DATA
  #Priors for mean and tau of missing covariates in the model
  mu.missingvis ~ dunif(-10, 10)
  sig.missingvis ~ dunif(0, 20)
  tau.missingvis <- pow(sig.missingvis, -2)
  
  
  # #DERIVED PARAMETERS##
  
  #BRAY-CURTIS DISSIMILARITY
  
  # for(i in 1:n.transects){
  #   for(t in 1:n.start[i]+1:n.end[i]){
  #     for(k in 1:n.species){
  #       
  #       #From the vegan package, Bray-Curtis is given as 
  #       #(sumfromspeciesitoS(abs(x[i,t-1] - x[i,t])))/
  #       #(sumfromspeciesitoS(x[i,t-1] + x[i,t]))
  #       #where x is the number of individuals of species i in 
  #       #a given timepoint (t-1 or t), the two communities we're comparing
  #       
  #       #species-level numerator and denominator values
  #       sp.nums[k,i,t] <- abs(N[k,i,t-1] - N[k,i,t])
  #       sp.denoms[k,i,t] <- N[k,i,t-1] + N[k,i,t]
  #       
  #     }
  #     
  #     #the numerator sums all of the species-level numerators
  #     Num[i,t] <- sum(sp.nums[,i,t])
  #     #the denominator sums all of the species-level denominators
  #     Denom1[i,t] <- sum(sp.denoms[,i,t])
  #     #not sure this is important anymore, but this is just 
  #     #when the denominator is ==0, which would also make the 
  #     #numerator equal to zero. Set it to 1 so it won't break
  #     Denom[i,t] <- ifelse(Denom1[i,t] == 0, 1, Denom1[i,t])
  #     
  #     #then bray-curtis is just the numerator over the denominator
  #     bray[i,t] <- Num[i,t]/Denom[i,t]
  #     #partitioning difference:
  #     #Maybe for later - but could just be good to have one 
  #     #type of output per dataset as a proof of concept
  #     
  #     #NOTE: I double-checked this equation for bray-curtis against
  #     #the frequnty used 1 - ((2C[i,j])/(S[i] + s[j])) and they
  #     #get the same ansewr. The version I was using before from
  #     #the partitioning paper also gives the same ansewr and could be 
  #     #good if we want to partition bray into components
  #    
  #   }
  # }
  
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