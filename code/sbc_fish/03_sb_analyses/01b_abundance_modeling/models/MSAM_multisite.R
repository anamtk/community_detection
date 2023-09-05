model{
  #MSAM for Kelp Fish
  #Ana Miller-ter Kuile
  #August 31, 2023
  
  # This model is a multi-site, dynamic (multi-year) multi-species model 
  # for fish in the kelp forest in the SB LTER. 
  # This script considers first year abundance and then any following year
  # abundance as dependent on the previous year abundance state 
  
  #This model also includes a covariance matrix among fish abundances, such
  #that we are aiming to incorporate some estimate of species interactions,
  #including competition, predation, facilitation, etc.
  
  # Right now, I'm modeling transect-level data, but I would really like to 
  # be able to add a random hierarchy of transects within a reef site - not
  # sure where to put this random effect into the model - might be a good
  # convo for a future meeting with Kiona
  
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
    la0[k] ~ dnorm(mu.a0, tau.a0)
    a0[k] <- ilogit(la0[k])
    
    #to make recursive portion (below), we need to know a mean species-level
    #lambda to populate the first year of species-year lambdas
    sp.llambda[k] ~ dnorm(mu.llambda, tau.llambda) #centered around community mean
    sp.lambda[k] <- ilogit(sp.mu.llambda[k])
    
    
  } #species
  
  #Species-year level priors
  #Create recursive (year to year dependence) and covariance between species
  #by centering all yearly lambdas by species around a mean value for each species
  #or around the lambda of the year before.
  
  #DOUBLE CHECK WITH KIONA - THAT these should be on lambda (ilogit) scale, not the llambda scale
  #year 1 lambda for each species
  lambda[1:n.species,1] ~ dmnorm(sp.lambda[1:n.species],omega[1:n.species,1:n.species])
  
  #years 2+ lambda for each species
  for(t in 2:n.years){
    lambda[1:n.species,t] ~ dmnorm(lambda[1:n.species, t-1],omega[1:n.species,1:n.species])
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
      #makign if else so it's never quite 0 and we can do division below
      A[i,t] <- sum(a[,i,t])
      #total number of individuals in only first time period
      #makign if else so it's never quite 0 and we can do division below
      B[i,t] <- ifelse(sum(b[,i,t]) == 0, 0.00001, sum(b[,i,t]))
      #total number of individuals in only second time period
      #makign if else so it's never quite 0 and we can do division below
      C[i,t] <- ifelse(sum(c[,i,t]) == 0, 0.00001, sum(c[,i,t]))
      
      #total bray-curtis (B+C)/(2A+B+C)
      #0 means the two sites have the same composition 
      #(that is they share all the same num of individuals), 
      #and 1 means the two sites do not share same
      #number of individuals.
      Denom1[i,t] <- (2*A[i,t]+B[i,t]+C[i,t])
      Denom[i,t] <- ifelse(Denom1[i,t] == 0, 1,Denom1[i,t])
      bray[i,t] <- (B[i,t] + C[i,t])/(Denom[i,t])
      # 
      # balanced variation in abundance:
      # #how much is dissimilarity shaped by
      # # individuals of one species being replaced by individuals
      # #of another species?
      # min[i,t] <- min(B[i,t], C[i,t])
      # denom[i,t] <- A[i,t] + min[i,t]
      # bray_b[i,t] <- min[i,t]/denom[i,t]
      #bray_balanced[i,t] <- min[i,t]/(A[i,t] + min[i,t])
      bray_balanced[i,t] <- min(B[i,t],C[i,t])/(A[i,t] + min(B[i,t],C[i,t]))
      # 
      # abundance gradient:
      # #how much is dissimilarity shaped by
      # # individuals that are lost without substitution?
      bray_gradient[i,t] <- bray[i,t] - bray_balanced[i,t]
    }
  }
  
  
  
}