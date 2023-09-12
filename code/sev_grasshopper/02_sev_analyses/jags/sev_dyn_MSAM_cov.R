model{
  #MSAM for Sevilleta Grasshoppers
  #Ana Miller-ter Kuile
  #September 11, 2023
  
  # This model is a multi-site, dynamic (multi-year) multi-species model 
  # for grasshoppers in the Sevilleta LTER
  # This script considers first year abundance and then any following year
  # abundance as dependent on the previous year abundance state 
  
  #This model also includes a covariance matrix among grasshopper
  # abundances, such that we are aiming to incorporate some estimate 
  # of species interactions, including competition, predation(?), 
  #facilitation, etc.
  
  # Right now this model has these aspects:
  ## 1. Multi-site, multi-species, multi-year
  ## 2. Species-level abundance parameter (lambda) varies by year
  ## 3. No covariates for the biological process (abundance)
  ## 4. Covariates on detection including season and life history
  
  
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
          a1.Rep[reprod[k]] #reproduction effect, not species dependent
        
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
    sp.lambda[k] <- ilogit(sp.llambda[k])
    
    
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
  
  #Prior for covariance matrix, omega
  # parameterize precision matrices with Wishart distributions
  #wishart uses R, which is a n.species x n.species matrix and 
  #degrees of freedom (second parameter) >= n.species
  #R is "data" you provide - will likely use cov() in R to create these data
  omega[1:n.species,1:n.species] ~ dwish(R[1:n.species,1:n.species], n.species)
  
  #ask Kiona: what is this sigma??? took from Shelby's code
  #This may be important if not converging with just omega -
  #look @ the JAGS user manual for more info on the other dmnorm distribution
  #and options
  #Sigma[1:n.species,1:n.species] <- inverse(omega[1:n.species,1:n.species])
  
  
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
  #categorical covariate of reproductive strategy:
  for(r in 2:3){
  a1.Rep[r] ~ dnorm(0, 1E-2)
  }
  
  a1.Rep[1] <- 0
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
      # NOTE: I double-checked this equation for bray-curtis against
      # the frequnty used 1 - ((2C[i,j])/(S[i] + s[j])) and they
      # get the same ansewr. The version I was using before from
      # the partitioning paper also gives the same ansewr and could be
      # good if we want to partition bray into components
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