model{
  #MSOM for Petrified forest Plants
  #Ana Miller-ter Kuile
  #September 11, 2023
  
  # This model is a multi-site, dynamic (multi-year) multi-species model 
  # for plants in the Petrified Forest NP datas
  # This script considers first year presence and then any following year
  # presence as dependent on the previous year presence state 
  
  #This model also includes a covariance matrix among plant abundances, such
  #that we are aiming to incorporate some estimate of species interactions,
  #including competition, facilitation, etc.
  
  # Right now this model has these aspects:
  ## 1. Multi-site, multi-species, multi-year
  ## 2. Species-level presence parameter (psi) varies by year
  ## 3. No covariates for the biological process (abundance)
  ## 4. Covariates on detection including what??? observer???
  
  
  for(k in 1:n.species){ #loop through all species in the community
    for(i in 1:n.quads){ #loop through each transect
      for(t in 1:n.yr[i]){ #loop through the start to end year for each transect
        
        #Biological process model:
        #latent presence
        
        z[k,i,t] ~ dbern(psi[k,i,t])
        #presence probability is dependent on species and year
        
        #Detection model:
        #with no covariates for detection:
        y[k,i,t] ~ dbin(p[k] * z[k,i,t], n.rep[i,t])
        
        #if we add covariates for detection, which I'll ask Megan abaout:
        #or add in cover classes here?
        # for(r in 1:n.rep[i,t]){ #for the number of surveys on each transect in each year
        #   # Observation model
        #   logit(p[k,i,t,r]) <- a0[k]  #species-level intercept
        #   
        #   #presence is binomial based on detection probability conditioned
        #     #on true abundance and the total number of reps as trials 
        #   
        #     y[k,i,t,r] ~ dbin(p[k,i,t,r] * z[k,i,t], n.rep[i,t])
        # } #reps
        
      }# years
    }#transects
    
    #SPECIES-LEVEL PRIORS:
    
    # #Detection intercept
    # la0[k] ~ dnorm(mu.a0, tau.a0)
    # a0[k] <- ilogit(la0[k])
    #Detection hierarchy
    lp[k] ~ dnorm(mu.lp, tau.lp)
    p[k] <- ilogit(lp[k])
    
    #to make recursive portion (below), we need to know a mean species-level
    #psi to populate the first year of species-year psis
    sp.lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi) #centered around community mean
    sp.psi[k] <- ilogit(sp.lpsi[k])
    
    
  } #species
  
  
  for(i in 1:n.quads){
    #Species-year level priors
    #Create recursive (year to year dependence) and covariance between species
    #by centering all yearly psis by species around a mean value for each species
    #or around the psis of the year before.
    
    #year 1 psi for each species
    psi[1:n.species,i,1] ~ dmnorm(sp.psi[1:n.species],omega[1:n.species,1:n.species])
    
    #years 2+ lambda for each species
    for(t in 2:n.yr[i]){
      psi[1:n.species,i,t]  ~ dmnorm(psi[1:n.species,i, t-1],omega[1:n.species,1:n.species])
    }
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
  psi.mean ~ dbeta(1,1)
  mu.lpsi <- logit(psi.mean)
  sig.lpsi ~ dunif(0, 10)
  tau.lpsi <- pow(sig.lpsi, -2)
  
  # #Detection intercept
  # a0.mean ~ dbeta(1,1)
  # mu.a0 <- logit(a0.mean)
  # tau.a0 <- pow(sig.a0, -2)
  # sig.a0 ~ dunif(0, 50)
  
  p.mean ~ dbeta(1, 1)
  mu.lp <- logit(p.mean)
  sig.lp ~ dunif(0, 5)
  tau.lp <- pow(sig.lp, -2)
  
  #covariate means
  # a1.Effort ~ dnorm(0, 1E-2)
  # a2.Size ~ dnorm(0, 1E-2)
  
  # #DERIVED PARAMETERS##
  
  #JACCARD DISSIMILARITY
  
  for(i in 1:n.quads){
    for(t in 2:n.yr[i]){
      for(k in 1:n.species){
        #is species k lost in site i between t and t+1?
        #if lost, value of a will be 1
        a[k,i,t] <- (z[k,i,t-1] == 1)*(z[k,i,t] == 0)
        #is species k gained in site i between t and t+1
        #if gained, value of b will be 1
        b[k,i,t] <- (z[k,i,t-1]== 0)*(z[k,i,t] == 1)
        #is species k shared in site i between t and t+1
        #if shared, value of c will be 1
        c[k,i,t] <- (z[k,i,t-1]==1)*(z[k,i,t]==1)
      }
      #for all years 2 onward:
      #total number of species lost
      A[i,t] <- sum(a[,i,t])
      #total number of species gained
      B[i,t] <- sum(b[,i,t])
      #total number of species shared
      C[i,t] <- sum(c[,i,t])
      
      # #total turnover is (A+B)/(A+B+C)
      tot_turnover[i,t] <- (A[i, t] + B[i, t])/
        (A[i, t] + B[i, t] + C[i, t])
      # #gain is B/(A+B+C)
      gain[i,t] <- (B[i, t])/
        (A[i, t] + B[i, t] + C[i, t])
      # #loss is A/(A+B+C)
      loss[i,t] <- (A[i, t])/
        (A[i, t] + B[i, t] + C[i, t])
      #
      # #Jaccard beta diversity is shared/total, so C/A+B+C
      jaccard[i,t] <- (C[i, t])/
        (A[i, t] + B[i, t] + C[i, t])
    }
  }

  
}