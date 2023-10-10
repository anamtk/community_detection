model{
  #MSOM for Petrified forest Plants
  #Ana Miller-ter Kuile
  #September 11, 2023
  
  # This model is a multi-site, dynamic (multi-year) multi-species model 
  # for plants in the Petrified Forest NP datas
  # This script considers first year presence and then any following year
  # presence as dependent on the previous year presence state, and values
  # estimating persistence and colonization
  
  # Right now this model has these aspects:
  ## 1. Multi-site, multi-species, multi-year
  ## 2. Species-level presence parameter (psi) varies by year
  ## 3. No covariates for the biological process (abundance)
  ## 4. Covariates on detection including what??? observer???
  
  
  for(k in 1:n.species){ #loop through all species in the community
    for(i in 1:n.quads){ #loop through each transect
      
      #Biological process model:
      
      #first year presence around species-level initial occupancy probability
      z[k,i,1] ~ dbern(psi1[k])
      
      #for all subsequent years
      for(t in 2:n.yr[i]){ #loop through the start to end year for each transect
        
        #Biological process model:
        #latent presence
        z[k,i,t] ~ dbern(psi[k,i,t])
        
        #psi is recursive based on extinction and colonization
        #it is dependent on species, site, and year
        #phi is the probbility a population persists if z[t-1] =1
        #phi changes by species and year
        #gamma is the probability of colonization if z[t-1] = 0
        #gamma changes by species and year
        
        #this first part will be = 0 if z=0
        psi[k,i,t] <- z[k,i,t-1]*phi[k,t-1] +
          #this second part will be =0 if z=1
          (1 - z[k,i,t-1])*gamma[k,t-1]
      }
      
      #detection is for all years - so we start a new loop
      for(t in 1:n.yr[i]){
        for(r in 1:n.rep[i,t]){
          #Detection model:
          #with no covariates for detection:
          y[k,i,t,r] ~ dbin(p[k] * z[k,i,t], n.rep[i,t])
          
        }
        
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
    } #quads
    
    
    #SPECIES-LEVEL PRIORS:
  
  #initial year occupancy species-level priors
  lpsi1[k] ~ dnorm(mu.lpsi, tau.lpsi)
  psi1[k] <- ilogit(lpsi[k])
  
  #persistence and colonization are currently
  #dependent on species and year
  for(t in 1:(n.yr[i]-1)){
    lphi[k,t] ~ dnorm(mu.lphi, tau.lphi)
    phi[k,t] <- ilogit(lphi[k,t])
    lgamma[k,t] ~ dnorm(mu.lgamma, tau.lgamma)
    gamma[k,t] <- ilogit(lgamma[k,t])
    
  } #persistence/colonization year loop
  
    #Detection hierarchy
    lp[k] ~ dnorm(mu.lp, tau.lp)
    p[k] <- ilogit(lp[k])
    
    #if we add in detection covariates
    # #Detection intercept
    # a0[k] ~ dnorm(mu.a0, tau.a0)
    
    }
    
  
  #Community-level hyperpriors
  #All species-level priors are centered around hyperpriors for 
  # the community for that variaable
  
  #initial abundance
  psi.mean ~ dbeta(1,1)
  mu.lpsi <- logit(psi.mean)
  sig.lpsi ~ dunif(0, 10)
  tau.lpsi <- pow(sig.lpsi, -2)
  
  #persistence and colonization
  phi.mean ~ dbeta(1,1)
  mu.lphi <- logit(phi.mean)
  sig.lphi ~ dunif(0, 10)
  tau.lphi <- pow(sig.lphi, -2)
  
  gamma.mean ~ dbeta(1,1)
  mu.lgamma <- logit(gamma.mean)
  sig.lgamma ~ dunif(0, 10)
  tau.lgamma <- pow(sig.lgamma, -2)
  
  #Detection community means
  p.mean ~ dbeta(1, 1)
  mu.lp <- logit(p.mean)
  sig.lp ~ dunif(0, 5)
  tau.lp <- pow(sig.lp, -2)
  
  #If we add detection covariates:
  # #Detection intercept
  # a0.mean ~ dbeta(1,1)
  # mu.a0 <- logit(a0.mean)
  # tau.a0 <- pow(sig.a0, -2)
  # sig.a0 ~ dunif(0, 50)
  
  #covariate means
  # a1.Effort ~ dnorm(0, 1E-3)
  # a2.Size ~ dnorm(0, 1E-3)
  
  # #DERIVED PARAMETERS##
  
  #JACCARD DISSIMILARITY
  
  for(i in 1:n.quads){
    for(t in 2:n.yr[i]){
      for(k in 1:n.species){
        #is species k lost in site i between t and t+1?
        #if lost, value of a will be 1
        b[k,i,t] <- (z[k,i,t-1] == 1)*(z[k,i,t] == 0)
        #is species k gained in site i between t and t+1
        #if gained, value of b will be 1
        c[k,i,t] <- (z[k,i,t-1]== 0)*(z[k,i,t] == 1)
        #is species k shared in site i between t and t+1
        #if shared, value of c will be 1
        a[k,i,t] <- (z[k,i,t-1]==1)*(z[k,i,t]==1)
      }
      #for all years 2 onward:
      #total number of species lost
      B[i,t] <- sum(a[,i,t])
      #total number of species gained
      C[i,t] <- sum(b[,i,t])
      #total number of species shared
      A[i,t] <- sum(c[,i,t])
      
      # #total turnover is (A+B)/(A+B+C)
      tot_turnover[i,t] <- (B[i, t] + C[i, t])/
        (A[i, t] + B[i, t] + C[i, t])
      # #gain is B/(A+B+C)
      gain[i,t] <- (C[i, t])/
        (A[i, t] + B[i, t] + C[i, t])
      # #loss is A/(A+B+C)
      loss[i,t] <- (B[i, t])/
        (A[i, t] + B[i, t] + C[i, t])
      #
      # #Jaccard beta diversity is shared/total, so C/A+B+C
      jaccard[i,t] <- (A[i, t])/
        (A[i, t] + B[i, t] + C[i, t])
    }
  }
  
  
}