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
      for(t in 1:n.yr[i]){
        
        z[k,i,t] ~ dbern(psi[k,i,t])
        #occupancy is dependent currently on species, quadrat, and time period
        #for that quadrat
        
        #Detection model:
        #with no covariates for detection:
        y[k,i,t] ~ dbin(p[k] * z[k,i,t], n.rep[i,t])
        #y is 1 if that species was ever observed in any repeat survey for that quad-year
        #y is 0 if that species was never observed in any repeat survey for that quad-year

        
        #if we add covariates for detection, which I'll ask Megan abaout:
        #or add in cover classes here?
        # for(r in 1:n.rep[i,t]){ #for the number of surveys on each transect in each year
        #   # Observation model
        #   logit(p[k,i,t,r]) <- a0[k]  #species-level intercept
        #   
        #   #presence is binomial based on detection probability conditioned
        #     #on true abundance and the total number of reps as trials 
        #   
        #   #in this case, y is dependent on repeat survey because of the
        #   # dependence of p on survey period
        #     y[k,i,t,r] ~ dbin(p[k,i,t,r] * z[k,i,t], n.rep[i,t])
        # } #reps
        
        #SPECIES-LEVEL PRIORS:
        
        #occupancy species-year-level priors
        lpsi[k,i,t] ~ dnorm(mu.lpsi, tau.lpsi)
        psi[k,i,t] <- ilogit(lpsi[k,i,t])
          
      }# years
    } #quads
    
    
    #SPECIES-LEVEL PRIORS:
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