model{
  #DYNAMIC MSOM FOR KONZA BIRDS
  #Ana Miller-ter Kuile
  #July 27, 2023
  
  # This script is a multi-site, dynamic (multi-year) multi-species model 
  # for birds in the Konza LTER
  # This script considers first year occupancy and then any following year
  # occupancy as dependent on the previous year occupancy state and 
  # persistence and colonization probabilities.
  
  # Right now, I'm modeling transect-level data
  
  # Right now this model has these aspects:
  ## 1. Multi-site, multi-species, multi-year
  ## 2. Persistence and colonization vary across years
  ## 3. No covariates for initial occupancy, persistence, or colonization
  ## 4. Covariates on detection including effort (time of survey) 
  
  
  for(k in 1:n.species){
    for(i in 1:n.transects){
      
      #BIOLOGICAL MODEL
      #year 1 true occupancy
      #year one occupancy probability is
      #based on initial occupancy probability, psi1 for species k in
      # site i
      #true occupancy is bernoulli around first-year
      #occupancy probability
      z[k,i,n.start[i]] ~ dbern(psi1[k])
      
      #year 2+ true occupancy
      for(t in (n.start[i]+1):n.end[i]){
        #occupancy dependent on persistence, phi
        # and colonization, gamma
        #persistence only applies if previous year z = 1
        # colonization only applies if previous year z = 0
        #right now, allowing persistence and colonization
        #to remain the same for a speciesxyear combo
        #but could add in site covariates if we wanted down the line
        z[k,i,t] ~ dbern(z[k,i, t-1]*phi[k,t] +
                           (1- z[k,i, t-1])*gamma[k,t])
      }#Year 2+ loop
      
      #OBSERVATION MODEL
      
      for(t in n.start[i]:n.end[i]){
        for(r in 1:n.rep[i,t]){
          # Observation model
          logit(p[k,i,t,r]) <- a0[k] + #species-level intercept
            a1.Effort*effort[i,t,r]  #effort effect
          
          #detection is bionmiaml based on detection probability times
          # true occupancy, with the number of trials being the 
          # nubmer of visits to that site that year
          y[k,i,t,r] ~ dbin(p[k,i,t,r] * z[k,i,t], n.rep[i,t])
          
        } #reps detection loop
      } #years detection model loop
    }#transects likelihood loop
    
    #Species-level PRIORS
    # Species level priors for occupancy, persistence, colonization,
    # and detection are centered around community-level priors for each
    # of these variables
    
    #occupancy
    
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi) 
    psi1[k] <- ilogit(lpsi[k])
    
    for(t in 2:n.years){
      #persistence
      lphi[k,t] ~ dnorm(mu.lphi, tau.lphi)
      phi[k,t] <- ilogit(lphi[k,t])
      
      #colonization
      lgamma[k,t] ~ dnorm(mu.lgamma, tau.lgamma)
      gamma[k,t] <- ilogit(lgamma[k,t])
      
    }
    
    #detection intercept and slopes
    la0[k] ~ dnorm(mu.a0, tau.a0)
    a0[k] <- ilogit(la0[k])
    
  } #species loop
  
  #missing data 
  #SOme data for visibility are missing, so we're imputing them
  for(i in 1:n.transects){
    for(t in n.start[i]:n.end[i]){
      for(r in 1:n.rep[i,t]){
        #missing data in the visibility column
        effort[i,t,r] ~ dnorm(mu.missingeff, tau.missingeff)
      }
    }
  }
  
  #Community-level hyperpriors
  #All species-level priors are centered around hyperpriors for 
  # the community for that variaable
  
  #initial occupancy
  psi.mean ~ dbeta(1,1)
  mu.lpsi <- logit(psi.mean)
  sd.lpsi ~ dunif(0, 10)
  tau.lpsi <- pow(sd.lpsi, -2)
  
  #persistence
  phi.mean ~ dbeta(1,1)
  mu.lphi <- logit(phi.mean)
  sd.lphi ~ dunif(0, 10)
  tau.lphi <- pow(sd.lphi, -2)
  
  #colonization
  gamma.mean ~ dbeta(1,1)
  mu.lgamma <- logit(gamma.mean)
  sd.lgamma ~ dunif(0, 10)
  tau.lgamma <- pow(sd.lgamma, -2)
  
  #Detection intercept
  a0.mean ~ dbeta(1,1)
  mu.a0 <- logit(a0.mean)
  tau.a0 <- pow(sig.a0, -2)
  sig.a0 ~ dunif(0, 10)
  
  #covariate means
  a1.Effort ~ dnorm(0, 1E-2)
  
  #PRIORS FOR IMPUTING MISSING DATA
  #Priors for mean and tau of missing covariates in the model
  mu.missingeff ~ dunif(-10, 10)
  sig.missingeff ~ dunif(0, 20)
  tau.missingeff <- pow(sig.missingeff, -2)
  
  
  # #DERIVED PARAMETERS##
  for(i in 1:n.transects){
    for(t in (n.start[i]+1):n.end[i]){
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