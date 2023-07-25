model{
  #DYNAMIC MSAM FOR KELP FISH
  #Ana Miller-ter Kuile
  #July 25, 2023
  
  # This script is a multi-site, dynamic (multi-year) multi-species model 
  # for fish in the kelp forest in the SB LTER. 
  # This script considers first year abundance and then any following year
  # abundance as dependent on the number of individuals 
  # surviving from the previous year and 
  # recruiting more individuals from the previous year
  
  # Right now, I'm modeling transect-level data, but I would really like to 
  # be able to add a random hierarchy of transects within a reef site - not
  # sure where to put this random effect into the model - might be a good
  # convo for a future meeting with Kiona
  
  # Right now this model has these aspects:
  ## 1. Multi-site, multi-species, multi-year, abundance data
  ## 2. Survival and Recruitment vary across years
  ## 3. No covariates for initial abundance, survival, or recruitment
  ## 4. Covariates on detection including average size of the fish species 
  ### and visibility during the survey period
  
  
  #Biological process model:
  #First year
  for(k in 1:n.species){
    for(i in 1:n.transects){
      #some initial latent abundance based on a psi term 
      N[k,i,n.start[i]] ~ dpois(lambda0[k])
      
      #for all subsequent years:
      for(t in (n.start[i]+1):n.end[i]){
        #latent abundance is now based on 
        #survival and recruitment from previous year
        N[k,i,t] <- S[k,i,t] + 
                    R[k,i,t]
        
        #individuals that survive are based on
        #survival probability for that species, phi,
        #and the number of "trials" (individuals) 
        #in the previous timepoint
        #Should phi be yearly?? 
        S[k,i,t] ~ dbin(phi[k,t], N[k,i,t-1])
        
        #recruited individuals are based on combination
        #of reproduction and immigration??
        #gamma is this combined "reproduction + immigration"
        #term
        #Should gamma be yearly?? 
        R[k,i,t] ~ dpois(lambda[k,t])
        
        #consider adding - some kind of site-level
        #random effect - so that these things are linked
        #at the site level?
      } #year loop for after first year
  
  #Observation process model:
  for(t in n.start[i]:n.end[i]){
    for(r in 1:n.rep[i,t]){
      # Observation model
      logit(p[k,i,t,r]) <- a0[k] + #species-level intercept
        a1.Vis*vis[i,t,r] + #visibility effect, not species dependent
        a2.Size*size[k] #size effect, not species dependent
      
      #abundance is binomial based on detection probability
      #and total true abundance at the site
      y[k,i,t,r] ~ dbin(p[k,i,t,r], N[k,i,t])
    
      } #reps detection loop
    } #years detection loop
  } #transects loop
  
  #SPECIES-LEVEL PRIORS
    #species level priors for initial abundance,
    #survival probability, and colonization/recruitment
    #number are centered around community-level priors
    #for each of these variables
    #initial abundance:
    llambda0[k] ~ dnorm(mu.llambda0, tau.llambda0)
    lambda0[k] <- ilogit(llambda0[k])
    
    #survival and recruitment
    for(t in 2:n.years){
      #survival
      lphi[k,t] ~ dnorm(mu.lphi, tau.lphi)
      phi[k,t] <- ilogit(lphi[k,t])
      #recruitment
      llambda[k,t] ~ dnorm(mu.llambda, tau.llambda)
      lambda[k,t] <- ilogit(llambda[k,t])
    }

    #Detection intercept and slopes
    la0[k] ~ dnorm(mu.a0, tau.a0)
    a0[k] <- ilogit(la0[k])
    
    
  } #species loop
  
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
  
  #Community-level hyperpriors
  #All species-level priors are centered around hyperpriors for 
  # the community for that variaable
  
  #initial abundance
  lambda0.mean ~ dbeta(1,1)
  mu.llambda0 <- logit(lambda0.mean)
  sig.llambda0 ~ dunif(0, 10)
  tau.llambda0 <- pow(sig.llambda0, -2)
  
  #survival
  phi.mean ~ dbeta(1,1)
  mu.lphi <- logit(phi.mean)
  sig.lphi ~ dunif(0, 10)
  tau.lphi <- pow(sig.lphi, -2)
  
  #recruitment
  lambda.mean ~ dbeta(1,1)
  mu.llambda <- logit(lambda.mean)
  sig.llambda ~ dunif(0, 10)
  tau.llambda <- pow(sig.llambda, -2)
  
  #Detection intercept
  a0.mean ~ dbeta(1,1)
  mu.a0 <- logit(a0.mean)
  tau.a0 <- pow(sig.a0, -2)
  sig.a0 ~ dunif(0, 10)
  
  #covariate means
  a1.Vis ~ dnorm(0, 1E-2)
  a2.Size ~ dnorm(0, 1E-2)
  
  #PRIORS FOR IMPUTING MISSING DATA
  #Priors for mean and tau of missing covariates in the model
  mu.missingvis ~ dunif(-10, 10)
  sig.missingvis ~ dunif(0, 20)
  tau.missingvis <- pow(sig.missingvis, -2)
  
  
}