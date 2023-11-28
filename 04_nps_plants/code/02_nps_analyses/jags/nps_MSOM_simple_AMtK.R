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
        
        logit(psi[k,i,t]) <- b0.species[k] +
          eps.site[Site.ID[i], k] +
          eps.year[Year.ID[t], k]
        
        #Detection model:
        #with no covariates for detection:
        # y[k,i,t] ~ dbin(p[k] * z[k,i,t], n.rep[i,t])
        #y is 1 if that species was ever observed in any repeat survey for that quad-year
        #y is 0 if that species was never observed in any repeat survey for that quad-year

        for(r in 1:n.rep[i,t]){ #for the number of surveys on each transect in each year
          # Observation model
          logit(p[k,i,t,r]) <- a0[k] + #species-level intercept
            a1.Cover*cover[k,i,t,r] + #proxy for abundance
            a2.LifeGroup[lifegroup[k]] #also potentially a proxy for abundance,
          #this is a categorical combination of the lifegroup and duration values

          #presence is binomial based on detection probability conditioned
            #on true abundance and the total number of reps as trials

          #in this case, y is dependent on repeat survey because of the
          # dependence of p on survey period
            y[k,i,t,r] ~ dbern(p[k,i,t,r] * z[k,i,t])
            
            y.rep[k,i,t,r] ~ dbern(p[k,i,t,r]*z[k,i,t])
            
            #MISSING DATA for cover imputation
            #cover[k,i,t,r] ~ dnorm(mu.missingcover, tau.missingcover)
        } #reps
        

          
      }# years
    } #quads
    
    
    #SPECIES-LEVEL PRIORS:
    #abundance model (biological process)
    #non-identifiable:
    b0.species[k] ~ dnorm(mu.b0species, tau.b0species)
    #identifiable species-level intercept:
    #track this for convergence
    b0.star[k] <- b0.species[k] + ave.eps.site[k] + ave.eps.year[k]
    
    #if we add in detection covariates
    # #Detection intercept
    a0[k] ~ dnorm(mu.a0, tau.a0)
    
    #"baseline" detection at all covariates == 0
    p0[k] <- ilogit(a0[k])
    
  }
  
  #SITE W/IN SPECIES RANDOM EFFECTS
  #sites nested within species (sites sum to zero w/in each species)
  for(s in 1:n.species){
    for(p in 1:n.sites){
      #non-identifiable random effect
      eps.site[p,s] ~ dnorm(0, tau.eps.site)
      #identifiable site random effect (monitor this)
      eps.site.star[p,s] <- eps.site[p,s] - ave.eps.site[s]
    }
    #mean site level random effects within each species
    ave.eps.site[s] <- mean(eps.site[,s])
  }
  
  #YEARS W/IN SPECIES RANDOM EFFECTS
  #sites nested within species (sites sum to zero w/in each species)
  for(s in 1:n.species){
    for(p in 1:n.years){
      #non-identifiable random effect
      eps.year[p,s] ~ dnorm(0, tau.eps.year)
      #identifiable year random effect (monitor this)
      eps.year.star[p,s] <- eps.year[p,s] - ave.eps.year[s]
    }
    #mean year level random effects within each species
    ave.eps.year[s] <- mean(eps.year[,s])
  }
  
  
  #Community-level hyperpriors
  #All species-level priors are centered around hyperpriors for 
  # the community for that variaable
  
  #If we add detection covariates:
  # #Detection intercept
  mu.a0 ~ dnorm(0, 0.001)
  tau.a0 <- pow(sig.a0, -2)
  sig.a0 ~ dunif(0, 50)
  
  #covariate means
  a1.Cover ~ dnorm(0, 1E-3)
  
  #species-level abundance
  mu.b0species ~ dnorm(0, 0.001)
  tau.b0species <- pow(sig.b0species, -2)
  sig.b0species ~ dunif(0,50)
  
  #site and year variances
  sig.eps.site ~ dunif(0, 10)
  tau.eps.site <- pow(sig.eps.site, -2)
  sig.eps.year ~ dunif(0, 10)
  tau.eps.year <- pow(sig.eps.year, -2)
  
  #categorical covariate of lifegroup
  #this is cell-referenced - so each level other than the first
  #gets a prior that looks like any other a prior, but the 
  #first level is set to 0, thus, all other level values are
  #in comparison to this baseline value. 
  #we set the baseline to be the group with the most observations because
  #this helps the model statistically 
  for(g in 2:n.groups){ #number of life groups
    a2.LifeGroup[g] ~ dnorm(0, 1E-3)
  }
  
  # SL: error with running when this is set to 0
  # added condition in regression for logit(p) so this is 0 when lifegroup = 1
  # but will have a value other than 0 in the posterior results
  a2.LifeGroup[1] <- 0
  
  #PRIORS FOR IMPUTING MISSING DATA
  #Priors for mean and tau of missing covariates in the model
  #mu.missingcover ~ dunif(-10, 10)
  #sig.missingcover ~ dunif(0, 20)
  #tau.missingcover <- pow(sig.missingcover, -2)
  
  
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