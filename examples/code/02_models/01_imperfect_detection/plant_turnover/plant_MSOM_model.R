model{

  #Multi-species occupancy model for NPS Plant dataset
  #September 2023

  #This is a multi-species occupancy model for understory plants in the 
  #NPS Southern Colorado Plateau Inventory and Monitoring Network
  #Data are here:
  #https://doi.org/10.57830/2300890

  #Aspects of the model:
  ## 1. Due to the structure of the data, this model differs from other
  ## examples in this paper because different sites were surveyed in different years
  ## so that indexing needed to be years and then sites within years vs. sites and then years
  ## 2. There are no covariates on occupancy (biological process model)
  ## 3. Occupancy depends on a species-level intercept and random effects of 
  ## site within species and year within species which are identifiable
  ## due to post-sweeping 
  ## 4. Detection is based on a cover class (converted to median % for 
  ## each class and life form (e.g., perennial vs. annual; grass vs. forb).
  ## 5. Unlike the other models in this paper, because of the coding of 
  ## sites within years, community dissimilarity is calcuated on posterior
  ## samples from this model in R, versus coded directly into this model.
  
  for(k in 1:n.species){ #loop through all species in the community
    for(t in 1:n.years){ #loop for all years in dataset
        for(i in 1:quads[t]){ #loop for all sites in the dataset for that year

          #Biological Process model#
          #latent occupancy is dependent on occupancy probability, psi, for
          #that species in that year at that site
          z[k,t,i] ~ dbern(psi[k,t, i])

          #occupancy probability is dependent on species and 
          #random effects of site within species and year
          #within species
          logit(psi[k,t,i]) <- b0.species[k] +
            eps.site[Site.ID[t,i], k] +
            eps.year[Year.ID[t], k]

          #Observation Model:#
          for(r in 1:n.rep[t,i]){ #the number of repeated surveys in year t for quadrat i

            #detection probability, p:
            logit(p[k,t,i,r]) <- a0[k] + #species-level intercept
              a1.Cover*cover[k] + #average cover for a species across all plots, proxy for "abundance/size"
              a2.LifeGroup[lifegroup[k]] #combo of lifespan nad type of plant

            #Observed data are bernoulli distributed around detection probability
            #times true occupancy for that species, year, site combo
            y[k,t,i,r] ~ dbern(p[k,t,i,r]*z[k,t,i])

            #replicated data to evaluate model fit
            y.rep[k,t,i,r] ~ dbern(p[k,t,i,r]*z[k,t,i])
            
          }
        }

    }

    #SPECIES-LEVEL PRIORS
    # #Species-level Detection intercept
    a0[k] ~ dnorm(mu.a0, tau.a0)
    
    #"baseline" detection for a given species at all covariates == 0
    p0[k] <- ilogit(a0[k])
    
    #abundance model (biological process)
    #non-identifiable:
    b0.species[k] ~ dnorm(mu.b0species, tau.b0species)
    #identifiable species-level intercept:
    #track this for convergence
    b0.star[k] <- b0.species[k] + ave.eps.site[k] + ave.eps.year[k]
    
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
  
  #species-level abundance
  mu.b0species ~ dnorm(0, 0.001)
  tau.b0species <- pow(sig.b0species, -2)
  sig.b0species ~ dunif(0,50)
  
  #site and year variances
  sig.eps.site ~ dunif(0, 10)
  tau.eps.site <- pow(sig.eps.site, -2)
  sig.eps.year ~ dunif(0, 10)
  tau.eps.year <- pow(sig.eps.year, -2)

  # #Detection intercept
  mu.a0 ~ dnorm(0, 0.001)
  tau.a0 <- pow(sig.a0, -2)
  sig.a0 ~ dunif(0, 50)

  #covariate means
  a1.Cover ~ dnorm(0, 1E-3)

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

  a2.LifeGroup[1] <- 0


}
