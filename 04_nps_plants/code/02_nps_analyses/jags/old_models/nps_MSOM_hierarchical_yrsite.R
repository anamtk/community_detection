model{
  
  #Hierarchical MSOM for NPS plants
  #Ana Miller-ter Kuile, Kiona Ogle, Shelby Lamm
  #April 10, 2024
  
  #This model is a multi-species occupancy model for understory plants
  # in the NPS Southern Colorado Plateau Inventory and Monitoring Network
  # Occupancy is based on species, year, and quadrat,with random effects
  # of site within species and year within species that are 
  # identifiable due to post-sweeping
  
  #This model has these aspects:
  ## 1. No covariates for the biological process (occupancy) model
  ## 2. Random effects on the biological process to account for species-specific
  ## site and year effects
  ## 3. Covariates on detection based on cover group
  ## 4. Two options for hierarchy based on species functional groups 
  ## and species taxonmoy
  
  for(k in 1:n.species){ #loop through all species in the community
    for(t in 1:n.years){ #loop through each sampling year
        for(i in 1:quads[t]){ #loop through the quadrats surveyed in that year
          
          #Biological process model
          #latent occupancy
          z[k,t,i] ~ dbern(psi[k,t, i])
          #occupancy probability parameter, psi, is based on species, year, and quadrat
          
          #Occupancy probability is based on baseline
          #species-level occupancy, as well as species within
          #site and species within year random effects
          logit(psi[k,t,i]) <- b0.species[k] +
            eps.site[Site.ID[t,i], k] +
            eps.year[Year.ID[t], k]
          
          #Observation model
          for(r in 1:n.rep[t,i]){ #for the number of survey reps in year t for quadrat i
            
            #Species-year-site-rep detection
            logit(p[k,t,i,r]) <- a0[k] + #species-level intercept
            a1.Cover*cover[k] #+ #average cover for a species across all plots, proxy for "abundance/size"
            #AMtK: Would probably remove this covariate for the hierarchical model
            #and instead use these catagories in the hierarchies for priors below
            #even if doing a taxonomic model - these may not be necessary, as they
            #co-vary with taxonomy
            #a2.LifeGroup[lifegroup[k]] #lifegroup
            
            #observed data are bernoulli distributed based on 
            #detection probability and latent occupancy
            y[k,t,i,r] ~ dbern(p[k,t,i,r]*z[k,t,i])
            
            #replicated data to evaluate model fit
            y.rep[k,t,i,r] ~ dbern(p[k,t,i,r]*z[k,t,i])
            
          }
        }

    }

  }
  

# Option A: Taxonomy ------------------------------------------------------
  
  #***#
  #AMtK: This is where I attempted to code in a hierarchy based 
  # on species within genus within family
  #I based this off of code on pages 10-11 of 
  #the supplemental model code for this paper:
  #https://projecteuclid.org/journals/bayesian-analysis/volume-8/issue-1/Feedback-and-Modularization-in-a-Bayesian-Metaanalysis-of-Tree-Traits/10.1214/13-BA806.full
  
  #SPECIES-LEVEL PRIORS:
  for(k in 1:n.species){
    #BIOLOGICAL Process Model Priors
    #non-identifiable:
    #original based on community-level hyperpriors
    #b0.species[k] ~ dnorm(mu.b0species, tau.b0species)
    
    #***#
    #AMtK: Is this correct?
    b0.species[k] ~ dnorm(mu.b0genus[genus[k]], tau.b0species)
    
    #***#
    #AMtK: would additional hierarchies change the
    #following post-sweeping code?
    
    #identifiable species-level intercept:
    #track this for convergence
    b0.star[k] <- b0.species[k] + ave.eps.site[k] + ave.eps.year[k]
    
    #OBSERVATION Process Model Priors
    #original based on community-level hyperpriors
    #a0[k] ~ dnorm(mu.a0, tau.a0)
    
    #***#
    #AMtK: Is this correct?
    a0[k] ~ dnorm(mu.a0genus[genus[k]], tau.a0)
    
    #"baseline" detection at all covariates == 0
    p0[k] <- ilogit(a0[k])
    
  }
  
  for(g in 1:n.genus){ #for all genuses
    
    #genus level hyperprios centered on the family
    mu.b0genus[g] ~ dnorm(mu.b0family[family[g]], tau.b0genus)
    
    mu.a0genus[g] ~ dnorm(mu.a0family[family[g]], tau.a0genus)
    
  }
  
  for(f in 1:n.family){
    #family level-hyperprior
    mu.b0family[f] ~ dnorm(0, 0.0001)
  }

# Option B: Lifegroup -----------------------------------------------------

  #***#
  #AMtK: This is where I attempted to code in a hierarchy based on
  #functional groups. I based this off of code on pages 10-11 of 
  #the supplemental model code for this paper:
  #https://projecteuclid.org/journals/bayesian-analysis/volume-8/issue-1/Feedback-and-Modularization-in-a-Bayesian-Metaanalysis-of-Tree-Traits/10.1214/13-BA806.full
  
  #SPECIES-LEVEL PRIORS:
  for(k in 1:n.species){
    #BIOLOGICAL Process Model Priors
    #non-identifiable:
    #original based on community-level hyperpriors
    #b0.species[k] ~ dnorm(mu.b0species, tau.b0species)
    
    #***#
    #AMtK: Is this correct?
    b0.species[k] ~ dnorm(mu.b0group[lifegroup[k]], tau.b0species)
    
    #***#
    #AMtK: would additional hierarchies change the
    #following post-sweeping code?
    
    #identifiable species-level intercept:
    #track this for convergence
    b0.star[k] <- b0.species[k] + ave.eps.site[k] + ave.eps.year[k]
  
    #OBSERVATION Process Model Priors
    #original based on community-level hyperpriors
    #a0[k] ~ dnorm(mu.a0, tau.a0)
    
    #***#
    #AMtK: Is this correct?
    a0[k] ~ dnorm(mu.a0group[lifegroup[k]], tau.a0)
  
    #"baseline" detection at all covariates == 0
    p0[k] <- ilogit(a0[k])
  
  }
  
  #***#
  #AMtK: then, would lifegroup level priors look like this?
  #LIFEGROUP-LEVEL PRIORS:
  for(g in 1:n.groups){ #number of life groups
    mu.b0group[g] ~ dnorm(0, 0.0001)
    mu.a0group[g] ~ dnorm(0, 0.0001)
  }
  

# Both --------------------------------------------------------------------
  
  #***#
  #AMtK: Would additional hiearchies need to be added to these
  #site and year random effects?
  
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
  

  #COMMUNITY-LEVEL HYPERPRIORS
  
  #***#
  #AMtK: Would species-level variances change at all with
  #additional hierarchies?
  
  #species-level abundance
  
  #original based on community-level hyperpriors:
  #mu.b0species ~ dnorm(0, 0.001)
  
  tau.b0species <- pow(sig.b0species, -2)
  sig.b0species ~ dunif(0,50)
  
  #site and year variances
  sig.eps.site ~ dunif(0, 10)
  tau.eps.site <- pow(sig.eps.site, -2)
  sig.eps.year ~ dunif(0, 10)
  tau.eps.year <- pow(sig.eps.year, -2)

  # #Detection intercept
  #original based on community-level hyperpriors:
  #mu.a0 ~ dnorm(0, 0.001)
  
  tau.a0 <- pow(sig.a0, -2)
  sig.a0 ~ dunif(0, 50)

  #COVARIATE PRIORS
  
  #prior for detection covariate of % cover group
  a1.Cover ~ dnorm(0, 1E-3)
  

# Additional variances for taxonomic model --------------------------------
  
  #***#
  #AMtK: Are these coded correctly for a hierarchical model?
  #I see in the paper I was referencing for above that 
  #species and genus level got folded Cauchy priors, would we 
  #want to do that here as well?
  #see page 11:
  #https://projecteuclid.org/journals/bayesian-analysis/volume-8/issue-1/Feedback-and-Modularization-in-a-Bayesian-Metaanalysis-of-Tree-Traits/10.1214/13-BA806.full

  sig.b0genus ~ dunif(0, 10)
  tau.b0genus <- pow(sig.b0genus, -2)
  sig.a0genus ~ dunif(0, 10)
  tau.a0genus <- pow(sig.b0genus, -2)


}