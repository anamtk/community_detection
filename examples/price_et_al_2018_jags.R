###############################################################################
#                                                                             #
# This file accompanies:                                                      #
#                                                                             #
# Price, Freytag, Bonner, Drayer, Muncy, Hutton and Barton (2018) Mountaintop #
#    removal mining alters stream salamander population dynamics. Diversity   #
#    and Distributions.                                                       #
#                                                                             #
###############################################################################

## Notes:
##    S indexes species
##    s indexes site
##    y indexes year
##    v indexes visit

model{
    ##### Likelihood #####

    ##### Occupancy #####

    for(S in 1:nspecies){
        ## 1) Mean parameters

        for(mtr in 0:1){
            ## a) Initial occupancy
            eta.psi0[S,mtr+1] <- beta.occ[1,mtr+1,S]
            logit(psi0.pop[S,mtr+1]) <- eta.psi0[S,mtr+1]

            ## b) Colonization
            eta.gamma[S,mtr+1] <- beta.occ[2,mtr+1,S]
            logit(gamma.pop[S,mtr+1]) <- eta.gamma[S,mtr+1]

            ## c) Survival
            eta.phi[S,mtr+1] <- beta.occ[3,mtr+1,S]
            logit(phi.pop[S,mtr+1]) <- eta.phi[S,mtr+1]
        }
        
        ## 2) Site specific parameters
        for(s in 1:nsite){
            ## Initial occupancy
            logit(psi[S,s,1]) <- eta.psi0[S,MTR[s]+1]

            Occupancy[S,s,1] ~ dbern(psi[S,s,1])
            
            ## Occupancy in subsequent years
            for(y in 2:nyear){
                ## Colonization
                logit(gamma[S,s,y-1]) <-  eta.gamma[S,MTR[s]+1]

                ## Survival
                logit(phi[S,s,y-1]) <-  eta.phi[S,MTR[s]+1]

                ## Occupancy
                psi[S,s,y] <- (1-Occupancy[S,s,y-1]) * gamma[S,s,y-1] +
                    Occupancy[S,s,y-1] * phi[S,s,y-1]
                
                Occupancy[S,s,y] ~ dbern(psi[S,s,y])
            }
        }
    }

    ## Abundance given occupancy
    for(S in 1:nspecies){
        log(lambda.pop[S,1]) <- beta.abund[1,S]
        log(lambda.pop[S,2]) <- beta.abund[2,S]
        
        for(s in 1:nsite){
            for(y in 1:nyear){
                log(lambda[S,s,y]) <- beta.abund[1,S] * (1-MTR[s]) +
                    beta.abund[2,S] * MTR[s]
                
                Abundance.tmp[S,s,y] ~ dpois(lambda[S,s,y])T(1,)

                Abundance[S,s,y] <- Abundance.tmp[S,s,y]*Occupancy[S,s,y]
            }
        }
    }

    ## Detection
    for(S in 1:nspecies){
        for(s in 1:nsite){
            for(y in 1:nyear){
                for(v in 1:nvisit[y]){
                    logit(p[S,s,y,v]) <- beta.det[1,S] +
                        beta.det[2,S] * CoverObjects[s,y] + 
                        beta.det[3,S] * Precip[s,y,v]
                }
            }
        }
    }

    ## Observations
    for(i in 1:nobs){
        for(S in 1:nspecies){
            Y[i,S] ~ dbinom(p[S,Site[i],Year[i],Visit[i]],
                            Abundance[S,Site[i],Year[i]])
        }
    }

    ##### Priors #####

    ## Parameters for half-t priors on variance
    df <- 3
    tau <- .25
    
    ## Occupancy
    for(i in 1:3){ # 1=Initial, 2=Colonization, 3=Survival
        for(k in 1:2){ # 1=Intercept, 2=MTR
            for(S in 1:nspecies){
                xi.occ[i,k,S] ~ dnorm(0,tau.beta.occ[i,k])
                beta.occ[i,k,S] <- mu.beta.occ[i,k] + alpha.beta.occ[i,k] * xi.occ[i,k,S]
            }
            
            mu.beta.occ[i,k] ~ dnorm(0,.36)
            tau.beta.occ[i,k] ~ dgamma(df/2,df/2/tau)
            sigma.beta.occ[i,k] <- abs(alpha.beta.occ[i,k])/sqrt(tau.beta.occ[i,k])
            alpha.beta.occ[i,k] ~ dnorm(0,1)
        }
    }
    
    ## Abundance
    for(i in 1:2){
        for(S in 1:nspecies){
            xi.abund[i,S] ~ dnorm(0,tau.beta.abund[i])
            beta.abund[i,S] <- mu.beta.abund[i] + alpha.beta.abund[i] * xi.abund[i,S]
        }
        
        mu.beta.abund[i] ~ dnorm(0,.001)
        tau.beta.abund[i] ~ dgamma(df/2,df/2/tau)
        sigma.beta.abund[i] <- abs(alpha.beta.abund[i])/sqrt(tau.beta.abund[i])
        alpha.beta.abund[i] ~ dnorm(0,1)
    }

    ## Detection
    for(i in 1:3){
        for(S in 1:nspecies){
            xi.det[i,S] ~ dnorm(0,tau.beta.det[i])
            beta.det[i,S] <- mu.beta.det[i] + alpha.beta.det[i] * xi.det[i,S]
        }

        mu.beta.det[i] ~ dnorm(0,.36)
        tau.beta.det[i] ~ dgamma(df/2,df/2/tau)
        sigma.beta.det[i] <- abs(alpha.beta.det[i])/sqrt(tau.beta.det[i])
        alpha.beta.det[i] ~ dnorm(0,1)
    }
    

    ##### Derived Values #####
    ## Percent occupancy
    for(S in 1:nspecies){
        for(s in 1:nsite){
            OccSum[S,s] <- sum(Occupancy[S,s,])
        }
        
        PercOcc[S,1] <- 100*inprod(OccSum[S,],(1-MTR[]))/(sum(1-MTR[])*nyear)
        PercOcc[S,2] <- 100*inprod(OccSum[S,],MTR[])/(sum(MTR[])*nyear)
    }

    ## Mean abundance given occupancy
    for(S in 1:nspecies){
        for(s in 1:nsite){
            AbundMean.tmp[S,s] <- mean(Abundance[S,s,])
        }
        
        AbundMean[S,1] <- inprod(AbundMean.tmp[S,],(1-MTR[]))/(sum(1-MTR[])*nyear)
        AbundMean[S,2] <- inprod(AbundMean.tmp[S,],MTR[])/(sum(MTR[])*nyear)
    }
}
    
