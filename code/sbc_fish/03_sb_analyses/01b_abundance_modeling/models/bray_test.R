model{

# #DERIVED PARAMETERS##

#BRAY-CURTIS DISSIMILARITY
for(i in 1:n.transects){
  for(t in (n.start[i]+1):n.end[i]){
    for(k in 1:n.species){
      # num individuals in both time periods per species
      a[k,i,t] <- min(N[k,i,t-1], N[k,i,t])
      # num individuals only in first time point
      b[k,i,t] <- N[k,i,t-1] - (min(N[k,i,t-1], N[k,i,t]))
      # num individuals only in second time point
      c[k,i,t] <- N[k,i,t] - (min(N[k,i,t-1], N[k,i,t]))
    }
    #for all years 2 onward:
    #total number of shared individuals across time periods
    A[i,t] <- sum(a[,i,t])
    #total number of individuals in only first time period
    B[i,t] <- sum(b[,i,t])
    #total number of individuals in only second time period
    C[i,t] <- sum(c[,i,t])
    
    #total bray-curtis (B+C)/(2A+B+C)
    bray[i,t] <- (B[i,t] + C[i,t])/(2*A[i,t]+B[i,t]+C[i,t])
    
    #how much is dissimilarity shaped by
    # individuals of one species being replaced by individuals
    #of another species?
    bray_balanced[i,t] <- min(B[i,t],C[i,t])/(A[i,t] + min(B[i,t],C[i,t]))
    
    #how much is dissimilarity shaped by
    # individuals that are lost without substitution?
    bray_gradient[i,t] <- bray[i,t] - bray_balanced[i,t]
  }
}
  
}



