posteriortimeparametersUNpenalized = function(c,Z, regz,N,D3 ) {
  
  numclust <- table(factor(c, levels = 1:K))
  activeclass<- which(numclust!=0)
  regz <- regz
  
  for (j in 1:length(activeclass)) {
    
    
    ### Part where I use the MONOMVN PACKAGE
    ## If the Cluster has just one member I draw the parameters from the prior 
    if (numclust[activeclass[j]] > 1){
        
        pos <- MCMCregress(That[c== activeclass[j]] ~ Z[c== activeclass[j],]) 
        value <- unlist(summary(pos, quantiles = c(0.5))[2])
        regz$beta0[activeclass[j]] <- value[1]
        indt <- 2:(D3+1)
        regz$betahat[activeclass[j],] <- as.vector(value[indt])
        regz$sigma2[activeclass[j]] <- as.numeric(value[D3+2])
      
      } else {
      
        regz$sigma2[activeclass[j]] <- rinvgamma(1, shape = 1, scale = 1)
        regz$beta0[activeclass[j]] <- rnorm(1, 0, sd = sig2.dat)
        for ( i in 1 :D3) {
          regz$betahat[activeclass[j], i] <- rinvgamma(1, shape = 1, scale = 1) 
        }
        
      }
    
  }
  
  
list('regz' = regz )
}
