Cmultiinit = function(Y,c,alphalist,N,D, sig2.dat) {
  
  alphalist <- alphalist
  betahat = matrix(data = NA, nrow = K, ncol = D)
  sigma2 <- rep(NA, K)
  beta0 <- rep(NA, K)
  
  
  source('priordraw.R')
  disclass <- table(factor(c, levels = 1:K))
  activeclass <- which(disclass!=0)
  
  ## To give some initial values to the parameters
  for ( j in 1:length(activeclass)){
    
    ## Approximating the Jeffery's prior by a Gamma distribution 
    sigma2[j] <- rinvgamma(1, shape = 1, scale = 1)
    beta0[j] <- rnorm(1, 0, sd = sig2.dat)
     for ( i in 1 :D) {
      betahat[j, i] <- rinvgamma(1, shape = 1, scale = 1) 
    }
    
  }
  
  
for ( i in 1:D){
  matrixte <- list(0)
  for ( j in 1:length(activeclass)){
  matrixte[[j]] <- count(Z[c == activeclass[j],i])
  }
      
  qert <- matrix(0, nrow = length(activeclass),ncol = nlevels(as.factor(Z[,i])))
  for (j in 1:length(activeclass)){
  qert[j,matrixte[[j]]$x] <- matrixte[[j]]$freq 
  }
  ## Updating the Cluster Specific Dirichlet Priors
  alphalist[[i]][1:length(activeclass),] <- alphalist[[i]][1:length(activeclass),] + qert
}
  
  
  list('alphalist'= alphalist, 'beta0'=beta0, 'betahat'= betahat, 'sigma2' =sigma2)  
}
