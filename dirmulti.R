dirmulti = function(cl, dmz, D3,N ) {
  
  logprob <- c(rep(0,D3))
  
  for ( i in 1:D3){
    
    tempA <- sum(dmz$alpha0[[i]][cl, ]) 
    tempN <- N
    tempk <-  ncol(dmz$alpha[[i]])
    logprk <- lgamma(tempA) - lgamma(tempN)
    for ( j in 1:tempk){
      logprk <- logprk + lgamma(dmz$alpha[[i]][cl,j ]) - lgamma(dmz$alpha0[[i]][cl,j ])
    }
    logprob[i] <- logprob[i] + logprk
  }
  
  
  
  return(sum(logprob))
  
  
}