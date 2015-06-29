posteriorDMparameters = function(c,Z,dmz,D3) {
  
  alphalist <- dmz$alpha
 
  disclass <- table(factor(c, levels = 1:K))
  activeclass <- which(disclass!=0)
  
  
  for ( i in 1:D3){
    matrixte <- list(0)
    for ( j in 1:length(activeclass)){
      matrixte[[j]] <- count(Z[c == activeclass[j],i])
    }
    
    qert <- matrix(0, nrow = length(activeclass),ncol = nlevels(as.factor(Z[,i])))
    for (j in 1:length(activeclass)){
      qert[j,matrixte[[j]]$x] <- matrixte[[j]]$freq 
    }
    ## Updating the Cluster Specific Dirichlet Priors
    alphalist[[i]][1:length(activeclass),] <- dmz$alpha0[[i]][1:length(activeclass),] + qert
  }
  
  
  list('alphalist'= alphalist)  
}
