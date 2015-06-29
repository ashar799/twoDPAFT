## This function generates categorical data

generatedata2 = function(N,D,F,p.dist) {
  
  
  N = N
  
  ## Number of Clusters
  F = F
  
  ## Distribution of the points within three clusters
  
  p.dist = p.dist
  
  ## Total Number of features D
  
  D = D
  
  ## Assuming the type of discretization is varying between 2and 5 for each variable type
  prob_list <- list(0)
  z_list <- list(0)
  
 
  
  ## The number of components for each variable  
  q = sample(c(1:3),size = D, replace =TRUE, prob = c(0.25,0.5,0.25))
    
  
  
  
  
 for (l in 1:F){
   
   ## A list which contains cluster specific cumaltive distributions
   marginal_list <- list(0)
   
   for ( i in 1:D){
     marginal_list[[i]] <- pnorm(sort(runif(q[i],min= -1,max =1)))
   }
   
   ## Correlation matrix 
   R <- diag(D) # Correlation matrix
   prob_list[[l]] <- marginal_list
   
   nt <- as.integer(N*p.dist[l])
   z_list[[l]] <- ordsample(nt, marginal_list, R)
   
   
 }
  
 Z<- c(0)
 for (i in 1:F){
   Z <- rbind(Z, z_list[[i]])
 } 
 Z <- Z[-1,]
 
 
 Z.rel.sc.list <- list(0)
 for ( i in 1:F){
   Z.rel.sc.list[[i]] <- scale(z_list[[i]], center = TRUE, scale = TRUE)
 }
 
 
 
 
 ## The Co-efficients have to be obtained from uniform distribution between [-3,3]
 beta.list <- list(0)
 half <- as.integer(D/2)
 ohalf <- D - half
 for ( i in 1:F){
   beta.list[[i]] <- c(runif(half, min = -3, max = -0.1), runif(ohalf, min = 0.1, max = 3))
 }
 
 ## The pure time is generated
 time.pur.list <- list(0)
 for ( i in 1:F){
   time.pur.list[[i]] <- t(beta.list[[i]]) %*% t((Z.rel.sc.list[[i]])) 
 }
 
 
 
 
 
 list('Z' = Z, 'beta' = beta.list, 'timepur'= time.pur.list)

}
