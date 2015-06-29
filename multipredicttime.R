multipredicttime = function(c, Yg, That, Time, beta0, betahatg, sigma2 ) {

Ytemp <- matrix(NA, nrow = N, ncol = Dg)
numclust <- table(factor(c, levels = 1:K))
activeclass<- which(numclust!=0)

prediction <- c(0)


for ( i in 1:length(activeclass)) {
  
  clust <- which(c == activeclass[i])
  
  Ytemp[clust,1:Dg] <- scale(Yg[clust,1:Dg], center = TRUE, scale = TRUE)
}

    
 for ( h in 1:N){  
   
   if (Time[h,2]==0) {
     prediction[h]<- rtruncnorm(1, a = Time[h,1], b = Inf, mean = beta0[c[h]] + betahatg[c[h],1:Dg ] %*% Ytemp[h,1:Dg] , sd = sqrt(sigma2[c[h]]) )
    
  } else {
    prediction[h] <- rnorm(1, mean = beta0[c[h]] + betahatg[c[h],1:Dg ] %*% Ytemp[h,1:Dg], sd = sqrt(sigma2[c[h]]))
}


}

list('predicttime' = prediction) 


}
