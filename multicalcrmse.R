multicalcrmse = function(c,N,Y1,D1,Y2,D2,Z,D3, reg1, reg2, reg3) {
  
  
  ## Let's calculate MSE between predicted times with the correct co-efficients and Real Time i.e. without censoring
  Ytemp1 <- matrix(NA, nrow = N, ncol = D1)
  Ytemp2 <- matrix(NA, nrow = N, ncol = D2)
  Ztemp <-  matrix(NA, nrow = N, ncol = D3)
  numclust <- table(factor(c.true, levels = 1:F))
  activeclass<- which(numclust!=0)
  for ( i in 1:length(activeclass)) {
    clust <- which(c.true == activeclass[i])
    Ytemp1[clust,1:D] <- scale(Y1[clust,1:D], center = TRUE, scale = TRUE)
    Ytemp2[clust,1:D] <- scale(Y2[clust,1:D], center = TRUE, scale = TRUE)
    Ztemp[clust,1:D] <-  scale(Y2[clust,1:D], center = TRUE, scale = TRUE)
  }
  
  
  
  
  time.predicted1 <- c(0)
  time.predicted2 <- c(0)
  time.predicred3 <- c(0)
  for ( i in 1:N){  
    time.predicted1[i] <-  rnorm(1, mean = time.cluster$Mu[c.true[i]] + t(beta.list[[c.true[i]]]) %*% Ytemp1[i,1:rel.D] +  t(beta.list2[[c.true[i]]]) %*% Ytemp2[i,1:rel.D], sd = sqrt(time.cluster$S[c.true[i]]))
  }
  
  source('calcrmse.R')
  er <- calcrmse(time.real,time.predicted)$rmse
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}