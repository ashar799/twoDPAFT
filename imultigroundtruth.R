############ Ground Truth on TRAINING DATA WITH INTEGRATIVE METHODS ###################################
##############################################################
###########

### iCLUSTER
### k CCA 
### kmean sparse CCA


imultigroundtruth = function(){
  
  ############## Using iCluster #######
  datas <- list(0)
  datas[[1]] <- Y1
  datas[[2]] <- Y2
  cv.fit <- tune.iCluster2(datas, k)
  fit <- iCluster2(datas, k= k, lambda= cv.fit$best.fit$lambda)
  randindexiCLUSTER <<- adjustedRandIndex(fit$clusters,c.true)
  
  ########### Using CCA ################
  fit.cc <- cc(Y1, Y2)
  y1 <- fit.cc$scores$xscores
  y2 <- fit.cc$scores$yscores
  f <- length(which(fit.cc$cor > 0.5))
  Y.CCA <- cbind(y1[,1:f],y2[,1:f])
  km <- kmeans(Y.CCA, centers =k, nstart =10)
  randindexCCA <<- adjustedRandIndex(c.true,as.factor(km$cluster))
  
  
  
  }
