## Burnin Iterations for the multi view DPMM 

burninmultiDPMM = function(){
  
source('priordraw.R')
param <- NA
paramtime1 <- NA
paramtime2 <- NA
cognate <- NA
hypercognate1 <- NA
hypercognate2 <- NA
loglike<- rep(0, iter)  




randy <- c(0)
likli <- c(0)

#################### BURNIN PHASE ###################################################
print("BURNIN...PHASE")
for (o in 1:iter.burnin) {
  
  
  ################## PARAMETERS OF THE DP Mixture Model ######################################################
  ## Updating the parameters based on the observations 
  source('posteriorGMMparametrs.R')
  param <- posteriorGMMparametrs(c,Y1,gmmx1$mu,gmmx1$S, alpha, K, gmmx1$epsilon, gmmx1$W, gmmx1$beta, gmmx1$ro,N,D1 )
  gmmx1$mu <- param$mean
  gmmx1$S <- param$precision
  param2 <- posteriorGMMparametrs(c,Y2,gmmx2$mu,gmmx2$S, alpha,K, gmmx2$epsilon, gmmx2$W, gmmx2$beta, gmmx2$ro,N,D2 )
  gmmx2$mu <- param2$mean
  gmmx2$S <- param2$precision
  
  
  source('posteriortimeparameterspenalized.R')
  paramtime2 <- posteriortimeparameterspenalized(c,Y2, That, regy2$lambda2, regy2$tau2, regy2$sigma2, regy2$beta0, regy2$betahat, K, gmmx2$epsilon, gmmx2$W, gmmx2$beta, gmmx2$ro, r, si, sig2.data,N, D2)
  regy2$beta0 <- paramtime2$beta0
  regy2$betahat <- paramtime2$betahat
  regy2$sigma2 <- paramtime2$sigma2
  regy2$lambda2 <- paramtime2$lambda2
  regy2$tau2 <- paramtime2$tau2
  
  
  paramtime1 <- posteriortimeparameterspenalized(c,Y1, That, regy1$lambda2, regy1$tau2, regy1$sigma2, regy1$beta0, regy1$betahat, K, gmmx1$epsilon, gmmx1$W, gmmx1$beta, gmmx1$ro, r, si, sig2.data,N, D1)
  regy1$beta0 <- paramtime1$beta0
  regy1$betahat <- paramtime1$betahat
  regy1$sigma2 <- paramtime1$sigma2
  regy1$lambda2 <- paramtime1$lambda2
  regy1$tau2 <- paramtime1$tau2
  
 
  
  
  ########################## THE HYPERPARAMETERS OF THE GMM #################################  
  source('posteriorhyperplus.R')  
  
  # Updating the hyper paramters for the first data set
  hypercognate <- posteriorhyperPLUS(c, Y1, gmmx1$mu, gmmx1$S, gmmx1$epsilon, gmmx1$W, gmmx1$beta, gmmx1$ro,D1 )
  gmmx1$epsilon <- hypercognate$epsilon
  tmpW <- hypercognate$W
  gmmx1$W <- matrix(as.matrix(tmpW),nrow = D1, ncol =D1)
  gmmx1$ro <- hypercognate$ro
  
  
  
  ##Updating the hyper parameter for the second data set
  hypercognate2 <- posteriorhyperPLUS(c, Y2, gmmx2$mu, gmmx2$S, gmmx2$epsilon, gmmx2$W, gmmx2$beta, gmmx2$ro,D2 )
  gmmx2$epsilon <- hypercognate2$epsilon
  tmpW2 <- hypercognate2$W
  gmmx2$W <- matrix(as.matrix(tmpW2),nrow = D2, ncol =D2)
  gmmx2$ro <- hypercognate2$ro
  
  
  ### Updating Beta parameter for the first view #################
#   source('posteriorbeta.R')
#   if( o%%10 == 0){
#     res <- try(posteriorbeta(c, gmmx1$beta, D1, gmmx1$S, gmmx1$W))
#     if (class(res) == "try-error"){
#       gmmx1$beta = gmmx1$beta
#     } else{
#       gmmx1$beta <- posteriorbeta(gmmx1$beta, D1, gmmx1$S, gmmx1$W)
#       
#     }
#   } 
#   ### Updating Beta parameter for the second view #################
#   source('posteriorbeta.R')
#   if( o%%10 == 0){
#     res <- try(posteriorbeta(c, gmmx2$beta, D2, gmmx2$S, gmmx2$W))
#     if (class(res) == "try-error"){
#       gmmx2$beta = gmmx2$beta
#     } else{
#       gmmx2$beta <- posteriorbeta(gmmx2$beta, D2, gmmx2$S, gmmx2$W)
#       
#     }
#   } 
#   
  
  
  ################# INDICATOR VARIABLE ##################################################################
  ## Updating the indicator variables and the parameters
  source('multiposteriorchineseAFT.R')  
  cognate <- multiposteriorchineseAFT(c,Y1,Y2,D1,D2,That, K, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2)
  c <- cognate$c
  gmmx1 <- cognate$gmmx1
  gmmx2 <- cognate$gmmx2
  regy1 <- cognate$regy1
  regy2 <- cognate$regy2
  
  
  
  ########################### The Concentration Parameter #################################################################
  
  
  source('posterioralpha.R') 
  # Updating the concentration parameter
  alpha <- posterioralpha(c, N, alpha, shape.alpha, rate.alpha)
  
  
  
  ######################## The Censored Times ###########################################################
#   source('multiupdatetime.R')
#   # Updating the Time Variable
#   ti <- NA
#   ti <- multiupdatetime(c, Y1, Y2, Time,That, regy1, regy2)
#   That <- ti$time
#   
#   
#   
  
  ##################### Print SOME Statistics #####################################################
  #randy[o] <- adjustedRandIndex(c.true,as.factor(c))
  #print(randy[o])
  likli[o] <- multiloglikelihood(c,Y1,Y2,D1,D2,That,K, beta, ro, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2)
  print(likli[o])
  print(o/iter.burnin)
  
} 

assign("alpha", alpha, envir = .GlobalEnv)
assign("gmmx1", gmmx1, envir = .GlobalEnv)
assign("gmmx2", gmmx2, envir = .GlobalEnv)
assign("regy1", regy1, envir = .GlobalEnv)
assign("regy2", regy2, envir = .GlobalEnv)
assign("c", c, envir = .GlobalEnv)
assign("randy.burnin", randy, envir = .GlobalEnv)
assign("likli.burnin", likli, envir = .GlobalEnv)

plot(likli, main = 'Burnin Iterations')

}

