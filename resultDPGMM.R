### This file takes the output of gene2pathways and runs my model on it

Y1 <- scale(data.gene.combined, center = TRUE, scale = TRUE)
Y2 <- scale(data.meth.combined, center = TRUE, scale = TRUE)
N <- nrow(Y1)
time <-  log(surv.time[index.patients])
censoring <- nestat[index.patients]
Time <- cbind(time, censoring)



bad.time.index = which(time == -Inf)
Y1 <- Y1[-bad.time.index,]
Y2 <- Y2[-bad.time.index,]
time <- time[-bad.time.index]
censoring <- censoring[-bad.time.index]
N <- N -1

K = as.integer(N/2)

surv.obj <- Surv(time, censoring)

############################# PARAMETERS for GIBB's SAMPLING ######################################
iter = 500
iter.burnin = 200
iter.thin  =5

## HYPER PRIORS
## Hyper parameters of the DP
shape.alpha <- 2
rate.alpha <- 1
## Hyperparameters for the GMM
beta  = max(D1,D2) + 1
ro = 0.5
## Initialize the c using chinese restaurant process
source('rchinese.R')
alpha  = rgamma(1, shape = shape.alpha, rate = rate.alpha )
c <-  rchinese(N,alpha)
f <- table(factor(c, levels = 1:max(c)))
#Sparsity controlling hyperparameter of the BAYESIAN LASSO MODEL
r =1
si = 1.78


### LETS MAKE A LIST "gmmx" to store parameters/hyperprameters for X and "regy" to store paameters for Regression Y
## For the First Data Set
gmmx1 <- list(0)
gmmx1$epsilon <-  as.vector(apply(Y1,2,mean))
gmmx1$W <- cov(Y1)
gmmx1$mu <- matrix(data = NA, nrow = K, ncol = D1)
gmmx1$S <-  array(data = NA, dim =c(K,D1,D1))

regy1 <- list(0)
regy1$lambda2 <- numeric(K)
regy1$tau2 = matrix(data = NA, nrow = K, ncol = D1)
regy1$betahat = matrix(data = NA, nrow = K, ncol = D1)
regy1$sigma2 <- rep(NA, K)
regy1$beta0 <- rep(NA, K)

## For the second data set
gmmx2 <- list(0)
gmmx2$epsilon <-  as.vector(apply(Y2,2,mean))
gmmx2$W <- cov(Y2)
gmmx2$mu <- matrix(data = NA, nrow = K, ncol = D2)
gmmx2$S <-  array(data = NA, dim =c(K,D2,D2))

regy2 <- list(0)
regy2$lambda2 <- numeric(K)
regy2$tau2 = matrix(data = NA, nrow = K, ncol = D2)
regy2$betahat = matrix(data = NA, nrow = K, ncol = D2)
regy2$sigma2 <- rep(NA, K)
regy2$beta0 <- rep(NA, K)


###### To initialize the parameters for all the data sets
That <-  time
####### We can use a simple Linear Model to get some estimates of the variance##########
Yg <- cbind(Y1,Y2)
Dg <- (D1 + D2)

## Fitting a linear model to the whole model
Ysc <- scale(Yg[1:N,1:Dg], center = TRUE, scale =TRUE)
lm.data <- lm(time ~ Ysc)
sig2.dat <-  var(lm.data$residuals)


## Set Some Initial Values for the Cluster Parameters
source('multiinit.R')

## For the first data set
cont1 <- multiinit(Y1,c,beta, gmmx1$W, gmmx1$epsilon, ro, r, si,N,D1, sig2.dat)
gmmx1$mu <- cont1$mu
gmmx1$S <- cont1$S
regy1$lambda2 <- cont1$lambda2
regy1$tau2 <- cont1$tau2
regy1$betahat <- cont1$betahat
regy1$sigma2 <- cont1$sigma2
regy1$beta0 <- cont1$beta0

## For the second data set
cont2 <- multiinit(Y2,c,beta, gmmx2$W, gmmx2$epsilon, ro, r, si,N,D2, sig2.dat)
gmmx2$mu <- cont2$mu
gmmx2$S <- cont2$S
regy2$lambda2 <- cont2$lambda2
regy2$tau2 <- cont2$tau2
regy2$betahat <- cont2$betahat
regy2$sigma2 <- cont2$sigma2
regy2$beta0 <- cont2$beta0





## Initialization part for the parmaters of AFT Model with k-means and Bayesian Lasso and Normal Bayesian Regression

lik = c(0)
c.init <- c
gmmx1.int <- gmmx1
gmmx2.int <- gmmx2
regy1.int <- regy1
regy2.int <- regy2
likint <- c(0)
bic <- c(0)


# ## Initialization
source('multilikelihood.R')
source('multikmeansBlasso.R')


km <- multikmeansBlasso(c.init,Y1,Y2,D1,D2,That,K, beta, ro, r, si,sig2.dat,gmmx1.int, gmmx2.int, regy1.int, regy2.int,surv.obj )
c <- km$c
gmmx1 <- km$gmmx1
gmmx2 <- km$gmmx2 
regy1 <- km$regy1
regy2 <- km$regy2
likint0 <- multiloglikelihood( c,Y1,Y2,D1,D2,That,K, beta, ro, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2)



### Initial p-value separation
logrank0 <- survdiff(surv.obj ~ c)




## Gibb's sampling 

source('priordraw.R')
cognate <- NA
param <- NA
paramtime <- NA
loglike<- rep(0, iter)  
timeparam <- NA
time.predicted <- c(0)
cindex <- c(0)



logrank <- c(0)
likli <- c(0)

o.init <- o
#################### BURNIN PHASE ###################################################
print("BURNIN...PHASE")
for (o in o.init:iter.burnin) {
  
  
  ################## PARAMETERS OF THE DP Mixture Model ######################################################
  ## Updating the parameters based on the observations 
  source('posteriorGMMparametrs.R')
  param <- posteriorGMMparametrs(c,Y1,gmmx1$mu,gmmx1$S, alpha, K, gmmx1$epsilon, gmmx1$W, beta, ro,N,D1 )
  gmmx1$mu <- param$mean
  gmmx1$S <- param$precision
  param2 <- posteriorGMMparametrs(c,Y2,gmmx2$mu,gmmx2$S, alpha,K, gmmx2$epsilon, gmmx2$W, beta, ro,N,D2 )
  gmmx2$mu <- param2$mean
  gmmx2$S <- param2$precision
  
  
  source('posteriortimeparameterspenalized.R')
  paramtime1 <- posteriortimeparameterspenalized(c,Y1, That, regy1$lambda2, regy1$tau2, regy1$sigma2, regy1$beta0, regy1$betahat, K, gmmx1$epsilon, gmmx1$W, beta, ro, r, si, sig2.data,N, D1)
  regy1$beta0 <- paramtime1$beta0
  regy1$betahat <- paramtime1$betahat
  regy1$sigma2 <- paramtime1$sigma2
  regy1$lambda2 <- paramtime1$lambda2
  regy1$tau2 <- paramtime1$tau2
  
  paramtime2 <- posteriortimeparameterspenalized(c,Y2, That, regy2$lambda2, regy2$tau2, regy2$sigma2, regy2$beta0, regy2$betahat, K, gmmx2$epsilon, gmmx2$W, beta, ro, r, si, sig2.data,N, D2)
  regy2$beta0 <- paramtime2$beta0
  regy2$betahat <- paramtime2$betahat
  regy2$sigma2 <- paramtime2$sigma2
  regy2$lambda2 <- paramtime2$lambda2
  regy2$tau2 <- paramtime2$tau2
  
  
  
  ########################## THE HYPERPARAMETERS OF THE GMM #################################  
  source('posteriorhyper.R')  
  
  # Updating the hyper paramters for the first data set
  hypercognate <- posteriorhyper (c, Y1, gmmx1$mu, gmmx1$S, gmmx1$epsilon, gmmx1$W, beta, ro,D1 )
  gmmx1$epsilon <- hypercognate$epsilon
  tmpW <- hypercognate$W
  gmmx1$W <- matrix(as.matrix(tmpW),nrow = D1, ncol =D1)
  
  ##Updating the hyper parameter for the second data set
  hypercognate2 <- posteriorhyper (c, Y2, gmmx2$mu, gmmx2$S, gmmx2$epsilon, gmmx2$W, beta, ro,D2 )
  gmmx2$epsilon <- hypercognate2$epsilon
  tmpW2 <- hypercognate2$W
  gmmx2$W <- matrix(as.matrix(tmpW2),nrow = D2, ncol =D2)
  
  
  
  
  ################# INDICATOR VARIABLE ##################################################################
  ## Updating the indicator variables and the parameters
  source('multiposteriorchineseAFT.R')  
  cognate <- multiposteriorchineseAFT(c,Y1,Y2,D1,D2,That, F,K, beta, ro, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2)
  c <- cognate$c
  gmmx1 <- cognate$gmmx1
  gmmx2 <- cognate$gmmx2
  regy1 <- cognate$regy1
  regy2 <- cognate$regy2
  
  
  
  ########################### The Concentration Parameter #################################################################
  
  
  source('posterioralpha.R') 
  # Updating the concentration parameter
  alpha <- posterioralpha(c, N, alpha, shape.alpha, rate.alpha)
  
  
  ### Calculating the number of active cluster and log rank #####################
  logrankt <- survdiff(surv.obj ~ c)
  G =  length(which(table(factor(c, levels = 1:K))!=0))
  logrank[o] <- 1 - pchisq(logrankt$chisq, (G-1))
  
  source('multilikelihood.R')
  likli[o] <- multiloglikelihood( c,Y1,Y2,D1,D2,That,K, beta, ro, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2)
  
  ##################### Print SOME Statistics #####################################################
  
  print(logrank[o])
  print(likli[o])
  print(o/iter.burnin)
  
} 







#####################   GIBBS SAMPLING WITH THINNING ######################################################

est.regy1 <- list(0)
est.regy2 <- list(0)
est.gmmx1 <- list(0)
est.gmmx2 <- list(0)
c.list <- list(0)





print("GIBB'S SAMPLING")
count = 1
for (o in 1:iter) {
  ################## PARAMETERS OF THE DP Mixture Model ######################################################
  ## Updating the parameters based on the observations 
  source('posteriorGMMparametrs.R')
  param <- posteriorGMMparametrs(c,Y1,gmmx1$mu,gmmx1$S, alpha, K, gmmx1$epsilon, gmmx1$W, beta, ro,N,D1 )
  gmmx1$mu <- param$mean
  gmmx1$S <- param$precision
  param2 <- posteriorGMMparametrs(c,Y2,gmmx2$mu,gmmx2$S, alpha,K, gmmx2$epsilon, gmmx2$W, beta, ro,N,D2 )
  gmmx2$mu <- param2$mean
  gmmx2$S <- param2$precision
  
  
  source('posteriortimeparameterspenalized.R')
  paramtime1 <- posteriortimeparameterspenalized(c,Y1, That, regy1$lambda2, regy1$tau2, regy1$sigma2, regy1$beta0, regy1$betahat, K, gmmx1$epsilon, gmmx1$W, beta, ro, r, si, sig2.data,N, D1)
  regy1$beta0 <- paramtime1$beta0
  regy1$betahat <- paramtime1$betahat
  regy1$sigma2 <- paramtime1$sigma2
  regy1$lambda2 <- paramtime1$lambda2
  regy1$tau2 <- paramtime1$tau2
  
  paramtime2 <- posteriortimeparameterspenalized(c,Y2, That, regy2$lambda2, regy2$tau2, regy2$sigma2, regy2$beta0, regy2$betahat, K, gmmx2$epsilon, gmmx2$W, beta, ro, r, si, sig2.data,N, D2)
  regy2$beta0 <- paramtime2$beta0
  regy2$betahat <- paramtime2$betahat
  regy2$sigma2 <- paramtime2$sigma2
  regy2$lambda2 <- paramtime2$lambda2
  regy2$tau2 <- paramtime2$tau2
  
  
  
  ########################## THE HYPERPARAMETERS OF THE GMM #################################  
  source('posteriorhyper.R')  
  
  # Updating the hyper paramters for the first data set
  hypercognate <- posteriorhyper (c, Y1, gmmx1$mu, gmmx1$S, gmmx1$epsilon, gmmx1$W, beta, ro,D1 )
  gmmx1$epsilon <- hypercognate$epsilon
  tmpW <- hypercognate$W
  gmmx1$W <- matrix(as.matrix(tmpW),nrow = D1, ncol =D1)
  
  ##Updating the hyper parameter for the second data set
  hypercognate2 <- posteriorhyper (c, Y2, gmmx2$mu, gmmx2$S, gmmx2$epsilon, gmmx2$W, beta, ro,D2 )
  gmmx2$epsilon <- hypercognate2$epsilon
  tmpW2 <- hypercognate2$W
  gmmx2$W <- matrix(as.matrix(tmpW2),nrow = D2, ncol =D2)
  
  
  
  
  ################# INDICATOR VARIABLE ##################################################################
  ## Updating the indicator variables and the parameters
  source('multiposteriorchineseAFT.R')  
  cognate <- multiposteriorchineseAFT(c,Y1,Y2,D1,D2,That, F,K, beta, ro, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2)
  c <- cognate$c
  gmmx1 <- cognate$gmmx1
  gmmx2 <- cognate$gmmx2
  regy1 <- cognate$regy1
  regy2 <- cognate$regy2
  
  
  
  ########################### The Concentration Parameter #################################################################
  
  
  source('posterioralpha.R') 
  # Updating the concentration parameter
  alpha <- posterioralpha(c, N, alpha, shape.alpha, rate.alpha)
  
  
  ##################### Print SOME Statistics #####################################################
  randy[o] <- adjustedRandIndex(c.true,as.factor(c))
  print(randy[o])
  likli[o] <- multiloglikelihood( c,Y1,Y2,D1,D2,That, F,K, beta, ro, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2)
  print(likli[o])
  print(o/iter.burnin)
  
  
  
  if(o%% iter.thin == 0 ){
    est.regy1[[count]] <- regy1
    est.regy2[[count]] <- regy2
    est.gmmx1[[count]] <- gmmx1
    est.gmmx2[[count]] <- gmmx2
    count <- count +1
  }
  
  
  
  
  print(o/iter) 
  #   print(loglike[o])
  #   print(cindex)
} 
