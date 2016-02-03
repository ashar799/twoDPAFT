###### THIS FILE COMBINES gene2pathway.R and resultDPGMM.R and creates ONE file ##########

###### This script calculates the Pathway level scores for both Gene Expression Data && Methylation Data   ####
rm(list = ls())
load('/home/bit/ashar/ExpressionSets/Verhark/IntegratedDataSurv.rda')


#### Libraries for gene2pathway ##############################################
library(globaltest)
library(org.Hs.eg.db)
library(KEGG.db)
library(IlluminaHumanMethylation450k.db)
library(survival)

## Defining the survival time object

status <- clinical[,5]
nestat <- c(0)
ind1 <- which(status == levels(status)[1] | status == levels(status)[3] )
ind2 <- which(status == levels(status)[2] | status == levels(status)[4] )
for ( i in 1:length(ind1)){
  nestat[ind1[i]] <- 1
}
for ( i in 1:length(ind2)){
  nestat[ind2[i]] <- 0
}
surv.time <- as.numeric(as.matrix(clinical[,4]))
surv.obj <- Surv(surv.time,nestat)


## SELECT COMMON PATIENTS
main.index   <- match(colnames(meth),colnames(gene.expr))
main.index <- main.index[!is.na(main.index)]
gene.expr.main <- t(gene.expr)[main.index,]


main.meth.index <- match(colnames(gene.expr),colnames(meth))
main.meth.index <- main.meth.index[!is.na(main.meth.index)]
meth.expr.main <- t(meth)[main.meth.index,]

###JUST CHECKING
rownames(gene.expr.main) == rownames(meth.expr.main)

## Make different survival objects
index.patients <- match(rownames(gene.expr.main),mypatients)
surv.obj.gene <- surv.obj[index.patients]
surv.obj.meth <-  surv.obj[index.patients]


## Gene Expression Data
## DO PATHWAY LEVEL ENRICHMENT
xx <- as.list(org.Hs.egALIAS2EG)
# Remove pathway identifiers that do not map to any entrez gene id
xx <- xx[!is.na(xx)]
glob1 <- gtKEGG( surv.obj.gene, gene.expr.main, annotation = 'org.Hs.eg.db', multtest = "BH",probe2entrez = xx)
## Select the significant ones p < 0.01
pathway.gene.names <- names(glob1)[p.value(glob1) < 0.005]
yy1 <- as.list(KEGGPATHID2EXTID)
pathway.gene.names <- pathway.gene.names[!is.na(pathway.gene.names)]
## Methylation Data


xx2 <- IlluminaHumanMethylation450kENTREZID
mapped_probes <- mappedkeys(xx2)
# Convert to a list
xx2 <- as.list(xx2[mapped_probes])
glob2 <- gtKEGG(surv.obj.meth, meth.expr.main, annotation = 'org.Hs.eg.db', multtest = "BH",probe2entrez = xx2)
pathway.meth.names <- names(glob2)[p.value(glob2) < 0.005]
pathway.meth.names <- pathway.meth.names[!is.na(pathway.meth.names)]
yy2 <- as.list(KEGGPATHID2EXTID)


######### TO GET GENE NAMES AND SUMMARIZE THE SCORES of GENES WITHIN PATHWAYS ##############################3

# For the reverse map:
# Convert the object to a list
dd <- as.list(org.Hs.egPATH2EG)
# Remove pathway identifiers that do not map to any entrez gene id
dd <- dd[!is.na(dd)]

path2entrez <- dd

path2entrez.subset <- path2entrez[pathway.gene.names]

library('biomaRt')
# listDatasets(ensembl)
ensembl = useMart('ensembl')
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

names.genes <- list(0)
for( i in 1:length(pathway.gene.names)){
  names.genes[[i]] = getBM(attributes = c("entrezgene","hgnc_symbol"), filters = "entrezgene",values =path2entrez.subset[i], mart=ensembl)
}

### Data frame FROM the EXPRESSION #################

data.subsets <- list(0)
for ( i in 1:length(pathway.gene.names) ){
  
  temp.names <- names.genes[[i]][,2][c(names.genes[[i]][,2])%in% colnames(gene.expr.main)]
  data.subsets[[i]] <- as.data.frame(gene.expr.main[,temp.names]) 
}

data.gene.combined <- matrix(0, nrow = nrow(gene.expr.main), ncol =length(pathway.gene.names) )
for ( i in 1:length(pathway.gene.names)){
  data.gene.combined[,i] <- as.vector(apply(data.subsets[[i]],1, mean))
}
rownames(data.gene.combined) <- rownames(gene.expr.main)
colnames(data.gene.combined) <- pathway.gene.names
########################################################################
#########################################################################
## Methylation Data
# For the reverse map:
# Convert the object to a list
dd2 <- as.list(IlluminaHumanMethylation450kPATH2PROBE)
# Remove pathway identifiers that do not map to any entrez gene id
dd2 <- dd2[!is.na(dd2)]

path2cg <- dd2

path2cg.subset <- path2cg[pathway.meth.names]



data.subsets <- list(0)

for ( i in 1:length(pathway.meth.names) ){
  inde <- colnames(meth.expr.main) %in% path2cg.subset[[i]]
  data.subsets[[i]] <- as.data.frame(meth.expr.main[,inde]) 
}


data.meth.combined <- matrix(0, nrow = nrow(meth.expr.main), ncol =length(pathway.meth.names) )
for ( i in 1:length(pathway.meth.names)){
  data.meth.combined[,i] <- as.vector(apply(data.subsets[[i]],1, mean))
}

rownames(data.meth.combined) <- rownames(meth.expr.main)
colnames(data.meth.combined) <- pathway.meth.names




######################################
## Principal componets of Pathway Gene Data
pc <- prcomp(data.gene.combined)
pc.pred <- predict(pc,newdata = data.gene.combined)
plot(pc.pred[,1], pc.pred[,2], pch = 19, col =c, main = "Gene Expression")
## Principal Components of Methylation Data
pc <- prcomp(data.meth.combined)
pc.pred <- predict(pc,newdata = data.meth.combined)
plot(pc.pred[,1], pc.pred[,2], pch = 19, col =c, main = "Methylations")

#################################################

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


D1 = ncol(Y1)
D2 = ncol(Y2)

Dmin = min(D1,D2)

Y1 <- Y1[,1:Dmin]
Y2 <- Y2[,1:Dmin]
D1 <- Dmin
D2 <- Dmin

################# Libraries ########################

library(MASS)
library(mixtools)
library(matrixcalc)
library(stats)
library(Runuran)
library(truncnorm)
library(Matrix)
library(MCMCpack)
library(psych)
library(VGAM)
library(MixSim)
library(statmod)
library(flexclust)
library(survcomp)
library(mixAK)
library(mclust)
library(monomvn)


setwd("/home/bit/ashar/Dropbox/Code/twoDPplusAFT")
############################# PARAMETERS for GIBB's SAMPLING ######################################
iter = 100
iter.burnin = 100
iter.thin  =5

## HYPER PRIORS
## Hyper parameters of the DP
shape.alpha <- 2
rate.alpha <- 1
# ## Hyperparameters for the GMM
# beta  = (D1 +D2)
# ro = 0.5
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
gmmx1$ro <- 0.5
gmmx1$beta <- D1

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
gmmx2$ro <- 0.5
gmmx2$beta <- D2


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
cont1 <- multiinit(Y1,c,gmmx1$beta, gmmx1$W, gmmx1$epsilon, gmmx1$ro, r, si,N,D1, sig2.dat)
gmmx1$mu <- cont1$mu
gmmx1$S <- cont1$S
regy1$lambda2 <- cont1$lambda2
regy1$tau2 <- cont1$tau2
regy1$betahat <- cont1$betahat
regy1$sigma2 <- cont1$sigma2
regy1$beta0 <- cont1$beta0

## For the second data set
source('multiinit.R')
cont2 <- multiinit(Y2,c,gmmx2$beta, gmmx2$W, gmmx2$epsilon, gmmx2$ro, r, si,N,D2, sig2.dat)
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


km <- multikmeansBlasso(c.init,Y1,Y2,D1,D2,That,K, r, si,sig2.dat,gmmx1.int, gmmx2.int, regy1.int, regy2.int,surv.obj )
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
  param <- posteriorGMMparametrs(c,Y1,gmmx1$mu,gmmx1$S, alpha, K, gmmx1$epsilon, gmmx1$W, gmmx1$beta, gmmx1$ro,N,D1 )
  gmmx1$mu <- param$mean
  gmmx1$S <- param$precision
  source('posteriorGMMparametrs.R')
  param2 <- posteriorGMMparametrs(c,Y2,gmmx2$mu,gmmx2$S, alpha,K, gmmx2$epsilon, gmmx2$W, gmmx2$beta, gmmx2$ro,N,D2 )
  gmmx2$mu <- param2$mean
  gmmx2$S <- param2$precision
  
  ########################## THE HYPERPARAMETERS OF THE GMM #################################  
  source('posteriorhyper.R')  
  # Updating the hyper paramters for the first data set
  hypercognate <- posteriorhyper (c, Y1, gmmx1$mu, gmmx1$S, gmmx1$epsilon, gmmx1$W, gmmx1$beta, gmmx1$ro,D1 )
  gmmx1$epsilon <- hypercognate$epsilon
  tmpW <- hypercognate$W
  gmmx1$W <- matrix(as.matrix(tmpW),nrow = D1, ncol =D1)
  gmmx1$ro <- hypercognate$ro
  
  ##Updating the hyper parameter for the second data set
  hypercognate2 <- posteriorhyper (c, Y2, gmmx2$mu, gmmx2$S, gmmx2$epsilon, gmmx2$W, gmmx2$beta, gmmx2$ro,D2 )
  gmmx2$epsilon <- hypercognate2$epsilon
  tmpW2 <- hypercognate2$W
  gmmx2$W <- matrix(as.matrix(tmpW2),nrow = D2, ncol =D2)
  gmmx2$ro <- hypercognate2$ro
  
  
  source('posteriorbeta.R')
  gmmx1$beta <- posteriorbeta(c, gmmx1$S, gmmx1$W)
  gmmx2$beta <- posteriorbeta(c, gmmx2$S, gmmx2$W)
  
  
  ############################# THE PARAMETERS of THE SURVIVAL TIME #########################################
  source('posteriortimeparameterspenalized.R')
  paramtime1 <- posteriortimeparameterspenalized(c,Y1, That, regy1$lambda2, regy1$tau2, regy1$sigma2, regy1$beta0, regy1$betahat, K, gmmx1$epsilon, gmmx1$W,  gmmx1$beta,gmmx1$ro,r, si, sig2.data,N, D1)
  regy1$beta0 <- paramtime1$beta0
  regy1$betahat <- paramtime1$betahat
  regy1$sigma2 <- paramtime1$sigma2
  regy1$lambda2 <- paramtime1$lambda2
  regy1$tau2 <- paramtime1$tau2
  
  paramtime2 <- posteriortimeparameterspenalized(c,Y2, That, regy2$lambda2, regy2$tau2, regy2$sigma2, regy2$beta0, regy2$betahat, K, gmmx2$epsilon, gmmx2$W, gmmx2$beta, gmmx2$ro, r, si, sig2.data,N, D2)
  regy2$beta0 <- paramtime2$beta0
  regy2$betahat <- paramtime2$betahat
  regy2$sigma2 <- paramtime2$sigma2
  regy2$lambda2 <- paramtime2$lambda2
  regy2$tau2 <- paramtime2$tau2
  
  
  
  
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


##################### GIBB'S SAMPLING PHASE #########################################
gmu.list <- list(0)
gbeta0.list <- list(0)
gbetahat.list <- list(0) 
gsigma2.list <- list(0)
glambda2.list <- list(0)
gtau2.list <- list(0)

mmu.list <- list(0)
mbeta0.list <- list(0)
mbetahat.list <- list(0) 
msigma2.list <- list(0)
mlambda2.list <- list(0)
mtau2.list <- list(0)





c.list <- list(0)
That.list <- list(0)
count <- 1
o.init =o
#################### BURNIN PHASE ###################################################
print("SAMPLING...PHASE")
for (o in o.init:iter) {
  
  
  ################## PARAMETERS OF THE DP Mixture Model ######################################################
  ## Updating the parameters based on the observations 
  source('posteriorGMMparametrs.R')
  param <- posteriorGMMparametrs(c,Y1,gmmx1$mu,gmmx1$S, alpha, K, gmmx1$epsilon, gmmx1$W, gmmx1$beta, gmmx1$ro,N,D1 )
  gmmx1$mu <- param$mean
  gmmx1$S <- param$precision
  source('posteriorGMMparametrs.R')
  param2 <- posteriorGMMparametrs(c,Y2,gmmx2$mu,gmmx2$S, alpha,K, gmmx2$epsilon, gmmx2$W, gmmx2$beta, gmmx2$ro,N,D2 )
  gmmx2$mu <- param2$mean
  gmmx2$S <- param2$precision
  
  ########################## THE HYPERPARAMETERS OF THE GMM #################################  
  source('posteriorhyper.R')  
  # Updating the hyper paramters for the first data set
  hypercognate <- posteriorhyper (c, Y1, gmmx1$mu, gmmx1$S, gmmx1$epsilon, gmmx1$W, gmmx1$beta, gmmx1$ro,D1 )
  gmmx1$epsilon <- hypercognate$epsilon
  tmpW <- hypercognate$W
  gmmx1$W <- matrix(as.matrix(tmpW),nrow = D1, ncol =D1)
  gmmx1$ro <- hypercognate$ro
  
  ##Updating the hyper parameter for the second data set
  hypercognate2 <- posteriorhyper (c, Y2, gmmx2$mu, gmmx2$S, gmmx2$epsilon, gmmx2$W, gmmx2$beta, gmmx2$ro,D2 )
  gmmx2$epsilon <- hypercognate2$epsilon
  tmpW2 <- hypercognate2$W
  gmmx2$W <- matrix(as.matrix(tmpW2),nrow = D2, ncol =D2)
  gmmx2$ro <- hypercognate2$ro
  
  
  source('posteriorbeta.R')
  gmmx1$beta <- posteriorbeta(c, gmmx1$S, gmmx1$W)
  gmmx2$beta <- posteriorbeta(c, gmmx2$S, gmmx2$W)
  
  
  ############################# THE PARAMETERS of THE SURVIVAL TIME #########################################
  source('posteriortimeparameterspenalized.R')
  paramtime1 <- posteriortimeparameterspenalized(c,Y1, That, regy1$lambda2, regy1$tau2, regy1$sigma2, regy1$beta0, regy1$betahat, K, gmmx1$epsilon, gmmx1$W,  gmmx1$beta,gmmx1$ro,r, si, sig2.data,N, D1)
  regy1$beta0 <- paramtime1$beta0
  regy1$betahat <- paramtime1$betahat
  regy1$sigma2 <- paramtime1$sigma2
  regy1$lambda2 <- paramtime1$lambda2
  regy1$tau2 <- paramtime1$tau2
  
  paramtime2 <- posteriortimeparameterspenalized(c,Y2, That, regy2$lambda2, regy2$tau2, regy2$sigma2, regy2$beta0, regy2$betahat, K, gmmx2$epsilon, gmmx2$W, gmmx2$beta, gmmx2$ro, r, si, sig2.data,N, D2)
  regy2$beta0 <- paramtime2$beta0
  regy2$betahat <- paramtime2$betahat
  regy2$sigma2 <- paramtime2$sigma2
  regy2$lambda2 <- paramtime2$lambda2
  regy2$tau2 <- paramtime2$tau2
  
  
  
  
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
  
  if(o%% thin == 0 ){
    gmu.list[[count]] <- gmmx1$mu
    gbeta0.list[[count]] <- regy1$beta0
    gbetahat.list[[count]] <- regy1$betahat  
    gsigma2.list[[count]] <- regy1$sigma2
    glambda2.list[[count]] <- regy1$lambda2
    gtau2.list[[count]] <- regy1$tau2
    
    mmu.list[[count]] <- gmmx2$mu
    mbeta0.list[[count]] <- regy2$beta0
    mbetahat.list[[count]] <- regy2$betahat  
    msigma2.list[[count]] <- regy2$sigma2
    mlambda2.list[[count]] <- regy2$lambda2
    mtau2.list[[count]] <- regy2$tau2
    
    c.list[[count]] <- c
    count <- count +1
  }
  
  
  
  
  
  
} 

################# OBTAINING POSTERIOR SAMPLES ########################################

count <- count -1

#### Taking the average values for the results #########################################

c.matrix <- matrix(NA, nrow = N, ncol = count)
for ( i in 1:count){
  c.matrix[,i] <- c.list[[i]]
}
c.final <- apply(c.matrix,1,median)

################## OBTAINING POSTERIOR BETA ESTIMATES #################################


glist.betahat <- list(0)
mlist.betahat <- list(0)
for ( i in 1:count){
  glist.betahat[[i]] <- (gbetahat.list[[i]][1:4,] != 0) +0
}

for ( i in 1:count){
  mlist.betahat[[i]] <- (mbetahat.list[[i]][1:4,] != 0) +0
}

### Probability of betahat of genes
final.betahat <- rbind(final.betahat1, final.betahat2, final.betahat3, final.betahat4)
rownames(final.betahat) = c("cluster_1","cluster_2","cluster_3","cluster_4")
yy <- as.list(KEGGPATHID2NAME)
colnames(final.betahat) <- as.character(unlist(yy[colnames(Y)])) 






#### Subset the betahat such that only non zero elements

subset.betahat <- which(apply(final.betahat,2,mean)!=0 )

## A new Data frame with only these pathways
Y.new <- Y[,subset.betahat]
colnames(Y.new) <- names(subset.betahat)
rownames(Y.new) <- c.final

ind1 <- which(c.final ==1)
ind2 <- which(c.final == 2)
ind3 <- which(c.final ==3)
ind4 <- which(c.final ==4)

Y.order <- matrix(0, nrow = N, ncol = ncol(Y.new))
Y.order <- rbind(Y.new[ind1,],Y.new[ind2,],Y.new[ind3,],Y.new[ind4,])


pdf("/home/bit/ashar/ExpressionSets/Verhark/PathwayHeatmap.pdf")
heatmap.2(t(Y.order), dendrogram = "none",col = , Rowv = FALSE, Colv = FALSE, ColSideColors = rownames(Y.order), cexRow =0.6, margins = c(5,12), trace = "none", main = "Pathway expression \n Verhaak data") 
dev.off()



heatmapdata <- as.data.frame(final.betahat)
pdf("/home/bit/ashar/ExpressionSets/Verhark/heatmap.pdf")
heatmap.2(t(as.matrix(heatmapdata)),dendrogram="none", col =cm.colors(90), margins=c(6,10), main = "Marginal Posterior prob. ", cexCol = 0.85, cexRow = 0.7, ylab ="KEGG pathways", trace = "none")
dev.off()


#### BAR PLOTS ##############
pdf("/home/bit/ashar/ExpressionSets/Verhark/Barplot.pdf")
par(mfrow=c(2,2))
barplot(final.betahat[1,], main = "Cluster1")
barplot(final.betahat[2,], main = "Cluster2")
barplot(final.betahat[3,], main = "Cluster3")
barplot(final.betahat[4,], main = "Cluster4")
dev.off()


surv.fit <- survfit(surv.ob ~ c.final)
logrank <- survdiff(surv.ob ~ c.final)
logrank_verhaak <- survdiff(surv.ob ~ c.true)






############### Drawing Some Plots ###################################################
t.grand <- exp(That)
surv.days <- Surv((t.grand),censoring)
fit <- survfit(surv.days ~ c.final)
pdf("TCGA_GBM_Integ.pdf")
par(mfrow=c(1,2))
pc <- prcomp(data.gene.combined)
pc.pred <- predict(pc,newdata = data.gene.combined)
plot(pc.pred[,1], pc.pred[,2], pch = 19, col =c.final, main = "Gene Expression")
## Principal Components of Methylation Data
pc <- prcomp(data.meth.combined)
pc.pred <- predict(pc,newdata = data.meth.combined)
plot(pc.pred[,1], pc.pred[,2], pch = 19, col =c.final, main = "Methylation")
dev.off()

pdf("TCGA_GBM_Survival.pdf")
plot(fit, col = c("black","red","green3","blue"), xlab = "Blue: 967d(23) \n Black 444d(116) ", main = "Survival Curves, N =194 \n Logrank pvalue 9.320254e-05")
dev.off()

