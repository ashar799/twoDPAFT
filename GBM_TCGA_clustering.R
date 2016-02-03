###### THIS FILE CHECKS SOME CLUSTERING FEATURES IN THE DATA ##########

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

main.mirna.index <- match(colnames(miRNA.expr),colnames(rownames(gene.expr.main)))
main.mirna.index <- main.mirna.index[!is.na(main.mirna.index)]
mirna.exprs.main <- t(miRNA.expr)[main.meth.index,]

  
  
###JUST CHECKING
rownames(gene.expr.main) == rownames(meth.expr.main)

### Convert Gene Expression and Methylation to Entrez ID levels
## Gene Expression Data
xx <- as.list(org.Hs.egALIAS2EG)
# Remove pathway identifiers that do not map to any entrez gene id
xx <- xx[!is.na(xx)]

xx2 <- IlluminaHumanMethylation450kENTREZID
mapped_probes <- mappedkeys(xx2)
# Convert to a list
xx2 <- as.list(xx2[mapped_probes])


### Calculate the variances of the features of Gene Expression
sd.gene <- apply(gene.expr.main,2,sd)
quantile <- quantile(sd.gene,probs = seq(0, 1, 0.05))
sum(sd.gene > quantile[20] + 0)
Y.gene <- gene.expr.main[,sd.gene > quantile[20]]


### For Methylation
sd.meth <- apply(meth.expr.main,2,sd)
quantile.meth <- quantile(sd.meth,probs = seq(0, 1, 0.02))
sum(sd.meth > quantile.meth[50] + 0)
Y.meth <- meth.expr.main[,sd.meth > quantile.meth[50]]

Y.mirna <- mirna.exprs.main

#################################################################################################
#################################################################################################
############ INDIVIDUAL CLUSTERING ##########################################################################
###############################################################################################
k.gene <- kmeans(Y.gene,centers =3, nstart = 10)
pc <- prcomp(Y.gene)
pc.pred <- predict(pc,newdata = Y.gene)
plot(pc.pred[,1], pc.pred[,2], pch = 19, col = k.gene$cluster)

k.meth <- kmeans(Y.meth,centers =3, nstart = 10)
pc <- prcomp(Y.meth)
pc.pred <- predict(pc,newdata = Y.meth)
plot(pc.pred[,1], pc.pred[,2], pch = 19, col = k.meth$cluster)

k.mirna <- kmeans(Y.mirna,centers =3, nstart = 10)
pc <- prcomp(Y.mirna)
pc.pred <- predict(pc,newdata = Y.mirna)
plot(pc.pred[,1], pc.pred[,2], pch = 19, col = k.mirna$cluster)

### Check the similarity of the clustering
require(mclust)
adjustedRandIndex(k.gene$cluster,k.meth$cluster)
adjustedRandIndex(k.gene$cluster,k.mirna$cluster)
adjustedRandIndex(k.meth$cluster,k.mirna$cluster)

### The most similar Clustering is from the Methylation and Gene Expression

#### Use ICluster #############################################
library(iCluster)
data <- list("gene" = Y.gene, "meth" = Y.meth, "miRNA" = Y.mirna)
#tune.iCluster2(data, k = 4)
#fit <- iCluster(datasets = data, k =3, lambda, scalar=FALSE, max.iter=50,epsilon=1e-3)

#### Let's just sparse CCA to analyze
library(PMA)
fit.cca <- MultiCCA.permute(xlist = data)
cca <- MultiCCA(data, penalty = fit.cca$bestpenalties, ncomponents = 5)
##plot the new points
Y.gene.cca <-  (Y.gene)%*% cca$ws[[1]] 
Y.meth.cca <-  (Y.meth)%*% cca$ws[[2]] 
Y.mirna.cca <-  (Y.mirna)%*% cca$ws[[3]] 

plot(Y.gene.cca[,1], Y.gene.cca[,2], pch = 19, col = k.meth.cca$cluster )
plot(Y.meth.cca[,1], Y.meth.cca[,2], pch = 19,col = k.meth.cca$cluster )
plot(Y.mirna.cca[,1], Y.mirna.cca[,2], pch = 19, col = k.meth.cca$cluster )

k.gene <- kmeans(Y.gene,centers =2, nstart = 10)
k.gene.cca <- kmeans(Y.gene.cca,centers =2, nstart = 10)
adjustedRandIndex(k.gene$cluster,k.gene.cca$cluster)

k.meth <- kmeans(Y.meth,centers =2, nstart = 10)
k.meth.cca <- kmeans(Y.meth.cca,centers =2, nstart = 10)
adjustedRandIndex(k.meth$cluster,k.meth.cca$cluster)

