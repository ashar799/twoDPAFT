###########################################################################################
### This applies TCGA data set for GBM on our DPMM model ###################################
############################################################################################
#### Copies from Simulation_Main from the one View Case ####################################
############################################################################################
rm(list = ls())

#### Load Data #########
load("/home/bit/ashar/ExpressionSets/TWO_VIEW/TCGA_GBM/FINAL/DataTCGA-GBM.RData")

Y1 <- relev$Y1
Y2 <- relev$Y2

Y1.test <- relev$Y1.test
Y2.test <- relev$Y2.test

pheno <- relev$pheno
pheno.test <- relev$pheno.test


## Number of points
N.train =  nrow(Y1)
N.test = nrow(Y1.test)

N <- N.train
## Number of Clusters
F = 4

######
D1 <- ncol(Y1)
D2 <- ncol(Y2)

####
time <- pheno$Survival
censoring <- pheno$Censoring

time.new <- pheno.test$Survival
censoring.new <- pheno.test$Censoring

###########
c.true <- pheno$Subtype
levels(c.true) <- c(1,2,3,4)

##########
c.true.new <- pheno.test$Subtype
levels(c.true.new) <- c(1,2,3,4)


############################# PARAMETERS for GIBB's SAMPLING ####
iter = 50
iter.burnin = 50
iter.thin  = 2
k = 4
K <-  as.integer(N)
Time <- cbind(time,censoring)


######################### Initialize the Parameters ################
source('initializemultiDPMM.R')
initializemultiDPMM()


########### Train the Model #########################################
source('burninmultiDPMM.R')
burninmultiDPMM()

########### Train the Model #########################################
source('gibbsmultiDPMM.R')
gibbsmultiDPMM()

##### Analyzing the Model #########################################
source('analyzemultiDPMM.R')
analyzemultiDPMM()


###### Log Rank Statistic for the Testing Data SEt ##############
####Fitting Survival Curves
surv.ob <- Surv(exp(time),censoring)
surv.fit <- survfit(surv.ob ~ c.final)
logrank <- survdiff(surv.ob ~ c.final)


#### Predicting Class ###########################################
source('predictmultiCLASS.R')
predictmultiCLASS(Y1.test, Y2.test, time.new,censoring.new)
## Check how much concordance is there
test.randindex <- adjustedRandIndex(apply(posteriorprob,1,which.max),c.true.new)
lr <- c(0)
for (j in 1:Nps){
  lr[j] <-  1 - pchisq(unlist(survdiff(surv.ob.new ~ c.new.list[[j]]))$chisq,df = length(table(c.new.list[[j]])) -1 )
}
### As the one which has the lowest p-value does not have 4 clusters
### I will take the onw which has 4 clusters and has kind of OK p-value
c.final.new <- c.new.list[[19]]



###### Predicting Survival Times ####################################
source('multipredictchineseAFTtime.R')
multipredictchineseAFTtime(Y1.test, Y2.test)
predicted.cindex <- survConcordance(Surv(exp(time.new),censoring.new) ~ exp(-post.time.avg))[1]


###### Log Rank Statistic for the Testing Data SEt ##############
####Fitting Survival Curves
surv.ob <- Surv(exp(time.new),censoring.new)
surv.fit <- survfit(surv.ob ~ c.final.new)
logrank <- survdiff(surv.ob ~ c.final.new)


##############################################################################################################################################
########################### VISUALIZATION ####################################################################################################
##############################################################################################################################################
##### Generating some plots #####################################
### Training Data Set ###########################################


pc1 <- prcomp(Y1)
pc.pred1 <- predict(pc1,newdata = Y1)
cb <- c("Orange","Green","Blue","Red")
p1 <- ggplot(as.data.frame(pc.pred1), aes(x=pc.pred1[,1], y= pc.pred1[,2], colour= c("Worst","GoodMod.","BadMod.","Best")[c.final])) + ggtitle(" TCGA GBM Gene Expression \n DPMM Clustering \n 45% Variance Explained") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") + scale_colour_manual(values=cb)

pc2 <- prcomp(Y2)
pc.pred2 <- predict(pc2,newdata = Y2)
p2 <- ggplot(as.data.frame(pc.pred2), aes(x=pc.pred2[,1], y= pc.pred2[,2], colour= c("Worst","GoodMod.","BadMod.","Best")[c.final])) + ggtitle(" TCGA GBM miRNA Expression \n DPMM Clustering \n 55% Variance Explained") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") + scale_colour_manual(values=cb)

surv.ob <- Surv(exp(time),censoring)
surv.fit <- survfit(surv.ob ~ c.final)
logrank <- survdiff(surv.ob ~ c.final)
p3 <- ggsurv(surv.fit, main = " TCGA -GBM DPMM \n Kaplan Meier Estimators \n p-value 6e-04 ") + ggplot2::guides(linetype = FALSE) + ggplot2::scale_colour_discrete(name = 'Classes',breaks = c(1,2,3,4),labels = c("Worst","GoodMod.","BadMod.","Best")) 

##########################################################
########## Creating a Ranking for the Points ####################
###### Creating a Ratio for the Points ########################
rank <- matrix(0, nrow = N, ncol =Nps)

for (j in 1:Nps){
  
  Y1.scaled <- matrix(0, nrow = N, ncol =D1)
  for ( v in 1:4){
    clust <- which(c.list[[j]] == v)
    Y1.scaled[clust,1:D1] <- scale(Y1[clust,1:D1], center = TRUE, scale = TRUE)
  }
  
  
  for ( i in 1:N){
    rank[i,j] <- dMVN(as.vector(t(Y1[i,1:D1])), mean = est.gmmx1[[j]]$mu[c.final[i],1:D1], Q= est.gmmx1[[j]]$S[c.final[i],1:D1,1:D1], log = TRUE) - dMVN(as.vector(t(Y1[i,1:D1])), mean = est.gmmx1[[j]]$mu[1,1:D1], Q= est.gmmx1[[j]]$S[1,1:D1,1:D1], log = TRUE)+  dnorm(x = That[i], mean = est.regy1[[j]]$beta0[c.final[i]] + est.regy1[[j]]$betahat[c.final[i],1:D1] %*% as.vector(t(Y1.scaled[i,1:D1])), sd = sqrt(est.regy1[[j]]$sigma2[c.final[i]]), log =TRUE) -  dnorm(x = That[i], mean = est.regy1[[j]]$betahat[1,1:D1] %*% as.vector(t(Y1.scaled[i,1:D1])), sd = sqrt(est.regy1[[j]]$sigma2[1]), log =TRUE) 
  }
}

avg.rank <- apply(rank,1,mean)
order.zo <- range01(avg.rank)

order.train <- sort(order.zo,index.return = TRUE, decreasing = TRUE)
Y1.order <- Y1[order.train$ix,]
c.final.order <- c.final[order.train$ix]


######## Reordering Again ################
order.2 <- c(which(c.final.order==1),which(c.final.order==3),which(c.final.order==2),which(c.final.order==4))
Y1.order.2 <- Y1.order[order.2,]
c.final.order.2 <- c.final.order[order.2]


##### Plotting the Heatmap #################
hmcols<-colorRampPalette(c("green","black","red"))(256)
heatmap.2(x = t(Y1.order.2), Rowv =TRUE ,Colv = FALSE, dendrogram = c("row"), scale ="row",ColSideColors =c("red","blue","orange","green")[c.final.order.2], labRow = colnames(Y1.order.2), labCol = NA, main = ' \n Training Set \n TCGA GBM Gene Expression Data', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")
legend("topright", legend = c("Worst","Good Moderate","Bad moderate","Best"),fill = c("Red","Blue","Orange","Green"), cex = 0.4)



#### Plotting mi-RNA expression ##############

rank2 <- matrix(0, nrow = N, ncol =Nps)

for (j in 1:Nps){
  
  Y2.scaled <- matrix(0, nrow = N, ncol =D2)
  for ( v in 1:4){
    clust <- which(c.list[[j]] == v)
    Y2.scaled[clust,1:D2] <- scale(Y2[clust,1:D2], center = TRUE, scale = TRUE)
  }
  
  
  for ( i in 1:N){
    rank2[i,j] <- dMVN(as.vector(t(Y2[i,1:D2])), mean = est.gmmx2[[j]]$mu[c.final[i],1:D2], Q= est.gmmx2[[j]]$S[c.final[i],1:D2,1:D2], log = TRUE) - dMVN(as.vector(t(Y2[i,1:D2])), mean = est.gmmx2[[j]]$mu[1,1:D2], Q= est.gmmx2[[j]]$S[1,1:D2,1:D2], log = TRUE)+  dnorm(x = That[i], mean = est.regy2[[j]]$beta0[c.final[i]] + est.regy2[[j]]$betahat[c.final[i],1:D2] %*% as.vector(t(Y2.scaled[i,1:D2])), sd = sqrt(est.regy2[[j]]$sigma2[c.final[i]]), log =TRUE) -  dnorm(x = That[i], mean = est.regy2[[j]]$betahat[1,1:D2] %*% as.vector(t(Y2.scaled[i,1:D2])), sd = sqrt(est.regy2[[j]]$sigma2[1]), log =TRUE) 
  }
}

avg.rank2 <- apply(rank2,1,mean)
order.zo2 <- range01(avg.rank2)

order.train2 <- sort(order.zo2,index.return = TRUE, decreasing = TRUE)
Y2.order <- Y2[order.train2$ix,]
c.final.order2 <- c.final[order.train2$ix]


######## Reordering Again ################
order.2 <- c(which(c.final.order2==1),which(c.final.order2==3),which(c.final.order2==2),which(c.final.order2==4))
Y2.order.2 <- Y2.order[order.2,]
c.final.order.22 <- c.final.order2[order.2]


##### Plotting the Heatmap #################
hmcols<-colorRampPalette(c("green","black","red"))(256)
heatmap.2(x = t(Y2.order.2), Rowv =TRUE ,Colv = FALSE, dendrogram = c("row"), scale ="row",ColSideColors =c("red","blue","orange","green")[c.final.order.22], labRow = colnames(Y2.order.2), labCol = NA, main = ' \n Training Set \n TCGA GBM mi-RNA Expression Data', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")
legend("topright", legend = c("Worst","Good Moderate","Bad moderate","Best"),fill = c("Red","Blue","Orange","Green"), cex = 0.4)


### Plotting Everything together #####
pdf('TCGA-GBM-Training.pdf')
hmcols<-colorRampPalette(c("green","black","red"))(256)
heatmap.2(x = t(Y1.order.2), Rowv =TRUE ,Colv = FALSE, dendrogram = c("row"), scale ="row",ColSideColors =c("red","blue","orange","green")[c.final.order.2], labRow = colnames(Y1.order.2), labCol = NA, main = ' \n Training Set \n TCGA GBM Gene Expression Data', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")
legend("topright", legend = c("Worst","Good Moderate","Bad moderate","Best"),fill = c("Red","Blue","Orange","Green"), cex = 0.4)
p1
heatmap.2(x = t(Y2.order.2), Rowv =TRUE ,Colv = FALSE, dendrogram = c("row"), scale ="row",ColSideColors =c("red","blue","orange","green")[c.final.order.2], labRow = colnames(Y2.order.2), labCol = NA, main = ' \n Training Set \n TCGA GBM mi-RNA Expression Data', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")
legend("topright", legend = c("Worst","Good Moderate","Bad moderate","Best"),fill = c("Red","Blue","Orange","Green"), cex = 0.4)
p2
p3
dev.off()


#### Analysing the signature ###############################
colnames(heatmapdata1) <- colnames(Y1)
colnames(heatmapdata2) <- colnames(Y2)
rownames(heatmapdata1) <- c("Worst","GoodModerate","BadModerate","Best")
rownames(heatmapdata2) <- c("Worst","GoodModerate","BadModerate","Best")

hmcols<-colorRampPalette(c("white","black"))(128)
pdf('Signature-TCGA.pdf')
par(mfrow=c(1,2))
heatmap.2(t(as.matrix(heatmapdata1))[,c(2:3)], Rowv =FALSE ,Colv = FALSE, col = hmcols, margins=c(6,10), main = "Posterior prob. for Selection \n  TCGA-Gene  Signature ", cexCol = 0.85, cexRow = 0.7)
heatmap.2(t(as.matrix(heatmapdata2))[,c(2:3)],Rowv = FALSE ,Colv = FALSE, col = hmcols, margins=c(6,10), main = "Posterior prob. for Selection \n TCGA-miRNA Signature ", cexCol = 0.85, cexRow = 0.7)
dev.off()

#############################################################################################
########################## REFERENCE VALUES FOR THE TABLE ##################################
############################################################################################
##############################################################################################
####################### FACTS ABOUT TESTING DATA WITH K-MEANS #################################
Y.combined <-  cbind(Y1,Y2)
Y.test.combined <- cbind(Y1.test,Y2.test)
km.labels <- kmeans(Y.combined, centers =k, nstart =10)$cluster
knn.labels <- knn(train = Y.combined, test = Y.test.combined, cl = km.labels, k = 4)
label.train <- km.labels
label.test <- knn.labels

surv <- Surv(exp(time),censoring)
survdiff(surv ~ km.labels)
surv.ob <- Surv(exp(time.new),censoring.new)
survdiff(surv.ob ~ knn.labels)

####### COMPARISON METHOD Their Labels for training + KNN (Penalized Cox) ###############################
linear.kkpcox.recovery <- c(0)
linear.kkpcox.prediction <- c(0)
for ( q in 1:F){
  ind <- which(label.train == q)
  ind.new <- which(label.test == q)
  reg.cox <- cv.glmnet(x = as.matrix(Y.combined[ind,]), y = Surv(exp(time[ind]),censoring[ind]), family ="cox")
  linear.kkpcox.recovery[ind] <- predict(object =reg.cox, newx = as.matrix(Y.combined[ind,]), s = 'lambda.1se')
  linear.kkpcox.prediction[ind.new] <- predict(object =reg.cox, newx = as.matrix(Y.test.combined[ind.new,]), s = 'lambda.1se')
}
recovCIndex.kkpcox <- as.numeric(survConcordance(Surv(exp(time), censoring) ~ linear.kkpcox.recovery)[1])
predCIndex.kkpcox <- as.numeric(survConcordance(Surv(exp(time.new), censoring.new) ~ linear.kkpcox.prediction)[1])

