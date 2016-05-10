### This model takes P Cells and Z cells microarray data and used DPMM to jointly cluster them
##### 38 Patients ############################################################################
#### Copies from Simulation_Main from the one View Case

load("/home/bit/ashar/ExpressionSets/TWO_VIEW/VIP/DataPandZCombined.RData")
Y1 <- relev$Y1
Y2 <- relev$Y2
time <- relev$time
N <- nrow(Y1)
D1 <- ncol(Y1)
D2 <- ncol(Y2)


censoring <- rep(1,N)



############################# PARAMETERS for GIBB's SAMPLING ####
iter = 100
iter.burnin = 50
iter.thin  = 5
k = 2
F=2
K <-  as.integer(N/2)
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

#############################################################################################
############ Generating some Plots ##########################
pc <- prcomp(Y1)
pc.pred <- predict(pc,newdata = Y1)
p1 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.final))) + ggtitle(" i-DPMM Clustering \n VIP patients (38) \n 60 Gene DPMM signature \n Center") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 

pc2 <- prcomp(Y2)
pc2.pred <- predict(pc2,newdata = Y2)
p2 <- ggplot(as.data.frame(pc2.pred), aes(x=pc2.pred[,1], y= pc2.pred[,2], colour= as.factor(c.final))) + ggtitle(" i-DPMM Clustering \n VIP patients (38) \n 60 Gene DPMM signature \n Periphery") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 

####Fitting Survival Curves
surv.ob <- Surv(exp(time),censoring)
surv.fit <- survfit(surv.ob ~ c.final)
logrank <- survdiff(surv.ob ~ c.final)

pdf("IntegratedVIPsurvival-curves.pdf")
plot(surv.fit, main = 'PF Survival Curves for VIP Patients \n p = 0.01 \n Clustering Based on P and Z Samples', xlab= 'Mean PFS Cluster 1: 285d \n Mean PFS Cluster2: 625d')
dev.off()



####### Calculating the rank
rank <- matrix(0, nrow = N, ncol =Nps)

for (j in 1:Nps){
  for ( i in 1:N){
    rank[i,j] <- dMVN(as.vector(t(Y1[i,1:D1])), mean = est.gmmx1[[j]]$mu[1,1:D1], Q= est.gmmx1[[j]]$S[1,1:D1,1:D1], log = TRUE) -  dMVN(as.vector(t(Y1[i,1:D1])), mean = est.gmmx1[[j]]$mu[2,1:D1], Q= est.gmmx1[[j]]$S[2,1:D1,1:D1], log = TRUE)
  }
}

avg.rank <- apply(rank,1,mean)
order.zo <- range01(avg.rank)

order.train <- sort(order.zo,index.return = TRUE, decreasing = TRUE)
Y1.order <- Y1[order.train$ix,]
c1.final.order <- c.final[order.train$ix]
rownames(Y1.order) <- rownames(Y1)[order.train$ix]

c1.ultimate <- as.factor(c1.final.order)
levels(c1.ultimate) <- c(1,2)

####### Calculating the rank
rank <- matrix(0, nrow = N, ncol =Nps)

for (j in 1:Nps){
  for ( i in 1:N){
    rank[i,j] <- dMVN(as.vector(t(Y2[i,1:D2])), mean = est.gmmx2[[j]]$mu[1,1:D2], Q= est.gmmx2[[j]]$S[1,1:D2,1:D2], log = TRUE) -  dMVN(as.vector(t(Y2[i,1:D2])), mean = est.gmmx2[[j]]$mu[2,1:D2], Q= est.gmmx2[[j]]$S[2,1:D2,1:D2], log = TRUE)
  }
}


avg.rank <- apply(rank,1,mean)
order.zo <- range01(avg.rank)

order.train <- sort(order.zo,index.return = TRUE, decreasing = TRUE)
Y2.order <- Y2[order.train$ix,]
c2.final.order <- c.final[order.train$ix]
rownames(Y2.order) <- rownames(Y2)[order.train$ix]

c2.ultimate <- as.factor(c2.final.order)
levels(c2.ultimate) <- c(1,2)


pdf('iDPMM.pdf')
hmcols<-colorRampPalette(c("green","black","red"))(256)
heatmap.2(x = t(Y1.order), Rowv =TRUE ,Colv = FALSE, dendrogram = c("row"), scale ="row",ColSideColors =c("red","black")[c1.ultimate], labRow = colnames(Y1), labCol = rownames(Y1.order), main = 'i DPMM clustering on VIP samples \n Microarray Center Cells \n 60 Gene DPMM signature', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")
legend("topright", legend = c("Cluster 1","Cluster 2"),fill = c("Red","Black"), cex = 0.4)
p1
heatmap.2(x = t(Y2.order), Rowv =TRUE ,Colv = FALSE, dendrogram = c("row"), scale ="row",ColSideColors =c("red","black")[c2.ultimate], labRow = colnames(Y2), labCol = rownames(Y2.order), main = 'i-DPMM clustering on VIP samples \n Microarray peripheral Cells \n 60 Gene DPMM signature', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")
legend("topright", legend = c("Cluster 1","Cluster 2"),fill = c("Red","Black"), cex = 0.4)
p2
dev.off()






