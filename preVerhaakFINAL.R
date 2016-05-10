############################################################################################
###################### This File Prepares the TCGA GBM data set for analysis #####################
##############################################################################################
rm(list = ls())
load("/home/bit/ashar/ownCloud/DPMM_RESULTS/ONE_VIEW/Verhark/Final/DataVerhaak.RData")
load("/home/bit/ashar/ExpressionSets/TWO_VIEW/TCGA_GBM/IntegratedDataSurv.rda")
load("/home/bit/ashar/ExpressionSets/TWO_VIEW/TCGA_GBM/mRNA.rda")
load("/home/bit/ashar/ExpressionSets/TWO_VIEW/TCGA_GBM/miRNA.rda")


pheno.train <- relev$pheno.train
pheno.test <- relev$pheno.test
signature.dpmm <- relev$signature.dpmm

patient.verhaak.train <- colnames(relev$Y.train.prelim)
patient.verhaak.test <- colnames(relev$Y.test.prelim)

patient.tcga <- rownames(mRNA)[rownames(mRNA) %in% rownames(miRNA)]

#### Checking how many training patients are common ##########
sum((patient.verhaak.train %in% patient.tcga)+0)


#### Checking how many testing patients are common ##########
sum((patient.verhaak.test %in% patient.tcga)+0)


####### TRAINING DATA ########################################
###############################################################
Y.gene.verk <- t(relev$Y.train.prelim[,patient.verhaak.train %in% patient.tcga])
Y.gene.tcga <- mRNA[patient.tcga[patient.tcga %in% patient.verhaak.train],]
Y.mirna.tcga <- miRNA[patient.tcga[patient.tcga %in% patient.verhaak.train],]

###### TESTING DATA ##########################################
###############################################################
Y.gene.verk.test <- t(relev$Y.test.prelim[,patient.verhaak.test %in% patient.tcga])
Y.gene.tcga.test <- mRNA[patient.tcga[patient.tcga %in% patient.verhaak.test],]
Y.mirna.tcga.test <- miRNA[patient.tcga[patient.tcga %in% patient.verhaak.test],]


######################### Just checking if the row names are the same ################################
rownames(Y.gene.tcga) == rownames(Y.mirna.tcga)
rownames(Y.gene.tcga.test) == rownames(Y.mirna.tcga.test)


########## LET US CHECK IF THE SAME INFORMATION IS CONTAINED IN VERHAAK AND TCGA DATA SET ############
order.verk <- match(rownames(Y.gene.tcga), rownames(Y.gene.verk))
Y.gene.verk2 <- Y.gene.verk[order.verk,]
rownames(Y.gene.verk2) == rownames(Y.gene.tcga)


########## LET US CHECK IF THE SAME INFORMATION IS CONTAINED IN VERHAAK AND TCGA DATA SET ############
order.verk.test <- match(rownames(Y.gene.tcga.test), rownames(Y.gene.verk.test))
Y.gene.verk2.test <- Y.gene.verk.test[order.verk.test,]
rownames(Y.gene.verk2.test) == rownames(Y.gene.tcga.test)


######## Lets Get the Pheno Data #####################################################################
pheno.train2 <- pheno.train[match(rownames(Y.gene.verk2),pheno.train$Patient),]
pheno.test2 <- pheno.test[match(rownames(Y.gene.verk2.test),pheno.test$Patient),]


#### Just Checking about the order of the phenotype #################################################
pheno.train2$Patient == rownames(Y.mirna.tcga)
pheno.test2$Patient == rownames(Y.mirna.tcga.test)


#################################################################################################################################################################
############################ PRE - FILTERED FEATURE SELECTION ####################################################################################################

##### GENE EXPRESSION #################
######## Getting survival times and status #####################
Time <- as.numeric(exp(pheno.train2[,3]))
status <- pheno.train2[,2]
surv.obj <- Surv(Time, status)

##### FITTING A UNIVARIATE COX REGRESSION MODEL ################################
pvalue.sig <- c(0)
pvalue.adj <- c(0)
X_train <- as.data.frame(Y.gene.tcga)

for ( i in 1:ncol(X_train)){
  q <- unlist(summary(coxph(surv.obj ~ X_train[,i], data = X_train)))
  pvalue.sig[i] <- q$logtest.pvalue  
}
pvalue.adj <- p.adjust(pvalue.sig, method = "fdr")
signature.dpmm.gene <- colnames(X_train)[(pvalue.sig < 0.005)]

###### MIRNA EXPRESSION #######################################################
##### FITTING A UNIVARIATE COX REGRESSION MODEL ################################
pvalue.sig <- c(0)
pvalue.adj <- c(0)
X_train <- as.data.frame(Y.mirna.tcga)

for ( i in 1:ncol(X_train)){
  q <- unlist(summary(coxph(surv.obj ~ X_train[,i], data = X_train)))
  pvalue.sig[i] <- q$logtest.pvalue  
}
pvalue.adj <- p.adjust(pvalue.sig, method = "fdr")
signature.dpmm.mirna <- colnames(X_train)[(pvalue.sig < 0.1)]


############### MAKING TRAINING AND TESTING DATA SETS ############
Y1 <- Y.gene.tcga[,signature.dpmm.gene]
Y2 <- Y.mirna.tcga[,signature.dpmm.mirna]
pheno <-  pheno.train2

Y1.test <- Y.gene.tcga.test[,signature.dpmm.gene]
Y2.test <- Y.mirna.tcga.test[,signature.dpmm.mirna]
pheno.test <- pheno.test2

#################### SAVING THE DATA SETS ################################################

relev <- list('Y1' =Y1, 'Y2' = Y2,'pheno' = pheno, 'Y1.test' = Y1.test, 'Y2.test' = Y2.test, 'pheno.test' = pheno.test)
save(relev, file = 'DataTCGA-GBM.RData')
