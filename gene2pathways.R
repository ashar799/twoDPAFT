library(globaltest)
library(org.Hs.eg.db)
library(KEGG.db)
library(IlluminaHumanMethylation450k.db)

###### This script calculates the Pathway level scores for both Gene Expression Data && Methylation Data   ####
rm(list = ls())
load('/home/bit/ashar/ExpressionSets/Verhark/IntegratedDataSurv.rda')

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


### Do some kind of variable selction to improve clustering

datag <- cbind(data.gene.combined, data.meth.combined)
library(clustvarsel)
cl <- clustvarsel(data = datag, G = 1:5, search = 'greedy', emModels1 = "E", direction = 'backward')


######################################
## Principal componets of Pathway Gene Data
pc <- prcomp(data.gene.combined)
pc.pred <- predict(pc,newdata = data.gene.combined)
plot(pc.pred[,1], pc.pred[,2], pch = 19, col =c)
## Principal Components of Methylation Data
pc <- prcomp(data.meth.combined)
pc.pred <- predict(pc,newdata = data.meth.combined)
plot(pc.pred[,1], pc.pred[,2], pch = 19)


