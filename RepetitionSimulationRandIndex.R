### This simulation calculates the RandIndices for the different methods
#### Copies from Simulation_Main from the one View Case
rm(list = ls())

RIground <- c(0)
CIground <- c(0)

for ( v in 1:20){
  
  
  ## Number of points
  N.test =  100
  N.train = 100
  
  N <- N.test
  ## Number of Clusters
  F = 2
  
  ## Distribution of the points within three clusters
  
  p.dist = c(0.5,0.5)
  
  ## Total Number of features D
  
  D1 = 20
  D2 = 20
  
  ## Total Percentage of irrelevant feature
  prob.noise.feature = 0.50
  
  ## Overlap between Cluster of molecular Data of the relevant features
  prob.overlap = 0.10

###### Get the Data #####################################

## Initialize the Training Data
source('simulatemultiDPMM.R')
simulatemultiDPMM()


############################# PARAMETERS for GIBB's SAMPLING ####
iter = 100
iter.burnin = 200
iter.thin  = 5
k = 2
K <-  as.integer(N)
F <- k
Time <- cbind(time,censoring)




########### Check Prediction Ground Truth
source('predictionGroundTruth.R')
predictionGroundTruth()


RIground[v] <- predRandIndex
CIground[v] <- predCIndex
}

boxplot(RIground, CIground)