## Analyzing the samples gotten from Gibbs' Sampler 
## Calculates Point Estimate for c
## Calculates summary statistics like RandIndex and C-Index
## Maybe Brier Scores

analyzemultiDPMM = function(){
  
  
  ############# TRAININIG DATA ###########################################
  #######################################################################
  ########################################################################
  #### This function calculates some important metrices for the TRAINING DATA Data
  #### Adjusted RandIndex
  #### C-Index
  #### Brier Scores for the Survival Curves
  #### Point Estimate of Clsuter Assignments
  ##### Brier Scores for feature Selection
  
    ########## ANLAYSING THE MCMC samples AND CALCULATING METRICES #######################################################
    
    Nps = as.integer(iter/ iter.thin)
    count <- Nps
    final.rand <- c(0)
    
    ############ The Matrices that will store the results #################################################
    final.rand <- c(0)
    cindex.final1 <- c(0)
    cindex.final2 <- c(0)   
    cindex.final <- c(0)  
    brier.final <- matrix(NA, nrow = Nps, ncol = F)
    
    ################ Begin Analysig the MCMC samples #######################################################
    
  
    
    for (h in 1:Nps){
      ### Adjusted Rand Indices
      #final.rand[h] <- adjustedRandIndex(c.list[[h]],as.factor(c.true))
      
      ### See C-Index (concordance index)
      surv.aft <- Surv(exp(time),censoring)
      
      ### Predict Time from the model
      source('predictlinear.R')
      tem.tim1 <- as.vector(unlist(predictlinear(c.list[[h]], Y1, est.regy1[[h]]$beta0, est.regy1[[h]]$betahat, est.regy1[[h]]$sigma2,D1)))
      tem.tim2 <- as.vector(unlist(predictlinear(c.list[[h]],Y2, est.regy2[[h]]$beta0, est.regy2[[h]]$betahat, est.regy2[[h]]$sigma2, D2)))
      
      cindex.final1[h] <-  survConcordance(surv.aft ~ exp(-tem.tim1))[[1]]
      cindex.final2[h] <-  survConcordance(surv.aft ~ exp(-tem.tim2))[[1]]
      
      source('multipredictlinear.R')
      tem.tim <- as.vector(unlist(multipredictlinear(c.list[[h]],  est.regy1[[h]], est.regy2[[h]] )))
      
      cindex.final[h] <-  survConcordance(surv.aft ~ exp(-tem.tim))[[1]]
    }
  
      
      
      ###############################################
      ###### Calculating POINT ESTIMATES ############
      ###############################################
    
    
       ##### Class Assignments ########################
        c.matrix <- matrix(NA, nrow = N, ncol = count)
        for ( i in 1:count){
        c.matrix[,i] <- c.list[[i]]
        }
        
      
    c.final <- apply(c.matrix,1,median)
    
    active <- as.numeric(rownames(table(c.final)))
    
    
    ############ Time Covariate Slopes FOR Relevant Clusters and Heatmap Plots ############
    list.betahat1 <- list(0)
    
    for ( i in 1:count){
      list.betahat1[[i]] <- (est.regy1[[i]]$betahat[active,] != 0) +0
    }
    
    Q <- length(active)
    matrix.betahat1 <- array(data = NA, dim =c(Q,count,D1))
    
    for ( z in 1:Q){
      for ( x  in 1:count){
        matrix.betahat1[z,x,] <- list.betahat1[[x]][z,]
      }
    }
    
    final.betahat1 <- apply(matrix.betahat1,c(1,3),mean)
    ### Probability of betahat of genes FOR ONE SIMULATION
    ##colnames(final.betahat) =  c(rep("relevant",rel.D),rep("irrelevant",irrel.D))
    heatmapdata1 <- as.data.frame(final.betahat1)
    #heatmap(as.matrix(heatmapdata1), col =cm.colors(180),main = "Posterior prob. \n for Selection for Data set 1 ", cexCol = 0.85, cexRow = 0.7)
     heatmap.2(t(as.matrix(heatmapdata1)),dendrogram="none", col =cm.colors(180), margins=c(6,10), main = "Posterior prob. \n for Selection for Data set 1 ", cexCol = 0.85, cexRow = 0.7, Rowv = FALSE)
    
    
    ########################## For the second data set ####################################################
    list.betahat2 <- list(0)
    
    for ( i in 1:count){
      list.betahat2[[i]] <- (est.regy2[[i]]$betahat[active,] != 0) +0
    }
    
    Q <- length(active)
    matrix.betahat2 <- array(data = NA, dim =c(Q,count,D2))
    
    for ( z in 1:Q){
      for ( x  in 1:count){
        matrix.betahat2[z,x,] <- list.betahat2[[x]][z,]
      }
    }
    
    final.betahat2 <- apply(matrix.betahat2,c(1,3),mean)
    heatmapdata2 <- as.data.frame(final.betahat2)
    
    #heatmap(as.matrix(heatmapdata2), col =cm.colors(180),main = "Posterior prob. \n for Selection for Data set 1 ", cexCol = 0.85, cexRow = 0.7)
    
    heatmap.2(t(as.matrix(heatmapdata2)),dendrogram="none", col =cm.colors(180), margins=c(6,10), main = "Posterior prob. \n for Selection for Data set 2 ", cexCol = 0.85, cexRow = 0.7, Rowv = FALSE)
    
    #final.rand <<- final.rand
    cindex.final1 <<- cindex.final1
    cindex.final2 <<- cindex.final2
    cindex.final <<- cindex.final
    c.final <<- c.final
    final.betahat1 <<- final.betahat1
    final.betahat2 <<- final.betahat2
    
    
  }
  