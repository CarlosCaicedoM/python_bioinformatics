# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 13:49:20 2021

@author: Carlos Caicedo-Montoya
"""
import pandas as pd
import numpy as np
from random import sample 
from numpy.random import randint

  



#ReadData

#Read_data
pan_matrix = pd.read_csv('pm.pfam.table2_file_Streptomyces_pangenomeR.RData.csv', 
                        low_memory=False)
pan_matrix.set_index('Genome', inplace=True)


pan_matrix = pan_matrix.T

# Transform it in a presence/absence matrix (1/0)
max_value= max(list(pan_matrix.max()))
presence_absence_matrix = pan_matrix.replace(list(range(1,max_value+1)), 1)
max(list(presence_absence_matrix.max()))

pan_matrix = presence_absence_matrix.T

cm = pan_matrix.iloc[randint(0, pan_matrix.shape[0], pan_matrix.shape[0]), :]
cm= cm.cumsum(axis=0, skipna=True)
ng = 121
nmat = np.zeros((pan_matrix.shape[0] - 1, 10))

nmat[:, 0] = 

cm2 = cm[(cm.iloc[1:ng,:] == 1) & (cm.iloc[0:(ng-1),:] == 0)]

cm.sum(axis = 0, [(cm.iloc[1:ng,:] == 1) & (cm.iloc[0:(ng-1),:] == 0)]) 

def heaps (pan_matrix, n_perm):
    ng = pan_matrix.shape[0]
    nmat = np.zeros((pan_matrix.shape[0] - 1, n_perm))
    for i in range(0, n_perm):
        cm = pan_matrix.iloc[randint(0, pan_matrix.shape[0], pan_matrix.shape[0]), :]
        cm= cm.cumsum(axis=0, skipna=True)
        nmat [:, i] = cm.sum(axis = 0, where=[(cm == 1)[2:ng,] & (cm == 0)[1:(ng-1),]]) 
        c

    
    
    
    
    
    
    
    
heaps <- function(pan.matrix, n.perm = 100){
  pan.matrix[which(pan.matrix > 0, arr.ind = T)] <- 1    #Convert in a presence/absence matrix (1/0)
  ng <- dim(pan.matrix)[1]      #Define the number of genomes
  nmat <- matrix(0, nrow = nrow(pan.matrix) - 1, ncol = n.perm) #Create a matrix of zeroes with 
                                                  #number of genomes -1 rows and n_perm columns
  for(i in 1:n.perm){
    cm <- apply(pan.matrix[sample(nrow(pan.matrix)),], 2, cumsum)  #get the cumulative sum 
                                                #of clusteer by randomly selected genomes
    nmat[,i] <- rowSums((cm == 1)[2:ng,] & (cm == 0)[1:(ng-1),])
    cat(i, "/", n.perm, "\r")
  }
  x <- rep((2:nrow(pan.matrix)), times = n.perm)
  y <- as.numeric(nmat)
  p0 <- c(mean(y[which(x == 2)] ), 1)
  fit <- optim(p0, objectFun, gr = NULL, x, y, method = "L-BFGS-B", lower = c(0, 0), upper = c(10000, 2))
  p.hat <- fit$par
  names(p.hat) <- c("Intercept", "alpha")
  return(p.hat)
}

objectFun <- function(p, x, y){
  y.hat <- p[1] * x^(-p[2])
  J <- sqrt(sum((y - y.hat)^2))/length(x)
  return(J)
}


