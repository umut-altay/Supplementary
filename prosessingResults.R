## This is the code script that can be used to extract the resulting values 
## of the analyses from the nested list "fullResults.RData", and then to 
## create the related figures and tables 

#Libraries 

library(reshape2)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(xtable)
library(viridis)
library(ggmap)

# 
# # # Result files
load("~/directory/fullResults.RData")  # fullResults.RData file for the Gaussian observation model
fittedGaussian = fitted
rm(fitted)

load("~/directory/fullResults.RData")  #fullResults.RData file for the Binomial observation model
fittedBinomial = fitted
rm(fitted)


nLikelihood = 1 #we ran Gaussian and binomial simulations seperately, so each 
                # "fullResults.RData" file is constructed on 1 likelihood
                # Use 2 If the input for the likelihood in Simulation.R is a 
                # matrix including values from both likelihoods.
nRange = 2  #number of spatial ranges
nScale = 2  # number of jittering schemes. we used 2 schemes (DHA and 4xDHS)
nSim = 50   #number of simulations per scenario
# 


#Extract each quantity from the fitted object, 
#remove the simulations which gave error and count the number of errors

Logscores1Binomial = list()
CRPS1Binomial = list()
coverage1Binomial = list()
Logscores2Binomial = list()
CRPS2Binomial = list()
coverage2Binomial = list()
errorBinomial = list()
factorBinomial = list()
rangeBinomial = list()
nuggetVarBinomial = list()
##########################
Logscores1Gaussian = list()
CRPS1Gaussian = list()
coverage1Gaussian = list()
Logscores2Gaussian = list()
CRPS2Gaussian = list()
coverage2Gaussian = list()
errorGaussian = list()
factorGaussian = list()
rangeGaussian = list()
nuggetVarGaussian = list()

for (i in 1:nLikelihood){
  temp1Logscores1Binomial = list()
  temp1CRPS1Binomial = list()
  temp1coverage1Binomial = list()
  temp1Logscores2Binomial = list()
  temp1CRPS2Binomial = list()
  temp1coverage2Binomial = list()
  error1Binomial = list()
  factor1Binomial = list()
  range1Binomial = list()
  nuggetVar1Binomial = list()
  #################################
  temp1Logscores1Gaussian = list()
  temp1CRPS1Gaussian = list()
  temp1coverage1Gaussian = list()
  temp1Logscores2Gaussian = list()
  temp1CRPS2Gaussian = list()
  temp1coverage2Gaussian = list()
  error1Gaussian = list()
  factor1Gaussian = list()
  range1Gaussian = list()
  nuggetVar1Gaussian = list()
  
  for (j in 1:nRange){
    temp2Logscores1Binomial = list()
    temp2CRPS1Binomial = list()
    temp2coverage1Binomial = list()
    temp2Logscores2Binomial = list()
    temp2CRPS2Binomial = list()
    temp2coverage2Binomial = list()
    error2Binomial = list()
    factor2Binomial = list()
    range2Binomial = list()
    nuggetVar2Binomial = list()
    #################################
    temp2Logscores1Gaussian = list()
    temp2CRPS1Gaussian = list()
    temp2coverage1Gaussian = list()
    temp2Logscores2Gaussian = list()
    temp2CRPS2Gaussian = list()
    temp2coverage2Gaussian = list()
    error2Gaussian = list()
    factor2Gaussian = list()
    range2Gaussian = list()
    nuggetVar2Gaussian = list()
    
    for (h in 1:nScale){
      temp3Logscores1Binomial = list()
      temp3CRPS1Binomial = list()
      temp3coverage1Binomial = list()
      temp3Logscores2Binomial = list()
      temp3CRPS2Binomial = list()
      temp3coverage2Binomial = list()
      factor3Binomial = list()
      range3Binomial = list()
      nuggetVar3Binomial = list()
      error3Binomial = 0
      ####################################
      temp3Logscores1Gaussian = list()
      temp3CRPS1Gaussian = list()
      temp3coverage1Gaussian = list()
      temp3Logscores2Gaussian = list()
      temp3CRPS2Gaussian = list()
      temp3coverage2Gaussian = list()
      factor3Gaussian = list()
      range3Gaussian = list()
      nuggetVar3Gaussian = list()
      error3Gaussian = 0
      
      for (l in 1:nSim){
        
        if(class(fittedBinomial[[1]][[1]][[i]][[j]][[h]][[l]])  == "try-error"){
          error3Binomial = error3Binomial + 1
          
        } else {
          
          temp3Logscores1Binomial[[l]] = try(fittedBinomial[[1]][[1]][[i]][[j]][[h]][[l]][["Logscores1"]])
          temp3CRPS1Binomial[[l]] =  try(fittedBinomial[[1]][[1]][[i]][[j]][[h]][[l]][["CRPSscores1"]])
          temp3coverage1Binomial[[l]] =  try(fittedBinomial[[1]][[1]][[i]][[j]][[h]][[l]][["coverage1"]])
          temp3Logscores2Binomial[[l]] =  try(fittedBinomial[[1]][[1]][[i]][[j]][[h]][[l]][["Logscores2"]])
          temp3CRPS2Binomial[[l]] =  try(fittedBinomial[[1]][[1]][[i]][[j]][[h]][[l]][["CRPSscores2"]])
          temp3coverage2Binomial[[l]] =  try(fittedBinomial[[1]][[1]][[i]][[j]][[h]][[l]][["coverage2"]])
        }
        
        
        if(class(fittedBinomial[[1]][[1]][[i]][[j]][[h]][[l]])  != "try-error"){
          if (h == 1) {factor3Binomial[[l]] = "1"                       # DHS jittering
          } else if (h == 2){factor3Binomial[[l]] = "4"                 # 4xDHS jittering
          } else if (h == 3){factor3Binomial[[l]] = ""                  # set this if you have more jittering schemes
          } else if (h == 4) {factor3Binomial[[l]] = ""                 # set this if you have more jittering schemes
          }
          
          if (j == 1) {range3Binomial[[l]] = "Range = 160 km"
          } else if (j == 2){range3Binomial[[l]] = "Range = 340 km"
          } else if (j == 3){range3Binomial[[l]] = ""                   # set this if you have more spatial ranges
          } 
          
          if (i == 1) {nuggetVar3Binomial[[l]] = ""                    # set this if you have nugget in binomial
          } else if (i == 2){nuggetVar3Binomial[[l]] = ""
          } else if (i == 3){nuggetVar3Binomial[[l]] = ""
          }
        }
        ########################################################################
        if(class(fittedGaussian[[1]][[1]][[i]][[j]][[h]][[l]])  == "try-error"){
          error3Gaussian = error3Gaussian + 1
          
        } else {
          
          temp3Logscores1Gaussian[[l]] = try(fittedGaussian[[1]][[1]][[i]][[j]][[h]][[l]][["Logscores1"]])
          temp3CRPS1Gaussian[[l]] =  try(fittedGaussian[[1]][[1]][[i]][[j]][[h]][[l]][["CRPSscores1"]])
          temp3coverage1Gaussian[[l]] =  try(fittedGaussian[[1]][[1]][[i]][[j]][[h]][[l]][["coverage1"]])
          temp3Logscores2Gaussian[[l]] =  try(fittedGaussian[[1]][[1]][[i]][[j]][[h]][[l]][["Logscores2"]])
          temp3CRPS2Gaussian[[l]] =  try(fittedGaussian[[1]][[1]][[i]][[j]][[h]][[l]][["CRPSscores2"]])
          temp3coverage2Gaussian[[l]] =  try(fittedGaussian[[1]][[1]][[i]][[j]][[h]][[l]][["coverage2"]])
        }
        
        
        if(class(fittedGaussian[[1]][[1]][[i]][[j]][[h]][[l]])  != "try-error"){
          if (h == 1) {factor3Gaussian[[l]] = "1"
          } else if (h == 2){factor3Gaussian[[l]] = "4"
          } else if (h == 3){factor3Gaussian[[l]] = ""                # set this if you have more jittering schemes
          } else if (h == 4) {factor3Gaussian[[l]] = ""               # set this if you have more jittering schemes
          }
          
          if (j == 1) {range3Gaussian[[l]] = "Range = 160 km"       
          } else if (j == 2){range3Gaussian[[l]] = "Range = 340 km"
          } else if (j == 3){range3Gaussian[[l]] = ""                 # set this if you have more spatial ranges
          } 
          
          if (i == 1) {nuggetVar3Gaussian[[l]] = "nuggetVar = 0.1"    # set this if you have more nugget variances
          } else if (i == 2){nuggetVar3Gaussian[[l]] = ""
          } else if (i == 3){nuggetVar3Gaussian[[l]] = ""
          }
        }
      ##########################################################################
      }
      temp2Logscores1Binomial[[h]] = unlist(temp3Logscores1Binomial)[!is.null(unlist(temp3Logscores1Binomial))]
      temp2CRPS1Binomial[[h]] = unlist(temp3CRPS1Binomial)[!is.null(unlist(temp3CRPS1Binomial))]
      temp2coverage1Binomial[[h]] = unlist(temp3coverage1Binomial)[!is.null(unlist(temp3coverage1Binomial))]
      temp2Logscores2Binomial[[h]] = unlist(temp3Logscores2Binomial)[!is.null(unlist(temp3Logscores2Binomial))]
      temp2CRPS2Binomial[[h]] = unlist(temp3CRPS2Binomial)[!is.null(unlist(temp3CRPS2Binomial))]
      temp2coverage2Binomial[[h]] = unlist(temp3coverage2Binomial)[!is.null(unlist(temp3coverage2Binomial))]
      error2Binomial[[h]] = error3Binomial
      factor2Binomial[[h]] = unlist(factor3Binomial)[!is.null(unlist(factor3Binomial))]
      range2Binomial[[h]] = unlist(range3Binomial)[!is.null(unlist(range3Binomial))]
      nuggetVar2Binomial[[h]] = unlist(nuggetVar3Binomial)[!is.null(unlist(nuggetVar3Binomial))]
      ##########################################################################
      temp2Logscores1Gaussian[[h]] = unlist(temp3Logscores1Gaussian)[!is.null(unlist(temp3Logscores1Gaussian))]
      temp2CRPS1Gaussian[[h]] = unlist(temp3CRPS1Gaussian)[!is.null(unlist(temp3CRPS1Gaussian))]
      temp2coverage1Gaussian[[h]] = unlist(temp3coverage1Gaussian)[!is.null(unlist(temp3coverage1Gaussian))]
      temp2Logscores2Gaussian[[h]] = unlist(temp3Logscores2Gaussian)[!is.null(unlist(temp3Logscores2Gaussian))]
      temp2CRPS2Gaussian[[h]] = unlist(temp3CRPS2Gaussian)[!is.null(unlist(temp3CRPS2Gaussian))]
      temp2coverage2Gaussian[[h]] = unlist(temp3coverage2Gaussian)[!is.null(unlist(temp3coverage2Gaussian))]
      error2Gaussian[[h]] = error3Gaussian
      factor2Gaussian[[h]] = unlist(factor3Gaussian)[!is.null(unlist(factor3Gaussian))]
      range2Gaussian[[h]] = unlist(range3Gaussian)[!is.null(unlist(range3Gaussian))]
      nuggetVar2Gaussian[[h]] = unlist(nuggetVar3Gaussian)[!is.null(unlist(nuggetVar3Gaussian))]
    }
    temp1Logscores1Binomial[[j]] = temp2Logscores1Binomial
    temp1CRPS1Binomial[[j]] = temp2CRPS1Binomial
    temp1coverage1Binomial[[j]] = temp2coverage1Binomial
    temp1Logscores2Binomial[[j]] = temp2Logscores2Binomial
    temp1CRPS2Binomial[[j]] = temp2CRPS2Binomial
    temp1coverage2Binomial[[j]] = temp2coverage2Binomial
    error1Binomial[[j]] = error2Binomial
    factor1Binomial[[j]] = factor2Binomial
    range1Binomial[[j]] = range2Binomial
    nuggetVar1Binomial[[j]] = nuggetVar2Binomial
    ############################################################################
    temp1Logscores1Gaussian[[j]] = temp2Logscores1Gaussian
    temp1CRPS1Gaussian[[j]] = temp2CRPS1Gaussian
    temp1coverage1Gaussian[[j]] = temp2coverage1Gaussian
    temp1Logscores2Gaussian[[j]] = temp2Logscores2Gaussian
    temp1CRPS2Gaussian[[j]] = temp2CRPS2Gaussian
    temp1coverage2Gaussian[[j]] = temp2coverage2Gaussian
    error1Gaussian[[j]] = error2Gaussian
    factor1Gaussian[[j]] = factor2Gaussian
    range1Gaussian[[j]] = range2Gaussian
    nuggetVar1Gaussian[[j]] = nuggetVar2Gaussian
  }
  Logscores1Binomial[[i]] = temp1Logscores1Binomial
  CRPS1Binomial[[i]] = temp1CRPS1Binomial
  coverage1Binomial[[i]] = temp1coverage1Binomial
  Logscores2Binomial[[i]] = temp1Logscores2Binomial
  CRPS2Binomial[[i]] = temp1CRPS2Binomial
  coverage2Binomial[[i]] = temp1coverage2Binomial
  errorBinomial[[i]] = error1Binomial
  factorBinomial[[i]] = factor1Binomial
  rangeBinomial[[i]] = range1Binomial
  nuggetVarBinomial[[i]] = nuggetVar1Binomial
  ##############################################################################
  Logscores1Gaussian[[i]] = temp1Logscores1Gaussian
  CRPS1Gaussian[[i]] = temp1CRPS1Gaussian
  coverage1Gaussian[[i]] = temp1coverage1Gaussian
  Logscores2Gaussian[[i]] = temp1Logscores2Gaussian
  CRPS2Gaussian[[i]] = temp1CRPS2Gaussian
  coverage2Gaussian[[i]] = temp1coverage2Gaussian
  errorGaussian[[i]] = error1Gaussian
  factorGaussian[[i]] = factor1Gaussian
  rangeGaussian[[i]] = range1Gaussian
  nuggetVarGaussian[[i]] = nuggetVar1Gaussian
}


#PLOTS

#Create the data frames first:     # use commented out lines when you have more input values

graphMainBinomial = data.frame(CRPS1Binomial = c(unlist(CRPS1Binomial[[1]][[1]][[1]]), unlist(CRPS1Binomial[[1]][[1]][[2]]), #unlist(CRPS1Binomial[[1]][[1]][[3]]), #unlist(CRPS1Binomial[[1]][[1]][[4]]),
                                 unlist(CRPS1Binomial[[1]][[2]][[1]]), unlist(CRPS1Binomial[[1]][[2]][[2]])), #unlist(CRPS1Binomial[[1]][[2]][[3]])),#,  unlist(CRPS1Binomial[[1]][[2]][[4]])),
                       # unlist(CRPS1Binomial[[1]][[3]][[1]]), unlist(CRPS1Binomial[[1]][[3]][[2]]), unlist(CRPS1Binomial[[1]][[3]][[3]]), 
                       # unlist(CRPS1Binomial[[2]][[1]][[1]]), unlist(CRPS1Binomial[[2]][[1]][[2]]), unlist(CRPS1Binomial[[2]][[1]][[3]]),
                       # unlist(CRPS1Binomial[[2]][[2]][[1]]), unlist(CRPS1Binomial[[2]][[2]][[2]]), unlist(CRPS1Binomial[[2]][[2]][[3]]),
                       # unlist(CRPS1Binomial[[2]][[3]][[1]]), unlist(CRPS1Binomial[[2]][[3]][[2]]), unlist(CRPS1Binomial[[2]][[3]][[3]]), 
                       # unlist(CRPS1Binomial[[3]][[1]][[1]]), unlist(CRPS1Binomial[[3]][[1]][[2]]), unlist(CRPS1Binomial[[3]][[1]][[3]]),
                       # unlist(CRPS1Binomial[[3]][[2]][[1]]), unlist(CRPS1Binomial[[3]][[2]][[2]]), unlist(CRPS1Binomial[[3]][[2]][[3]]),
                       # unlist(CRPS1Binomial[[3]][[3]][[1]]), unlist(CRPS1Binomial[[3]][[3]][[2]]), unlist(CRPS1Binomial[[3]][[3]][[3]])), 
                       CRPS2Binomial = c(unlist(CRPS2Binomial[[1]][[1]][[1]]), unlist(CRPS2Binomial[[1]][[1]][[2]]), #unlist(CRPS2Binomial[[1]][[1]][[3]]),# unlist(CRPS2Binomial[[1]][[1]][[4]]),
                                 unlist(CRPS2Binomial[[1]][[2]][[1]]), unlist(CRPS2Binomial[[1]][[2]][[2]])),# #unlist(CRPS2Binomial[[1]][[2]][[3]])), #unlist(CRPS2Binomial[[1]][[2]][[4]])),
                       # unlist(CRPS2Binomial[[1]][[3]][[1]]), unlist(CRPS2Binomial[[1]][[3]][[2]]), unlist(CRPS2Binomial[[1]][[3]][[3]]), 
                       # unlist(CRPS2Binomial[[2]][[1]][[1]]), unlist(CRPS2Binomial[[2]][[1]][[2]]), unlist(CRPS2Binomial[[2]][[1]][[3]]),
                       # unlist(CRPS2Binomial[[2]][[2]][[1]]), unlist(CRPS2Binomial[[2]][[2]][[2]]), unlist(CRPS2Binomial[[2]][[2]][[3]]),
                       # unlist(CRPS2Binomial[[2]][[3]][[1]]), unlist(CRPS2Binomial[[2]][[3]][[2]]), unlist(CRPS2Binomial[[2]][[3]][[3]]), 
                       # unlist(CRPS2Binomial[[3]][[1]][[1]]), unlist(CRPS2Binomial[[3]][[1]][[2]]), unlist(CRPS2Binomial[[3]][[1]][[3]]),
                       # unlist(CRPS2Binomial[[3]][[2]][[1]]), unlist(CRPS2Binomial[[3]][[2]][[2]]), unlist(CRPS2Binomial[[3]][[2]][[3]]),
                       # unlist(CRPS2Binomial[[3]][[3]][[1]]), unlist(CRPS2Binomial[[3]][[3]][[2]]), unlist(CRPS2Binomial[[3]][[3]][[3]])),
                       Logscores1Binomial = c(unlist(Logscores1Binomial[[1]][[1]][[1]]), unlist(Logscores1Binomial[[1]][[1]][[2]]), #unlist(Logscores1Binomial[[1]][[1]][[3]]), #unlist(Logscores1Binomial[[1]][[1]][[4]]),
                                      unlist(Logscores1Binomial[[1]][[2]][[1]]), unlist(Logscores1Binomial[[1]][[2]][[2]])),# unlist(Logscores1Binomial[[1]][[2]][[3]])),# unlist(Logscores1Binomial[[1]][[2]][[4]])),
                       # unlist(Logscores1Binomial[[1]][[3]][[1]]), unlist(Logscores1Binomial[[1]][[3]][[2]]), unlist(Logscores1Binomial[[1]][[3]][[3]]), 
                       # unlist(Logscores1Binomial[[2]][[1]][[1]]), unlist(Logscores1Binomial[[2]][[1]][[2]]), unlist(Logscores1Binomial[[2]][[1]][[3]]),
                       # unlist(Logscores1Binomial[[2]][[2]][[1]]), unlist(Logscores1Binomial[[2]][[2]][[2]]), unlist(Logscores1Binomial[[2]][[2]][[3]]),
                       # unlist(Logscores1Binomial[[2]][[3]][[1]]), unlist(Logscores1Binomial[[2]][[3]][[2]]), unlist(Logscores1Binomial[[2]][[3]][[3]]), 
                       # unlist(Logscores1Binomial[[3]][[1]][[1]]), unlist(Logscores1Binomial[[3]][[1]][[2]]), unlist(Logscores1Binomial[[3]][[1]][[3]]),
                       # unlist(Logscores1Binomial[[3]][[2]][[1]]), unlist(Logscores1Binomial[[3]][[2]][[2]]), unlist(Logscores1Binomial[[3]][[2]][[3]]),
                       # unlist(Logscores1Binomial[[3]][[3]][[1]]), unlist(Logscores1Binomial[[3]][[3]][[2]]), unlist(Logscores1Binomial[[3]][[3]][[3]])),
                       Logscores2Binomial = c(unlist(Logscores2Binomial[[1]][[1]][[1]]), unlist(Logscores2Binomial[[1]][[1]][[2]]), #unlist(Logscores2Binomial[[1]][[1]][[3]]), #unlist(Logscores2Binomial[[1]][[1]][[4]]),
                                      unlist(Logscores2Binomial[[1]][[2]][[1]]), unlist(Logscores2Binomial[[1]][[2]][[2]])),# unlist(Logscores2Binomial[[1]][[2]][[3]])), #unlist(Logscores2Binomial[[1]][[2]][[4]])),
                       # unlist(Logscores2Binomial[[1]][[3]][[1]]), unlist(Logscores2Binomial[[1]][[3]][[2]]), unlist(Logscores2Binomial[[1]][[3]][[3]]), 
                       # unlist(Logscores2Binomial[[2]][[1]][[1]]), unlist(Logscores2Binomial[[2]][[1]][[2]]), unlist(Logscores2Binomial[[2]][[1]][[3]]),
                       # unlist(Logscores2Binomial[[2]][[2]][[1]]), unlist(Logscores2Binomial[[2]][[2]][[2]]), unlist(Logscores2Binomial[[2]][[2]][[3]]),
                       # unlist(Logscores2Binomial[[2]][[3]][[1]]), unlist(Logscores2Binomial[[2]][[3]][[2]]), unlist(Logscores2Binomial[[2]][[3]][[3]]), 
                       # unlist(Logscores2Binomial[[3]][[1]][[1]]), unlist(Logscores2Binomial[[3]][[1]][[2]]), unlist(Logscores2Binomial[[3]][[1]][[3]]),
                       # unlist(Logscores2Binomial[[3]][[2]][[1]]), unlist(Logscores2Binomial[[3]][[2]][[2]]), unlist(Logscores2Binomial[[3]][[2]][[3]]),
                       # unlist(Logscores2Binomial[[3]][[3]][[1]]), unlist(Logscores2Binomial[[3]][[3]][[2]]), unlist(Logscores2Binomial[[3]][[3]][[3]])),
                       coverage1Binomial = c(unlist(coverage1Binomial[[1]][[1]][[1]]), unlist(coverage1Binomial[[1]][[1]][[2]]), #unlist(coverage1Binomial[[1]][[1]][[3]]), #unlist(coverage1Binomial[[1]][[1]][[4]]),
                                     unlist(coverage1Binomial[[1]][[2]][[1]]), unlist(coverage1Binomial[[1]][[2]][[2]])),# unlist(coverage1Binomial[[1]][[2]][[3]])), #unlist(coverage1Binomial[[1]][[2]][[4]])),
                       # unlist(coverage1Binomial[[1]][[3]][[1]]), unlist(coverage1Binomial[[1]][[3]][[2]]), unlist(coverage1Binomial[[1]][[3]][[3]]), 
                       # unlist(coverage1Binomial[[2]][[1]][[1]]), unlist(coverage1Binomial[[2]][[1]][[2]]), unlist(coverage1Binomial[[2]][[1]][[3]]),
                       # unlist(coverage1Binomial[[2]][[2]][[1]]), unlist(coverage1Binomial[[2]][[2]][[2]]), unlist(coverage1Binomial[[2]][[2]][[3]]),
                       # unlist(coverage1Binomial[[2]][[3]][[1]]), unlist(coverage1Binomial[[2]][[3]][[2]]), unlist(coverage1Binomial[[2]][[3]][[3]]), 
                       # unlist(coverage1Binomial[[3]][[1]][[1]]), unlist(coverage1Binomial[[3]][[1]][[2]]), unlist(coverage1Binomial[[3]][[1]][[3]]),
                       # unlist(coverage1Binomial[[3]][[2]][[1]]), unlist(coverage1Binomial[[3]][[2]][[2]]), unlist(coverage1Binomial[[3]][[2]][[3]]),
                       # unlist(coverage1Binomial[[3]][[3]][[1]]), unlist(coverage1Binomial[[3]][[3]][[2]]), unlist(coverage1Binomial[[3]][[3]][[3]])),
                       coverage2Binomial = c(unlist(coverage2Binomial[[1]][[1]][[1]]), unlist(coverage2Binomial[[1]][[1]][[2]]), #unlist(coverage2Binomial[[1]][[1]][[3]]), #unlist(coverage2Binomial[[1]][[1]][[4]]),
                                     unlist(coverage2Binomial[[1]][[2]][[1]]), unlist(coverage2Binomial[[1]][[2]][[2]])),# unlist(coverage2Binomial[[1]][[2]][[3]])), #unlist(coverage2Binomial[[1]][[2]][[4]])),
                       # unlist(coverage2Binomial[[1]][[3]][[1]]), unlist(coverage2Binomial[[1]][[3]][[2]]), unlist(coverage2Binomial[[1]][[3]][[3]]), 
                       # unlist(coverage2Binomial[[2]][[1]][[1]]), unlist(coverage2Binomial[[2]][[1]][[2]]), unlist(coverage2Binomial[[2]][[1]][[3]]),
                       # unlist(coverage2Binomial[[2]][[2]][[1]]), unlist(coverage2Binomial[[2]][[2]][[2]]), unlist(coverage2Binomial[[2]][[2]][[3]]),
                       # unlist(coverage2Binomial[[2]][[3]][[1]]), unlist(coverage2Binomial[[2]][[3]][[2]]), unlist(coverage2Binomial[[2]][[3]][[3]]), 
                       # unlist(coverage2Binomial[[3]][[1]][[1]]), unlist(coverage2Binomial[[3]][[1]][[2]]), unlist(coverage2Binomial[[3]][[1]][[3]]),
                       # unlist(coverage2Binomial[[3]][[2]][[1]]), unlist(coverage2Binomial[[3]][[2]][[2]]), unlist(coverage2Binomial[[3]][[2]][[3]]),
                       # unlist(coverage2Binomial[[3]][[3]][[1]]), unlist(coverage2Binomial[[3]][[3]][[2]]), unlist(coverage2Binomial[[3]][[3]][[3]])),
                       factorBinomial = c(unlist(factorBinomial[[1]][[1]][[1]]), unlist(factorBinomial[[1]][[1]][[2]]),# unlist(factorBinomial[[1]][[1]][[3]]), #unlist(factorBinomial[[1]][[1]][[4]]),
                                  unlist(factorBinomial[[1]][[2]][[1]]), unlist(factorBinomial[[1]][[2]][[2]])),# unlist(factorBinomial[[1]][[2]][[3]])), #unlist(factorBinomial[[1]][[2]][[4]])),
                       # unlist(factorBinomial[[1]][[3]][[1]]), unlist(factorBinomial[[1]][[3]][[2]]), unlist(factorBinomial[[1]][[3]][[3]]), 
                       # unlist(factorBinomial[[2]][[1]][[1]]), unlist(factorBinomial[[2]][[1]][[2]]), unlist(factorBinomial[[2]][[1]][[3]]),
                       # unlist(factorBinomial[[2]][[2]][[1]]), unlist(factorBinomial[[2]][[2]][[2]]), unlist(factorBinomial[[2]][[2]][[3]]),
                       # unlist(factorBinomial[[2]][[3]][[1]]), unlist(factorBinomial[[2]][[3]][[2]]), unlist(factorBinomial[[2]][[3]][[3]]), 
                       # unlist(factorBinomial[[3]][[1]][[1]]), unlist(factorBinomial[[3]][[1]][[2]]), unlist(factorBinomial[[3]][[1]][[3]]),
                       # unlist(factorBinomial[[3]][[2]][[1]]), unlist(factorBinomial[[3]][[2]][[2]]), unlist(factorBinomial[[3]][[2]][[3]]),
                       # unlist(factorBinomial[[3]][[3]][[1]]), unlist(factorBinomial[[3]][[3]][[2]]), unlist(factorBinomial[[3]][[3]][[3]])),
                       rangeBinomial = c(unlist(rangeBinomial[[1]][[1]][[1]]), unlist(rangeBinomial[[1]][[1]][[2]]),# unlist(rangeBinomial[[1]][[1]][[3]]), #unlist(rangeBinomial[[1]][[1]][[4]]),
                                 unlist(rangeBinomial[[1]][[2]][[1]]), unlist(rangeBinomial[[1]][[2]][[2]])),# unlist(rangeBinomial[[1]][[2]][[3]])), #unlist(rangeBinomial[[1]][[2]][[4]])),
                       # unlist(rangeBinomial[[1]][[3]][[1]]), unlist(rangeBinomial[[1]][[3]][[2]]), unlist(rangeBinomial[[1]][[3]][[3]]), 
                       # unlist(rangeBinomial[[2]][[1]][[1]]), unlist(rangeBinomial[[2]][[1]][[2]]), unlist(rangeBinomial[[2]][[1]][[3]]),
                       # unlist(rangeBinomial[[2]][[2]][[1]]), unlist(rangeBinomial[[2]][[2]][[2]]), unlist(rangeBinomial[[2]][[2]][[3]]),
                       # unlist(rangeBinomial[[2]][[3]][[1]]), unlist(rangeBinomial[[2]][[3]][[2]]), unlist(rangeBinomial[[2]][[3]][[3]]), 
                       # unlist(rangeBinomial[[3]][[1]][[1]]), unlist(rangeBinomial[[3]][[1]][[2]]), unlist(rangeBinomial[[3]][[1]][[3]]),
                       # unlist(rangeBinomial[[3]][[2]][[1]]), unlist(rangeBinomial[[3]][[2]][[2]]), unlist(rangeBinomial[[3]][[2]][[3]]),  
                       # unlist(rangeBinomial[[3]][[3]][[1]]), unlist(rangeBinomial[[3]][[3]][[2]]), unlist(rangeBinomial[[3]][[3]][[3]])),
                       nuggetVarBinomial = c(unlist(nuggetVarBinomial[[1]][[1]][[1]]), unlist(nuggetVarBinomial[[1]][[1]][[2]]),# unlist(nuggetVarBinomial[[1]][[1]][[3]]), #unlist(nuggetVarBinomial[[1]][[1]][[4]]),
                                     unlist(nuggetVarBinomial[[1]][[2]][[1]]), unlist(nuggetVarBinomial[[1]][[2]][[2]])),# unlist(nuggetVarBinomial[[1]][[2]][[3]])), #unlist(nuggetVarBinomial[[1]][[2]][[4]])),
                       # unlist(nuggetVarBinomial[[1]][[3]][[1]]), unlist(nuggetVarBinomial[[1]][[3]][[2]]), unlist(nuggetVarBinomial[[1]][[3]][[3]]), 
                       # unlist(nuggetVarBinomial[[2]][[1]][[1]]), unlist(nuggetVarBinomial[[2]][[1]][[2]]), unlist(nuggetVarBinomial[[2]][[1]][[3]]),
                       # unlist(nuggetVarBinomial[[2]][[2]][[1]]), unlist(nuggetVarBinomial[[2]][[2]][[2]]), unlist(nuggetVarBinomial[[2]][[2]][[3]]),
                       # unlist(nuggetVarBinomial[[2]][[3]][[1]]), unlist(nuggetVarBinomial[[2]][[3]][[2]]), unlist(nuggetVarBinomial[[2]][[3]][[3]]), 
                       # unlist(nuggetVarBinomial[[3]][[1]][[1]]), unlist(nuggetVarBinomial[[3]][[1]][[2]]), unlist(nuggetVarBinomial[[3]][[1]][[3]]),
                       # unlist(nuggetVarBinomial[[3]][[2]][[1]]), unlist(nuggetVarBinomial[[3]][[2]][[2]]), unlist(nuggetVarBinomial[[3]][[2]][[3]]),
                       # unlist(nuggetVarBinomial[[3]][[3]][[1]]), unlist(nuggetVarBinomial[[3]][[3]][[2]]), unlist(nuggetVarBinomial[[3]][[3]][[3]])),
                       errorBinomial = c(rep(unlist(errorBinomial[[1]][[1]][[1]]), length(unlist(CRPS1Binomial[[1]][[1]][[1]]))), rep(unlist(errorBinomial[[1]][[1]][[2]]), length(unlist(CRPS1Binomial[[1]][[1]][[2]]))), #rep(unlist(errorBinomial[[1]][[1]][[3]]), length(unlist(CRPS1Binomial[[1]][[1]][[3]]))),
                                 #rep(unlist(errorBinomial[[1]][[1]][[4]]), length(unlist(CRPS1Binomial[[1]][[1]][[4]]))),
                                 rep(unlist(errorBinomial[[1]][[2]][[1]]), length(unlist(CRPS1Binomial[[1]][[2]][[1]]))), rep(unlist(errorBinomial[[1]][[2]][[2]]), length(unlist(CRPS1Binomial[[1]][[2]][[2]])))))#, rep(unlist(errorBinomial[[1]][[2]][[3]]), length(unlist(CRPS1Binomial[[1]][[2]][[3]])))))#,
#rep(unlist(errorBinomial[[1]][[2]][[4]]), length(unlist(CRPS1Binomial[[1]][[2]][[4]])))))
# rep(unlist(errorBinomial[[1]][[3]][[1]]), length(unlist(CRPS1Binomial[[1]][[3]][[1]]))), rep(unlist(errorBinomial[[1]][[3]][[2]]), length(unlist(CRPS1Binomial[[1]][[3]][[2]]))), rep(unlist(errorBinomial[[1]][[3]][[3]]), length(unlist(CRPS1Binomial[[1]][[3]][[3]]))), 
# rep(unlist(errorBinomial[[2]][[1]][[1]]), length(unlist(CRPS1Binomial[[2]][[1]][[1]]))), rep(unlist(errorBinomial[[2]][[1]][[2]]), length(unlist(CRPS1Binomial[[2]][[1]][[2]]))), rep(unlist(errorBinomial[[1]][[3]][[3]]), length(unlist(CRPS1Binomial[[2]][[1]][[3]]))),
# rep(unlist(errorBinomial[[2]][[2]][[1]]), length(unlist(CRPS1Binomial[[2]][[2]][[1]]))), rep(unlist(errorBinomial[[2]][[2]][[2]]), length(unlist(CRPS1Binomial[[2]][[2]][[2]]))), rep(unlist(errorBinomial[[2]][[2]][[3]]), length(unlist(CRPS1Binomial[[2]][[2]][[3]]))),
# rep(unlist(errorBinomial[[2]][[3]][[1]]), length(unlist(CRPS1Binomial[[2]][[3]][[1]]))), rep(unlist(errorBinomial[[2]][[3]][[2]]), length(unlist(CRPS1Binomial[[2]][[3]][[2]]))), rep(unlist(errorBinomial[[2]][[2]][[3]]), length(unlist(CRPS1Binomial[[2]][[3]][[3]]))), 
# rep(unlist(errorBinomial[[3]][[1]][[1]]), length(unlist(CRPS1Binomial[[3]][[1]][[1]]))), rep(unlist(errorBinomial[[3]][[1]][[2]]), length(unlist(CRPS1Binomial[[3]][[1]][[2]]))), rep(unlist(errorBinomial[[3]][[1]][[3]]), length(unlist(CRPS1Binomial[[3]][[1]][[3]]))),
# rep(unlist(errorBinomial[[3]][[2]][[1]]), length(unlist(CRPS1Binomial[[3]][[2]][[1]]))), rep(unlist(errorBinomial[[3]][[2]][[2]]), length(unlist(CRPS1Binomial[[3]][[2]][[2]]))), rep(unlist(errorBinomial[[3]][[2]][[3]]), length(unlist(CRPS1Binomial[[3]][[2]][[3]]))),
# rep(unlist(errorBinomial[[3]][[3]][[1]]), length(unlist(CRPS1Binomial[[3]][[3]][[1]]))), rep(unlist(errorBinomial[[3]][[3]][[2]]), length(unlist(CRPS1Binomial[[3]][[3]][[2]]))), rep(unlist(errorBinomial[[3]][[3]][[3]]), length(unlist(CRPS1Binomial[[3]][[3]][[3]])))))
#


#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


graphMainGaussian = data.frame(CRPS1Gaussian = c(unlist(CRPS1Gaussian[[1]][[1]][[1]]), unlist(CRPS1Gaussian[[1]][[1]][[2]]), #unlist(CRPS1Gaussian[[1]][[1]][[3]]), #unlist(CRPS1Gaussian[[1]][[1]][[4]]),
                                                 unlist(CRPS1Gaussian[[1]][[2]][[1]]), unlist(CRPS1Gaussian[[1]][[2]][[2]])), #unlist(CRPS1Gaussian[[1]][[2]][[3]])),#,  unlist(CRPS1Gaussian[[1]][[2]][[4]])),
                               # unlist(CRPS1Gaussian[[1]][[3]][[1]]), unlist(CRPS1Gaussian[[1]][[3]][[2]]), unlist(CRPS1Gaussian[[1]][[3]][[3]]), 
                               # unlist(CRPS1Gaussian[[2]][[1]][[1]]), unlist(CRPS1Gaussian[[2]][[1]][[2]]), unlist(CRPS1Gaussian[[2]][[1]][[3]]),
                               # unlist(CRPS1Gaussian[[2]][[2]][[1]]), unlist(CRPS1Gaussian[[2]][[2]][[2]]), unlist(CRPS1Gaussian[[2]][[2]][[3]]),
                               # unlist(CRPS1Gaussian[[2]][[3]][[1]]), unlist(CRPS1Gaussian[[2]][[3]][[2]]), unlist(CRPS1Gaussian[[2]][[3]][[3]]), 
                               # unlist(CRPS1Gaussian[[3]][[1]][[1]]), unlist(CRPS1Gaussian[[3]][[1]][[2]]), unlist(CRPS1Gaussian[[3]][[1]][[3]]),
                               # unlist(CRPS1Gaussian[[3]][[2]][[1]]), unlist(CRPS1Gaussian[[3]][[2]][[2]]), unlist(CRPS1Gaussian[[3]][[2]][[3]]),
                               # unlist(CRPS1Gaussian[[3]][[3]][[1]]), unlist(CRPS1Gaussian[[3]][[3]][[2]]), unlist(CRPS1Gaussian[[3]][[3]][[3]])), 
                               CRPS2Gaussian = c(unlist(CRPS2Gaussian[[1]][[1]][[1]]), unlist(CRPS2Gaussian[[1]][[1]][[2]]), #unlist(CRPS2Gaussian[[1]][[1]][[3]]),# unlist(CRPS2Gaussian[[1]][[1]][[4]]),
                                                 unlist(CRPS2Gaussian[[1]][[2]][[1]]), unlist(CRPS2Gaussian[[1]][[2]][[2]])),# #unlist(CRPS2Gaussian[[1]][[2]][[3]])), #unlist(CRPS2Gaussian[[1]][[2]][[4]])),
                               # unlist(CRPS2Gaussian[[1]][[3]][[1]]), unlist(CRPS2Gaussian[[1]][[3]][[2]]), unlist(CRPS2Gaussian[[1]][[3]][[3]]), 
                               # unlist(CRPS2Gaussian[[2]][[1]][[1]]), unlist(CRPS2Gaussian[[2]][[1]][[2]]), unlist(CRPS2Gaussian[[2]][[1]][[3]]),
                               # unlist(CRPS2Gaussian[[2]][[2]][[1]]), unlist(CRPS2Gaussian[[2]][[2]][[2]]), unlist(CRPS2Gaussian[[2]][[2]][[3]]),
                               # unlist(CRPS2Gaussian[[2]][[3]][[1]]), unlist(CRPS2Gaussian[[2]][[3]][[2]]), unlist(CRPS2Gaussian[[2]][[3]][[3]]), 
                               # unlist(CRPS2Gaussian[[3]][[1]][[1]]), unlist(CRPS2Gaussian[[3]][[1]][[2]]), unlist(CRPS2Gaussian[[3]][[1]][[3]]),
                               # unlist(CRPS2Gaussian[[3]][[2]][[1]]), unlist(CRPS2Gaussian[[3]][[2]][[2]]), unlist(CRPS2Gaussian[[3]][[2]][[3]]),
                               # unlist(CRPS2Gaussian[[3]][[3]][[1]]), unlist(CRPS2Gaussian[[3]][[3]][[2]]), unlist(CRPS2Gaussian[[3]][[3]][[3]])),
                               Logscores1Gaussian = c(unlist(Logscores1Gaussian[[1]][[1]][[1]]), unlist(Logscores1Gaussian[[1]][[1]][[2]]), #unlist(Logscores1Gaussian[[1]][[1]][[3]]), #unlist(Logscores1Gaussian[[1]][[1]][[4]]),
                                                      unlist(Logscores1Gaussian[[1]][[2]][[1]]), unlist(Logscores1Gaussian[[1]][[2]][[2]])),# unlist(Logscores1Gaussian[[1]][[2]][[3]])),# unlist(Logscores1Gaussian[[1]][[2]][[4]])),
                               # unlist(Logscores1Gaussian[[1]][[3]][[1]]), unlist(Logscores1Gaussian[[1]][[3]][[2]]), unlist(Logscores1Gaussian[[1]][[3]][[3]]), 
                               # unlist(Logscores1Gaussian[[2]][[1]][[1]]), unlist(Logscores1Gaussian[[2]][[1]][[2]]), unlist(Logscores1Gaussian[[2]][[1]][[3]]),
                               # unlist(Logscores1Gaussian[[2]][[2]][[1]]), unlist(Logscores1Gaussian[[2]][[2]][[2]]), unlist(Logscores1Gaussian[[2]][[2]][[3]]),
                               # unlist(Logscores1Gaussian[[2]][[3]][[1]]), unlist(Logscores1Gaussian[[2]][[3]][[2]]), unlist(Logscores1Gaussian[[2]][[3]][[3]]), 
                               # unlist(Logscores1Gaussian[[3]][[1]][[1]]), unlist(Logscores1Gaussian[[3]][[1]][[2]]), unlist(Logscores1Gaussian[[3]][[1]][[3]]),
                               # unlist(Logscores1Gaussian[[3]][[2]][[1]]), unlist(Logscores1Gaussian[[3]][[2]][[2]]), unlist(Logscores1Gaussian[[3]][[2]][[3]]),
                               # unlist(Logscores1Gaussian[[3]][[3]][[1]]), unlist(Logscores1Gaussian[[3]][[3]][[2]]), unlist(Logscores1Gaussian[[3]][[3]][[3]])),
                               Logscores2Gaussian = c(unlist(Logscores2Gaussian[[1]][[1]][[1]]), unlist(Logscores2Gaussian[[1]][[1]][[2]]), #unlist(Logscores2Gaussian[[1]][[1]][[3]]), #unlist(Logscores2Binomial[[1]][[1]][[4]]),
                                                      unlist(Logscores2Gaussian[[1]][[2]][[1]]), unlist(Logscores2Gaussian[[1]][[2]][[2]])),# unlist(Logscores2Gaussian[[1]][[2]][[3]])), #unlist(Logscores2Binomial[[1]][[2]][[4]])),
                               # unlist(Logscores2Gaussian[[1]][[3]][[1]]), unlist(Logscores2Gaussian[[1]][[3]][[2]]), unlist(Logscores2Gaussian[[1]][[3]][[3]]), 
                               # unlist(Logscores2Gaussian[[2]][[1]][[1]]), unlist(Logscores2Gaussian[[2]][[1]][[2]]), unlist(Logscores2Gaussian[[2]][[1]][[3]]),
                               # unlist(Logscores2Gaussian[[2]][[2]][[1]]), unlist(Logscores2Gaussian[[2]][[2]][[2]]), unlist(Logscores2Gaussian[[2]][[2]][[3]]),
                               # unlist(Logscores2Gaussian[[2]][[3]][[1]]), unlist(Logscores2Gaussian[[2]][[3]][[2]]), unlist(Logscores2Gaussian[[2]][[3]][[3]]), 
                               # unlist(Logscores2Gaussian[[3]][[1]][[1]]), unlist(Logscores2Gaussian[[3]][[1]][[2]]), unlist(Logscores2Gaussian[[3]][[1]][[3]]),
                               # unlist(Logscores2Gaussian[[3]][[2]][[1]]), unlist(Logscores2Gaussian[[3]][[2]][[2]]), unlist(Logscores2Gaussian[[3]][[2]][[3]]),
                               # unlist(Logscores2Gaussian[[3]][[3]][[1]]), unlist(Logscores2Gaussian[[3]][[3]][[2]]), unlist(Logscores2Gaussian[[3]][[3]][[3]])),
                               coverage1Gaussian = c(unlist(coverage1Gaussian[[1]][[1]][[1]]), unlist(coverage1Gaussian[[1]][[1]][[2]]), #unlist(coverage1Gaussian[[1]][[1]][[3]]), #unlist(coverage1Gaussian[[1]][[1]][[4]]),
                                                     unlist(coverage1Gaussian[[1]][[2]][[1]]), unlist(coverage1Gaussian[[1]][[2]][[2]])),# unlist(coverage1Gaussian[[1]][[2]][[3]])), #unlist(coverage1Gaussian[[1]][[2]][[4]])),
                               # unlist(coverage1Gaussian[[1]][[3]][[1]]), unlist(coverage1Gaussian[[1]][[3]][[2]]), unlist(coverage1Gaussian[[1]][[3]][[3]]), 
                               # unlist(coverage1Gaussian[[2]][[1]][[1]]), unlist(coverage1Gaussian[[2]][[1]][[2]]), unlist(coverage1Gaussian[[2]][[1]][[3]]),
                               # unlist(coverage1Gaussian[[2]][[2]][[1]]), unlist(coverage1Gaussian[[2]][[2]][[2]]), unlist(coverage1Gaussian[[2]][[2]][[3]]),
                               # unlist(coverage1Gaussian[[2]][[3]][[1]]), unlist(coverage1Gaussian[[2]][[3]][[2]]), unlist(coverage1Gaussian[[2]][[3]][[3]]), 
                               # unlist(coverage1Gaussian[[3]][[1]][[1]]), unlist(coverage1Gaussian[[3]][[1]][[2]]), unlist(coverage1Gaussian[[3]][[1]][[3]]),
                               # unlist(coverage1Gaussian[[3]][[2]][[1]]), unlist(coverage1Gaussian[[3]][[2]][[2]]), unlist(coverage1Gaussian[[3]][[2]][[3]]),
                               # unlist(coverage1Gaussian[[3]][[3]][[1]]), unlist(coverage1Gaussian[[3]][[3]][[2]]), unlist(coverage1Gaussian[[3]][[3]][[3]])),
                               coverage2Gaussian = c(unlist(coverage2Gaussian[[1]][[1]][[1]]), unlist(coverage2Gaussian[[1]][[1]][[2]]), #unlist(coverage2Gaussian[[1]][[1]][[3]]), #unlist(coverage2Gaussian[[1]][[1]][[4]]),
                                                     unlist(coverage2Gaussian[[1]][[2]][[1]]), unlist(coverage2Gaussian[[1]][[2]][[2]])),# unlist(coverage2Gaussian[[1]][[2]][[3]])), #unlist(coverage2Gaussian[[1]][[2]][[4]])),
                               # unlist(coverage2Gaussian[[1]][[3]][[1]]), unlist(coverage2Gaussian[[1]][[3]][[2]]), unlist(coverage2Gaussian[[1]][[3]][[3]]), 
                               # unlist(coverage2Gaussian[[2]][[1]][[1]]), unlist(coverage2Gaussian[[2]][[1]][[2]]), unlist(coverage2Gaussian[[2]][[1]][[3]]),
                               # unlist(coverage2Gaussian[[2]][[2]][[1]]), unlist(coverage2Gaussian[[2]][[2]][[2]]), unlist(coverage2Gaussian[[2]][[2]][[3]]),
                               # unlist(coverage2Gaussian[[2]][[3]][[1]]), unlist(coverage2Gaussian[[2]][[3]][[2]]), unlist(coverage2Gaussian[[2]][[3]][[3]]), 
                               # unlist(coverage2Gaussian[[3]][[1]][[1]]), unlist(coverage2Gaussian[[3]][[1]][[2]]), unlist(coverage2Gaussian[[3]][[1]][[3]]),
                               # unlist(coverage2Gaussian[[3]][[2]][[1]]), unlist(coverage2Gaussian[[3]][[2]][[2]]), unlist(coverage2Gaussian[[3]][[2]][[3]]),
                               # unlist(coverage2Gaussian[[3]][[3]][[1]]), unlist(coverage2Gaussian[[3]][[3]][[2]]), unlist(coverage2Gaussian[[3]][[3]][[3]])),
                               factorGaussian = c(unlist(factorGaussian[[1]][[1]][[1]]), unlist(factorGaussian[[1]][[1]][[2]]),# unlist(factorGaussian[[1]][[1]][[3]]), #unlist(factorGaussian[[1]][[1]][[4]]),
                                                  unlist(factorGaussian[[1]][[2]][[1]]), unlist(factorGaussian[[1]][[2]][[2]])),# unlist(factorGaussian[[1]][[2]][[3]])), #unlist(factorGaussian[[1]][[2]][[4]])),
                               # unlist(factorGaussian[[1]][[3]][[1]]), unlist(factorGaussian[[1]][[3]][[2]]), unlist(factorGaussian[[1]][[3]][[3]]), 
                               # unlist(factorGaussian[[2]][[1]][[1]]), unlist(factorGaussian[[2]][[1]][[2]]), unlist(factorGaussian[[2]][[1]][[3]]),
                               # unlist(factorGaussian[[2]][[2]][[1]]), unlist(factorGaussian[[2]][[2]][[2]]), unlist(factorGaussian[[2]][[2]][[3]]),
                               # unlist(factorGaussian[[2]][[3]][[1]]), unlist(factorGaussian[[2]][[3]][[2]]), unlist(factorGaussian[[2]][[3]][[3]]), 
                               # unlist(factorGaussian[[3]][[1]][[1]]), unlist(factorGaussian[[3]][[1]][[2]]), unlist(factorGaussian[[3]][[1]][[3]]),
                               # unlist(factorGaussian[[3]][[2]][[1]]), unlist(factorGaussian[[3]][[2]][[2]]), unlist(factorGaussian[[3]][[2]][[3]]),
                               # unlist(factorGaussian[[3]][[3]][[1]]), unlist(factorGaussian[[3]][[3]][[2]]), unlist(factorGaussian[[3]][[3]][[3]])),
                               rangeGaussian = c(unlist(rangeGaussian[[1]][[1]][[1]]), unlist(rangeGaussian[[1]][[1]][[2]]),# unlist(rangeGaussian[[1]][[1]][[3]]), #unlist(rangeGaussian[[1]][[1]][[4]]),
                                                 unlist(rangeGaussian[[1]][[2]][[1]]), unlist(rangeGaussian[[1]][[2]][[2]])),# unlist(rangeGaussian[[1]][[2]][[3]])), #unlist(rangeGaussian[[1]][[2]][[4]])),
                               # unlist(rangeGaussian[[1]][[3]][[1]]), unlist(rangeGaussian[[1]][[3]][[2]]), unlist(rangeGaussian[[1]][[3]][[3]]), 
                               # unlist(rangeGaussian[[2]][[1]][[1]]), unlist(rangeGaussian[[2]][[1]][[2]]), unlist(rangeGaussian[[2]][[1]][[3]]),
                               # unlist(rangeGaussian[[2]][[2]][[1]]), unlist(rangeGaussian[[2]][[2]][[2]]), unlist(rangeGaussian[[2]][[2]][[3]]),
                               # unlist(rangeGaussian[[2]][[3]][[1]]), unlist(rangeGaussian[[2]][[3]][[2]]), unlist(rangeGaussian[[2]][[3]][[3]]), 
                               # unlist(rangeGaussian[[3]][[1]][[1]]), unlist(rangeGaussian[[3]][[1]][[2]]), unlist(rangeGaussian[[3]][[1]][[3]]),
                               # unlist(rangeGaussian[[3]][[2]][[1]]), unlist(rangeGaussian[[3]][[2]][[2]]), unlist(rangeGaussian[[3]][[2]][[3]]),  
                               # unlist(rangeGaussian[[3]][[3]][[1]]), unlist(rangeGaussian[[3]][[3]][[2]]), unlist(rangeGaussian[[3]][[3]][[3]])),
                               nuggetVarGaussian = c(unlist(nuggetVarGaussian[[1]][[1]][[1]]), unlist(nuggetVarGaussian[[1]][[1]][[2]]),# unlist(nuggetVarGaussian[[1]][[1]][[3]]), #unlist(nuggetVarGaussian[[1]][[1]][[4]]),
                                                     unlist(nuggetVarGaussian[[1]][[2]][[1]]), unlist(nuggetVarGaussian[[1]][[2]][[2]])),# unlist(nuggetVarGaussian[[1]][[2]][[3]])), #unlist(nuggetVarGaussian[[1]][[2]][[4]])),
                               # unlist(nuggetVarGaussian[[1]][[3]][[1]]), unlist(nuggetVarGaussian[[1]][[3]][[2]]), unlist(nuggetVarGaussian[[1]][[3]][[3]]), 
                               # unlist(nuggetVarGaussian[[2]][[1]][[1]]), unlist(nuggetVarGaussian[[2]][[1]][[2]]), unlist(nuggetVarGaussian[[2]][[1]][[3]]),
                               # unlist(nuggetVarGaussian[[2]][[2]][[1]]), unlist(nuggetVarGaussian[[2]][[2]][[2]]), unlist(nuggetVarGaussian[[2]][[2]][[3]]),
                               # unlist(nuggetVarGaussian[[2]][[3]][[1]]), unlist(nuggetVarGaussian[[2]][[3]][[2]]), unlist(nuggetVarGaussian[[2]][[3]][[3]]), 
                               # unlist(nuggetVarGaussian[[3]][[1]][[1]]), unlist(nuggetVarGaussian[[3]][[1]][[2]]), unlist(nuggetVarGaussian[[3]][[1]][[3]]),
                               # unlist(nuggetVarGaussian[[3]][[2]][[1]]), unlist(nuggetVarGaussian[[3]][[2]][[2]]), unlist(nuggetVarGaussian[[3]][[2]][[3]]),
                               # unlist(nuggetVarGaussian[[3]][[3]][[1]]), unlist(nuggetVarGaussian[[3]][[3]][[2]]), unlist(nuggetVarGaussian[[3]][[3]][[3]])),
                               errorGaussian = c(rep(unlist(errorGaussian[[1]][[1]][[1]]), length(unlist(CRPS1Gaussian[[1]][[1]][[1]]))), rep(unlist(errorGaussian[[1]][[1]][[2]]), length(unlist(CRPS1Gaussian[[1]][[1]][[2]]))), #rep(unlist(errorGaussian[[1]][[1]][[3]]), length(unlist(CRPS1Gaussian[[1]][[1]][[3]]))),
                                                 #rep(unlist(errorGaussian[[1]][[1]][[4]]), length(unlist(CRPS1Gaussian[[1]][[1]][[4]]))),
                                                 rep(unlist(errorGaussian[[1]][[2]][[1]]), length(unlist(CRPS1Gaussian[[1]][[2]][[1]]))), rep(unlist(errorGaussian[[1]][[2]][[2]]), length(unlist(CRPS1Gaussian[[1]][[2]][[2]])))))#, rep(unlist(errorGaussian[[1]][[2]][[3]]), length(unlist(CRPS1Gaussian[[1]][[2]][[3]])))))#,
#rep(unlist(errorGaussian[[1]][[2]][[4]]), length(unlist(CRPS1Gaussian[[1]][[2]][[4]])))))
# rep(unlist(errorGaussian[[1]][[3]][[1]]), length(unlist(CRPS1Gaussian[[1]][[3]][[1]]))), rep(unlist(errorGaussian[[1]][[3]][[2]]), length(unlist(CRPS1Gaussian[[1]][[3]][[2]]))), rep(unlist(errorGaussian[[1]][[3]][[3]]), length(unlist(CRPS1Gaussian[[1]][[3]][[3]]))), 
# rep(unlist(errorGaussian[[2]][[1]][[1]]), length(unlist(CRPS1Gaussian[[2]][[1]][[1]]))), rep(unlist(errorGaussian[[2]][[1]][[2]]), length(unlist(CRPS1Gaussian[[2]][[1]][[2]]))), rep(unlist(errorGaussian[[1]][[3]][[3]]), length(unlist(CRPS1Gaussian[[2]][[1]][[3]]))),
# rep(unlist(errorGaussian[[2]][[2]][[1]]), length(unlist(CRPS1Gaussian[[2]][[2]][[1]]))), rep(unlist(errorGaussian[[2]][[2]][[2]]), length(unlist(CRPS1Gaussian[[2]][[2]][[2]]))), rep(unlist(errorGaussian[[2]][[2]][[3]]), length(unlist(CRPS1Gaussian[[2]][[2]][[3]]))),
# rep(unlist(errorGaussian[[2]][[3]][[1]]), length(unlist(CRPS1Gaussian[[2]][[3]][[1]]))), rep(unlist(errorGaussian[[2]][[3]][[2]]), length(unlist(CRPS1Gaussian[[2]][[3]][[2]]))), rep(unlist(errorGaussian[[2]][[2]][[3]]), length(unlist(CRPS1Gaussian[[2]][[3]][[3]]))), 
# rep(unlist(errorGaussian[[3]][[1]][[1]]), length(unlist(CRPS1Gaussian[[3]][[1]][[1]]))), rep(unlist(errorGaussian[[3]][[1]][[2]]), length(unlist(CRPS1Gaussian[[3]][[1]][[2]]))), rep(unlist(errorGaussian[[3]][[1]][[3]]), length(unlist(CRPS1Gaussian[[3]][[1]][[3]]))),
# rep(unlist(errorGaussian[[3]][[2]][[1]]), length(unlist(CRPS1Gaussian[[3]][[2]][[1]]))), rep(unlist(errorGaussian[[3]][[2]][[2]]), length(unlist(CRPS1Gaussian[[3]][[2]][[2]]))), rep(unlist(errorGaussian[[3]][[2]][[3]]), length(unlist(CRPS1Gaussian[[3]][[2]][[3]]))),
# rep(unlist(errorGaussian[[3]][[3]][[1]]), length(unlist(CRPS1Gaussian[[3]][[3]][[1]]))), rep(unlist(errorGaussian[[3]][[3]][[2]]), length(unlist(CRPS1Gaussian[[3]][[3]][[2]]))), rep(unlist(errorGaussian[[3]][[3]][[3]]), length(unlist(CRPS1Gaussian[[3]][[3]][[3]])))))




#   THE PLOTS


# Percent differences in CRPS (Gaussian)

graphInputGaussian_CRPSdifference = data.frame(range = graphMainGaussian$range, factor = graphMainGaussian$factor,
                                               difference = 100*(graphMainGaussian$CRPS2- graphMainGaussian$CRPS1)/graphMainGaussian$CRPS1)

graphInputGaussian_CRPSdifference$factor[graphInputGaussian_CRPSdifference$factor==1]<-"DHS"
graphInputGaussian_CRPSdifference$factor[graphInputGaussian_CRPSdifference$factor==4]<-"4xDHS"

g1 <- ggplot(graphInputGaussian_CRPSdifference, aes(x = fct_inorder(factor), fill = factor, y = difference)) +
  facet_grid(~ range)+geom_boxplot(position = "dodge2", outlier.alpha = 0.1, alpha=0.2) + #ggtitle("CRPS") +
  theme_bw() + scale_fill_grey(start = 0, end = .9)+ theme(legend.position = "none") + theme(axis.title.x=element_blank()) +
  #theme(axis.title.y=element_blank()) + 
  theme(axis.text.y = element_text(size = 30), axis.text.x = element_text(size = 30)) + 
  theme(text = element_text(size=rel(3.5)), strip.text.x = element_text(size=rel(12)), strip.text.y = element_text(size=rel(5))) +
  theme(axis.title.x=element_text(size = rel(10), margin=margin(30,0,0,0))) + xlab("Jittering Scheme") +
  theme(axis.title.y=element_text(size = rel(10), angle = 0))+ylab("Relative (%)")+theme(axis.title.y = element_text(angle=90))+
  theme(panel.margin.x=unit(0, "lines")) + geom_hline(yintercept = 0, linetype="dashed", color = "red", size = 1)

g1


ggsave("CRPSDiffGaussTRUE.png", path = "~/directory")

# Pair wise absolute differences in Log-Scores  (Gaussian)

graphInputGaussian_Logscoresdifference = data.frame(range = graphMainGaussian$range, factor = graphMainGaussian$factor,
                                                       difference = graphMainGaussian$Logscores2- graphMainGaussian$Logscores1)


graphInputGaussian_Logscoresdifference$factor[graphInputGaussian_Logscoresdifference$factor==1]<-"DHS"
graphInputGaussian_Logscoresdifference$factor[graphInputGaussian_Logscoresdifference$factor==4]<-"4xDHS"


g2 <- ggplot(graphInputGaussian_Logscoresdifference, aes(x = fct_inorder(factor), fill = factor, y = difference)) +
  facet_grid(~ range)+geom_boxplot(position = "dodge2", outlier.alpha = 0.1, alpha=0.2) + #ggtitle("Log-Score") +
  theme_bw() + scale_fill_grey(start = 0, end = .9)+ theme(legend.position = "none") + theme(axis.title.x=element_blank()) +
  #theme(axis.title.y=element_blank()) + 
  theme(axis.text.y = element_text(size = 30), axis.text.x = element_text(size = 30)) + 
  theme(text = element_text(size=rel(3.5)), strip.text.x = element_text(size=rel(12)), strip.text.y = element_text(size=rel(5))) +
  theme(axis.title.x=element_text(size = rel(10), margin=margin(30,0,0,0))) + xlab("Jittering Scheme") +
  theme(axis.title.y=element_text(size = rel(10), angle = 0))+ylab("Absolute")+theme(axis.title.y = element_text(angle=90))+
  theme(panel.margin.x=unit(0, "lines")) + geom_hline(yintercept = 0, linetype="dashed", color = "red", size = 1)

g2

ggsave("LogScoreDiffGaussTRUE.png", path = "~/directory")


# Percent differences in CRPS (Binomial)

graphInputBinomial_CRPSdifference = data.frame(range = graphMainBinomial$range, factor = graphMainBinomial$factor,
                                               difference = 100*(graphMainBinomial$CRPS2- graphMainBinomial$CRPS1)/graphMainBinomial$CRPS1)

graphInputBinomial_CRPSdifference$factor[graphInputBinomial_CRPSdifference$factor==1]<-"DHS"
graphInputBinomial_CRPSdifference$factor[graphInputBinomial_CRPSdifference$factor==4]<-"4xDHS"

g3 <- ggplot(graphInputBinomial_CRPSdifference, aes(x = fct_inorder(factor), fill = factor, y = difference)) +
  facet_grid(~ range)+geom_boxplot(position = "dodge2", outlier.alpha = 0.1, alpha=0.2) + #ggtitle("CRPS") +
  theme_bw() + scale_fill_grey(start = 0, end = .9)+ theme(legend.position = "none") + theme(axis.title.x=element_blank()) +
  #theme(axis.title.y=element_blank()) + 
  theme(axis.text.y = element_text(size = 30), axis.text.x = element_text(size = 30)) + 
  theme(text = element_text(size=rel(3.5)), strip.text.x = element_text(size=rel(12)), strip.text.y = element_text(size=rel(5))) +
  theme(axis.title.x=element_text(size = rel(10), margin=margin(30,0,0,0))) + xlab("Jittering Scheme") +
  theme(axis.title.y=element_text(size = rel(10), angle = 0))+ylab("Relative (%)")+theme(axis.title.y = element_text(angle=90))+
  theme(panel.margin.x=unit(0, "lines")) + geom_hline(yintercept = 0, linetype="dashed", color = "red", size = 1)

g3

ggsave("CRPSDiffBinomTRUE.png", path = "~/directory")

# Pair wise absolute differences in Log-Scores  (Binomial)


graphInputBinomial_Logscoresdifference = data.frame(range = graphMainBinomial$range, factor = graphMainBinomial$factor,
                                                    difference = graphMainBinomial$Logscores2- graphMainBinomial$Logscores1)

graphInputBinomial_Logscoresdifference$factor[graphInputBinomial_Logscoresdifference$factor==1]<-"DHS"
graphInputBinomial_Logscoresdifference$factor[graphInputBinomial_Logscoresdifference$factor==4]<-"4xDHS"


g4 <- ggplot(graphInputBinomial_Logscoresdifference, aes(x = fct_inorder(factor), fill = factor, y = difference)) +
  facet_grid(~ range)+geom_boxplot(position = "dodge2", outlier.alpha = 0.1, alpha=0.2) + #ggtitle("Log-Score") +
  theme_bw() + scale_fill_grey(start = 0, end = .9)+ theme(legend.position = "none") + theme(axis.title.x=element_blank()) +
  #theme(axis.title.y=element_blank()) + 
  theme(axis.text.y = element_text(size = 30), axis.text.x = element_text(size = 30)) + 
  theme(text = element_text(size=rel(3.5)), strip.text.x = element_text(size=rel(12)), strip.text.y = element_text(size=rel(5))) +
  theme(axis.title.x=element_text(size = rel(10), margin=margin(30,0,0,0))) + xlab("Jittering Scheme") +
  theme(axis.title.y=element_text(size = rel(10), angle = 0))+ylab("Absolute")+theme(axis.title.y = element_text(angle=90))+
  theme(panel.margin.x=unit(0, "lines")) + geom_hline(yintercept = 0, linetype="dashed", color = "red", size = 1)


g4

ggsave("LogScoreDiffBinomTRUE.png", path = "~/directory")

################################################################################

#BOXPLOTS OF CRPS AND LOG-SCORES

graphMainBinomial_160 = graphMainBinomial[graphMainBinomial$range=="Range = 160 km",]
graphMainBinomial_340 = graphMainBinomial[graphMainBinomial$range=="Range = 340 km",]

graphMainGaussian_160 = graphMainGaussian[graphMainGaussian$range=="Range = 160 km",]
graphMainGaussian_340 = graphMainGaussian[graphMainGaussian$range=="Range = 340 km",]
# 

######

#Gaussian
#Plots for Range = 160 km

#CRPS

graphInput_160 = data.frame(range = graphMainGaussian_160$range, factor = graphMainGaussian_160$factor,
                            CRPS1 = graphMainGaussian_160$CRPS1, CRPS2 = graphMainGaussian_160$CRPS2)


graphInput_160$factor[graphInput_160$factor==1]<-"DHS"
graphInput_160$factor[graphInput_160$factor==4]<-"4xDHS"

colnames(graphInput_160)=c("range", "factor", "S", "J")
#convert it to the long format
dat <- melt(graphInput_160, id.vars = c("range", "factor"))


g1 <- ggplot(data = dat, aes(x = fct_inorder(factor), y = value, group = interaction(factor, variable), fill = variable))+ #+ scale_y_continuous(trans='log10')+
geom_boxplot(position = "dodge2", outlier.alpha = 0.1, alpha=0.2, show.legend = FALSE) + #ggtitle("CRPS")+
facet_grid(. ~ range)+
theme_bw() + scale_x_discrete(limits = c("DHS", "4xDHS"))+
scale_fill_grey(start = 0, end = .9) +
theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15)) + 
  theme(text = element_text(size=rel(2)), strip.text.x = element_text(size=rel(10)), strip.text.y = element_text(size=rel(5))) +
  theme(axis.title.x=element_text(size = rel(8))) + xlab("Jittering Scheme") +
  theme(axis.title.y=element_text(size = rel(8), angle = 0))+ylab("CRPS")+theme(axis.title.y = element_text(angle=90))+
  theme(panel.margin.x=unit(0, "lines")) + #geom_hline(yintercept = 0, linetype="dashed", color = "red", size = 1) +
  theme(legend.text=element_text(size=15))+ #guides(title = labs("CRPS"), title.vjust=3)+
  theme(legend.title = element_text(size = rel(7))) +theme(legend.title = element_blank())
  

g1


#Log Scores
graphInput_160 = data.frame(range = graphMainGaussian_160$range, factor = graphMainGaussian_160$factor,
                            Logscores1 = graphMainGaussian_160$Logscores1, Logscores2 = graphMainGaussian_160$Logscores2)



graphInput_160$factor[graphInput_160$factor==1]<-"DHS"
graphInput_160$factor[graphInput_160$factor==4]<-"4xDHS"

colnames(graphInput_160)=c("range", "factor", "S", "J")
#convert it to the long format
dat <- melt(graphInput_160, id.vars = c("range", "factor"))


g2 <-  ggplot(data = dat, aes(x = fct_inorder(factor), y = value, group = interaction(factor, variable), fill = variable))+ #+ scale_y_continuous(trans='log10')+
  geom_boxplot(position = "dodge2", outlier.alpha = 0.1, alpha=0.2, show.legend = FALSE) + #ggtitle("CRPS")+
  facet_grid(. ~ range)+
  theme_bw() + scale_x_discrete(limits = c("DHS", "4xDHS"))+
  scale_fill_grey(start = 0, end = .9) +
  theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15)) + 
  theme(text = element_text(size=rel(2)), strip.text.x = element_text(size=rel(10)), strip.text.y = element_text(size=rel(5))) +
  theme(axis.title.x=element_text(size = rel(8))) + xlab("Jittering Scheme") +
  theme(axis.title.y=element_text(size = rel(8), angle = 0))+ylab("Log-Score")+theme(axis.title.y = element_text(angle=90))+
  theme(panel.margin.x=unit(0, "lines")) + #geom_hline(yintercept = 0, linetype="dashed", color = "red", size = 1) +
  theme(legend.text=element_text(size=15))+ #guides(title = labs("CRPS"), title.vjust=3)+
  theme(legend.title = element_text(size = rel(7))) +theme(legend.title = element_blank()) 
  
  

g2

#Plots for Range = 340 km


#CRPS
graphInput_340 = data.frame(range = graphMainGaussian_340$range, factor = graphMainGaussian_340$factor,
                            CRPS1 = graphMainGaussian_340$CRPS1, CRPS2 = graphMainGaussian_340$CRPS2)

graphInput_340$factor[graphInput_340$factor==1]<-"DHS"
graphInput_340$factor[graphInput_340$factor==4]<-"4xDHS"

colnames(graphInput_340)=c("range", "factor", "S", "J")

#convert it to the long format
dat <- melt(graphInput_340, id.vars = c("range", "factor"))


g3 <- ggplot(data = dat, aes(x = fct_inorder(factor), y = value, group = interaction(factor, variable), fill = variable)) + #+ scale_y_continuous(trans='log10')+
geom_boxplot(position = "dodge2", outlier.alpha = 0.1, alpha=0.2, show.legend = FALSE)+ #+ ggtitle("CRPS")+
facet_grid(. ~ range)+
theme_bw() + scale_x_discrete(limits = c("DHS", "4xDHS"))+
  scale_fill_grey(start = 0, end = .9) +
  theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15)) + 
  theme(text = element_text(size=rel(2)), strip.text.x = element_text(size=rel(10)), strip.text.y = element_text(size=rel(5))) +
  theme(axis.title.x=element_text(size = rel(8))) + xlab("Jittering Scheme") +
  theme(axis.title.y=element_text(size = rel(8), angle = 0))+ylab("CRPS")+theme(axis.title.y = element_text(angle=90))+
  theme(panel.margin.x=unit(0, "lines")) + #geom_hline(yintercept = 0, linetype="dashed", color = "red", size = 1) +
  theme(legend.text=element_text(size=15))+ #guides(title = labs("CRPS"), title.vjust=3)+
  theme(legend.title = element_text(size = rel(7))) +theme(legend.title = element_blank())

g3


#Log Scores
graphInput_340 = data.frame(range = graphMainGaussian_340$range, factor = graphMainGaussian_340$factor,
                            Logscores1 = graphMainGaussian_340$Logscores1, Logscores2 = graphMainGaussian_340$Logscores2)


graphInput_340$factor[graphInput_340$factor==1]<-"DHS"
graphInput_340$factor[graphInput_340$factor==4]<-"4xDHS"

colnames(graphInput_340)=c("range", "factor", "S", "J")


#convert it to the long format
dat <- melt(graphInput_340, id.vars = c("range", "factor"))


g4 <- ggplot(data = dat, aes(x = fct_inorder(factor), y = value, group = interaction(factor, variable), fill = variable))+
geom_boxplot(position = "dodge2", outlier.alpha = 0.1, alpha=0.2, show.legend = FALSE) +#+ ggtitle("Log Scores")+
facet_grid(. ~ range)+
theme_bw() + scale_x_discrete(limits = c("DHS", "4xDHS"))+
  scale_fill_grey(start = 0, end = .9) +
  theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15)) + 
  theme(text = element_text(size=rel(2)), strip.text.x = element_text(size=rel(10)), strip.text.y = element_text(size=rel(5))) +
  theme(axis.title.x=element_text(size = rel(8))) + xlab("Jittering Scheme") +
  theme(axis.title.y=element_text(size = rel(8), angle = 0))+ylab("Log-Score")+theme(axis.title.y = element_text(angle=90))+
  theme(panel.margin.x=unit(0, "lines")) + #geom_hline(yintercept = 0, linetype="dashed", color = "red", size = 1) +
  theme(legend.text=element_text(size=15))+ #guides(title = labs("CRPS"), title.vjust=3)+
  theme(legend.title = element_text(size = rel(7))) +theme(legend.title = element_blank()) 
g4

#Place the plots into a matrix


plot<-plot_grid(g1, g2, g3, g4, align='vh', vjust=1, scale = 1)
plot
ggsave("TRUE_GaussOnly.png", path = "~/directory")

################################################################################


#Binomial
#Plots for Range = 160 km

#CRPS

graphInput_160 = data.frame(range = graphMainBinomial_160$range, factor = graphMainBinomial_160$factor,
                            CRPS1 = graphMainBinomial_160$CRPS1, CRPS2 = graphMainBinomial_160$CRPS2)


graphInput_160$factor[graphInput_160$factor==1]<-"DHS"
graphInput_160$factor[graphInput_160$factor==4]<-"4xDHS"

colnames(graphInput_160)=c("range", "factor", "S", "J")


#convert it to the long format
dat <- melt(graphInput_160, id.vars = c("range", "factor"))

g5 <- ggplot(data = dat, aes(x = fct_inorder(factor), y = value, group = interaction(factor, variable), fill = variable)) +#+ scale_y_continuous(trans='log10')+
geom_boxplot(position = "dodge2", outlier.alpha = 0.1, alpha=0.2, show.legend = FALSE) + #ggtitle("CRPS")+
facet_grid(. ~ range)+
theme_bw() + scale_x_discrete(limits = c("DHS", "4xDHS"))+
  scale_fill_grey(start = 0, end = .9) +
  theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15)) + 
  theme(text = element_text(size=rel(2)), strip.text.x = element_text(size=rel(10)), strip.text.y = element_text(size=rel(5))) +
  theme(axis.title.x=element_text(size = rel(8))) + xlab("Jittering Scheme") +
  theme(axis.title.y=element_text(size = rel(8), angle = 0))+ylab("CRPS")+theme(axis.title.y = element_text(angle=90))+
  theme(panel.margin.x=unit(0, "lines")) + #geom_hline(yintercept = 0, linetype="dashed", color = "red", size = 1) +
  theme(legend.text=element_text(size=15))+ #guides(title = labs("CRPS"), title.vjust=3)+
  theme(legend.title = element_text(size = rel(7))) +theme(legend.title = element_blank())
g5

#Log Scores
graphInput_160 = data.frame(range = graphMainBinomial_160$range, factor = graphMainBinomial_160$factor,
                            Logscores1 = graphMainBinomial_160$Logscores1, Logscores2 = graphMainBinomial_160$Logscores2)


graphInput_160$factor[graphInput_160$factor==1]<-"DHS"
graphInput_160$factor[graphInput_160$factor==4]<-"4xDHS"

colnames(graphInput_160)=c("range", "factor", "S", "J")

#convert it to the long format
dat <- melt(graphInput_160, id.vars = c("range", "factor"))


g6 <- ggplot(data = dat, aes(x = fct_inorder(factor), y = value, group = interaction(factor, variable), fill = variable))+
geom_boxplot(position = "dodge2", outlier.alpha = 0.1, alpha=0.2, show.legend = FALSE) +# ggtitle("Log Scores")
facet_grid(. ~ range)+
theme_bw() + scale_x_discrete(limits = c("DHS", "4xDHS"))+
  scale_fill_grey(start = 0, end = .9) +
  theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15)) + 
  theme(text = element_text(size=rel(2)), strip.text.x = element_text(size=rel(10)), strip.text.y = element_text(size=rel(5))) +
  theme(axis.title.x=element_text(size = rel(8))) + xlab("Jittering Scheme") +
  theme(axis.title.y=element_text(size = rel(8), angle = 0))+ylab("Log-Score")+theme(axis.title.y = element_text(angle=90))+
  theme(panel.margin.x=unit(0, "lines")) + #geom_hline(yintercept = 0, linetype="dashed", color = "red", size = 1) +
  theme(legend.text=element_text(size=15))+ #guides(title = labs("CRPS"), title.vjust=3)+
  theme(legend.title = element_text(size = rel(7))) +theme(legend.title = element_blank()) 
g6

#Plots for Range = 340 km


#CRPS
graphInput_340 = data.frame(range = graphMainBinomial_340$range, factor = graphMainBinomial_340$factor,
                            CRPS1 = graphMainBinomial_340$CRPS1, CRPS2 = graphMainBinomial_340$CRPS2)

graphInput_340$factor[graphInput_340$factor==1]<-"DHS"
graphInput_340$factor[graphInput_340$factor==4]<-"4xDHS"

colnames(graphInput_340)=c("range", "factor", "S", "J")

#convert it to the long format
dat <- melt(graphInput_340, id.vars = c("range", "factor"))

g7 <- ggplot(data = dat, aes(x = fct_inorder(factor), y = value, group = interaction(factor, variable), fill = variable))+ #+ scale_y_continuous(trans='log10')
geom_boxplot(position = "dodge2", outlier.alpha = 0.1, alpha=0.2, show.legend = FALSE)+ #+ ggtitle("CRPS")+
facet_grid(. ~ range)+
theme_bw() + scale_x_discrete(limits = c("DHS", "4xDHS"))+
  scale_fill_grey(start = 0, end = .9) +
  theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15)) + 
  theme(text = element_text(size=rel(2)), strip.text.x = element_text(size=rel(10)), strip.text.y = element_text(size=rel(5))) +
  theme(axis.title.x=element_text(size = rel(8))) + xlab("Jittering Scheme") +
  theme(axis.title.y=element_text(size = rel(8), angle = 0))+ylab("CRPS")+theme(axis.title.y = element_text(angle=90))+
  theme(panel.margin.x=unit(0, "lines")) + #geom_hline(yintercept = 0, linetype="dashed", color = "red", size = 1) +
  theme(legend.text=element_text(size=15))+ #guides(title = labs("CRPS"), title.vjust=3)+
  theme(legend.title = element_text(size = rel(7))) +theme(legend.title = element_blank())
g7


#Log Scores
graphInput_340 = data.frame(range = graphMainBinomial_340$range, factor = graphMainBinomial_340$factor,
                            Logscores1 = graphMainBinomial_340$Logscores1, Logscores2 = graphMainBinomial_340$Logscores2)

graphInput_340$factor[graphInput_340$factor==1]<-"DHS"
graphInput_340$factor[graphInput_340$factor==4]<-"4xDHS"

colnames(graphInput_340)=c("range", "factor", "S", "J")

#convert it to the long format
dat <- melt(graphInput_340, id.vars = c("range", "factor"))


g8 <- ggplot(data = dat, aes(x = fct_inorder(factor), y = value, group = interaction(factor, variable), fill = variable))+
geom_boxplot(position = "dodge2", outlier.alpha = 0.1, alpha=0.2, show.legend = FALSE) +#+ ggtitle("Log Scores")
facet_grid(. ~ range)+
  theme_bw() + scale_x_discrete(limits = c("DHS", "4xDHS"))+
  scale_fill_grey(start = 0, end = .9) +
  theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15)) + 
  theme(text = element_text(size=rel(2)), strip.text.x = element_text(size=rel(10)), strip.text.y = element_text(size=rel(5))) +
  theme(axis.title.x=element_text(size = rel(8))) + xlab("Jittering Scheme") +
  theme(axis.title.y=element_text(size = rel(8), angle = 0))+ylab("Log-Score")+theme(axis.title.y = element_text(angle=90))+
  theme(panel.margin.x=unit(0, "lines")) + #geom_hline(yintercept = 0, linetype="dashed", color = "red", size = 1) +
  theme(legend.text=element_text(size=15))+ #guides(title = labs("CRPS"), title.vjust=3)+
  theme(legend.title = element_text(size = rel(7))) +theme(legend.title = element_blank()) 
g8

#Place the plots into a matrix

plot<-plot_grid(g5, g6, g7, g8, align='vh', vjust=1, scale = 1)
plot

ggsave("TRUE_BinomialOnly.png", path = "~/directory")


########################  TABLES  ##############################################

# 
# # # Result files
load("~/directory/fullResults.RData")  # fullResults.RData file for the Gaussian observation model
fittedGaussian = fitted
rm(fitted)

load("~/directory/fullResults.RData")  #fullResults.RData file for the Binomial observation model
fittedBinomial = fitted
rm(fitted)

nLikelihood = 1
nRange = 2
nScale = 2
nSim = 50



#Tables for the Bias and Credible Intervals of the parameter estimates
# 

# B : stands for the relative bias (except for the intercept. It is just bias for the intercept, not the relative bias.)
#CI stands for the length of the credible interval

# #Standard model with Gaussian Simulations
B_range_Gauss_standard = list()
CI_range_Gauss_standard = list()

B_nuggVar_Gauss_standard = list()
CI_nuggVar_Gauss_standard = list()

B_margVar_Gauss_standard = list()
CI_margVar_Gauss_standard = list()

B_alpha_Gauss_standard = list()
CI_alpha_Gauss_standard = list()
# ################################################################################
# #Jittering Accounted model with Gaussian Simulations
# 
B_range_Gauss_Jitt = list()
CI_range_Gauss_Jitt = list()

B_nuggVar_Gauss_Jitt = list()
CI_nuggVar_Gauss_Jitt = list()

B_margVar_Gauss_Jitt = list()
CI_margVar_Gauss_Jitt = list()

B_alpha_Gauss_Jitt = list()
CI_alpha_Gauss_Jitt = list()


#Standard model with Binomial Simulations
B_range_Binom_standard = list()
CI_range_Binom_standard = list()
# 
# B_nuggVar_Binom_standard = list()
# CI_nuggVar_Binom_standard = list()

B_margVar_Binom_standard = list()
CI_margVar_Binom_standard = list()

B_alpha_Binom_standard = list()
CI_alpha_Binom_standard = list()
################################################################################
#Jittering Accounted model with Binomial Simulations

B_range_Binom_Jitt = list()
CI_range_Binom_Jitt = list()

# B_nuggVar_Binom_Jitt = list()
# CI_nuggVar_Binom_Jitt = list()

B_margVar_Binom_Jitt = list()
CI_margVar_Binom_Jitt = list()
# 
B_alpha_Binom_Jitt = list()
CI_alpha_Binom_Jitt = list()

for (i in 1:nLikelihood){
  #Standard model with Gaussian Simulations
  # 
  B_range_Gauss_standardTemp1 = list()
  CI_range_Gauss_standardTemp1 = list()

  B_nuggVar_Gauss_standardTemp1 = list()
  CI_nuggVar_Gauss_standardTemp1 = list()

  B_margVar_Gauss_standardTemp1 = list()
  CI_margVar_Gauss_standardTemp1 = list()

  B_alpha_Gauss_standardTemp1 = list()
  CI_alpha_Gauss_standardTemp1 = list()
  # 
  # #Jittering Accounted model with Gaussian Simulations
  # 
  B_range_Gauss_JittTemp1 = list()
  CI_range_Gauss_JittTemp1 = list()

  B_nuggVar_Gauss_JittTemp1 = list()
  CI_nuggVar_Gauss_JittTemp1 = list()

  B_margVar_Gauss_JittTemp1 = list()
  CI_margVar_Gauss_JittTemp1 = list()

  B_alpha_Gauss_JittTemp1 = list()
  CI_alpha_Gauss_JittTemp1 = list()
  # 
  #Standard model with Binomial Simulations
  
   B_range_Binom_standardTemp1 = list()
   CI_range_Binom_standardTemp1 = list()
  #
   # B_nuggVar_Binom_standardTemp1 = list()
   # CI_nuggVar_Binom_standardTemp1 = list()
  #
   B_margVar_Binom_standardTemp1 = list()
   CI_margVar_Binom_standardTemp1 = list()
  #
   B_alpha_Binom_standardTemp1 = list()
   CI_alpha_Binom_standardTemp1 = list()

  #Jittering Accounted model with Binomial Simulations

   B_range_Binom_JittTemp1 = list()
   CI_range_Binom_JittTemp1 = list()
  #
   # B_nuggVar_Binom_JittTemp1 = list()
   # CI_nuggVar_Binom_JittTemp1 = list()
  #
   B_margVar_Binom_JittTemp1 = list()
   CI_margVar_Binom_JittTemp1 = list()
  #
   B_alpha_Binom_JittTemp1 = list()
   CI_alpha_Binom_JittTemp1 = list()
  #
  for (j in 1:nRange){
    #Standard model with Gaussian Simulations
    
    B_range_Gauss_standardTemp2 = list()
    CI_range_Gauss_standardTemp2 = list()

    B_nuggVar_Gauss_standardTemp2 = list()
    CI_nuggVar_Gauss_standardTemp2 = list()

    B_margVar_Gauss_standardTemp2 = list()
    CI_margVar_Gauss_standardTemp2 = list()

    B_alpha_Gauss_standardTemp2 = list()
    CI_alpha_Gauss_standardTemp2 = list()
    # 
    # #Jittering Accounted model with Gaussian Simulations
    # 
    B_range_Gauss_JittTemp2 = list()
    CI_range_Gauss_JittTemp2 = list()

    B_nuggVar_Gauss_JittTemp2 = list()
    CI_nuggVar_Gauss_JittTemp2 = list()

    B_margVar_Gauss_JittTemp2 = list()
    CI_margVar_Gauss_JittTemp2 = list()

    B_alpha_Gauss_JittTemp2 = list()
    CI_alpha_Gauss_JittTemp2 = list()
    
    #Standard model with Binomial Simulations
    
    B_range_Binom_standardTemp2 = list()
    CI_range_Binom_standardTemp2 = list()

    # B_nuggVar_Binom_standardTemp2 = list()
    # CI_nuggVar_Binom_standardTemp2 = list()

    B_margVar_Binom_standardTemp2 = list()
    CI_margVar_Binom_standardTemp2 = list()

    B_alpha_Binom_standardTemp2 = list()
    CI_alpha_Binom_standardTemp2 = list()
    
    #Jittering Accounted model with Binomial Simulations
    
    B_range_Binom_JittTemp2 = list()
    CI_range_Binom_JittTemp2 = list()

    # B_nuggVar_Binom_JittTemp2 = list()
    # CI_nuggVar_Binom_JittTemp2 = list()

    B_margVar_Binom_JittTemp2 = list()
    CI_margVar_Binom_JittTemp2 = list()

    B_alpha_Binom_JittTemp2 = list()
    CI_alpha_Binom_JittTemp2 = list()
    # 
    for (h in 1:nScale){
      #Standard model with Gaussian Simulations
      # 
      B_range_Gauss_standardTemp3 = list()
      CI_range_Gauss_standardTemp3 = list()

      B_nuggVar_Gauss_standardTemp3 = list()
      CI_nuggVar_Gauss_standardTemp3 = list()

      B_margVar_Gauss_standardTemp3 = list()
      CI_margVar_Gauss_standardTemp3 = list()

      B_alpha_Gauss_standardTemp3 = list()
      CI_alpha_Gauss_standardTemp3 = list()

      # #Jittering Accounted model with Gaussian Simulations
      # 
      B_range_Gauss_JittTemp3 = list()
      CI_range_Gauss_JittTemp3 = list()

      B_nuggVar_Gauss_JittTemp3 = list()
      CI_nuggVar_Gauss_JittTemp3 = list()

      B_margVar_Gauss_JittTemp3 = list()
      CI_margVar_Gauss_JittTemp3 = list()

      B_alpha_Gauss_JittTemp3 = list()
      CI_alpha_Gauss_JittTemp3 = list()
      
      #Standard model with Binomial Simulations
      
      B_range_Binom_standardTemp3 = list()
      CI_range_Binom_standardTemp3 = list()

      # B_nuggVar_Binom_standardTemp3 = list()
      # CI_nuggVar_Binom_standardTemp3 = list()

      B_margVar_Binom_standardTemp3 = list()
      CI_margVar_Binom_standardTemp3 = list()

      B_alpha_Binom_standardTemp3 = list()
      CI_alpha_Binom_standardTemp3 = list()
      
      #Jittering Accounted model with Binomial Simulations
      
      B_range_Binom_JittTemp3 = list()
      CI_range_Binom_JittTemp3 = list()

      # B_nuggVar_Binom_JittTemp3 = list()
      # CI_nuggVar_Binom_JittTemp3 = list()

      B_margVar_Binom_JittTemp3 = list()
      CI_margVar_Binom_JittTemp3 = list()

      B_alpha_Binom_JittTemp3 = list()
      CI_alpha_Binom_JittTemp3 = list()
      
      for (l in 1:nSim){
        #the bias here represents the relative bias except ((estimate-true)/true) for the one for the intercept
        if(class(fittedGaussian[[1]][[1]][[i]][[j]][[h]][[l]])  != "try-error"){
          if (j ==1){
            B_range_Gauss_standardTemp3[[l]] = (fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_standard"]][1,1] - 160)/160
            B_range_Gauss_JittTemp3[[l]] = (fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_JittAcc"]][1,1] - 160)/160
          } else {
            B_range_Gauss_standardTemp3[[l]] = (fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_standard"]][1,1] - 340)/340
            B_range_Gauss_JittTemp3[[l]] = (fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_JittAcc"]][1,1] - 340)/340
          }
          #Standard model with Gaussian Simulations
          CI_range_Gauss_standardTemp3[[l]] = fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_standard"]][1,4]
          
          B_nuggVar_Gauss_standardTemp3[[l]] = (fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_standard"]][3,1] - 0.1)/0.1
          CI_nuggVar_Gauss_standardTemp3[[l]] = fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_standard"]][3,4]
          
          B_margVar_Gauss_standardTemp3[[l]] = (fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_standard"]][2,1] - 1)/1
          CI_margVar_Gauss_standardTemp3[[l]] = fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_standard"]][2,4]
          
          B_alpha_Gauss_standardTemp3[[l]] = fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_standard"]][4,1] - 1
          CI_alpha_Gauss_standardTemp3[[l]] = fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_standard"]][4,4]
          # 
          # #Jittering Accounted model with Gaussian Simulations
          # 
          CI_range_Gauss_JittTemp3[[l]] = fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_JittAcc"]][1,4]
          
          B_nuggVar_Gauss_JittTemp3[[l]] = (fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_JittAcc"]][3,1] - 0.1)/0.1
          CI_nuggVar_Gauss_JittTemp3[[l]] = fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_JittAcc"]][3,4]
          
          B_margVar_Gauss_JittTemp3[[l]] = (fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_JittAcc"]][2,1] - 1)/1
          CI_margVar_Gauss_JittTemp3[[l]] = fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_JittAcc"]][2,4]
          
          B_alpha_Gauss_JittTemp3[[l]] = fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_JittAcc"]][4,1] - 1
          CI_alpha_Gauss_JittTemp3[[l]] = fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_JittAcc"]][4,4]
        }
        
        if(class(fittedBinomial[[1]][[1]][[i]][[j]][[h]][[l]])  != "try-error"){
          if (j ==1){
            B_range_Binom_standardTemp3[[l]] = (fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_standard"]][1,1] - 160)/160
            B_range_Binom_JittTemp3[[l]] = (fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_JittAcc"]][1,1] - 160)/160
          } else {
            B_range_Binom_standardTemp3[[l]] = (fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_standard"]][1,1] - 340)/340
            B_range_Binom_JittTemp3[[l]] = (fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_JittAcc"]][1,1] - 340)/340
          }
          #Standard model with Binomial Simulations
          CI_range_Binom_standardTemp3[[l]] = fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_standard"]][1,4]
          
          # B_nuggVar_Binom_standardTemp3[[l]] = ((fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_standard"]][3,1])^2 - 0.1)/0.1
          # CI_nuggVar_Binom_standardTemp3[[l]] = fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_standard"]][3,4]
          
          B_margVar_Binom_standardTemp3[[l]] = (fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_standard"]][2,1] - 1)/1
          CI_margVar_Binom_standardTemp3[[l]] = fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_standard"]][2,4]
          
          B_alpha_Binom_standardTemp3[[l]] = fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_standard"]][4,1] - 0
          CI_alpha_Binom_standardTemp3[[l]] = fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_standard"]][4,4]
          
          #Jittering Accounted model with Binomial Simulations
          
          CI_range_Binom_JittTemp3[[l]] = fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_JittAcc"]][1,4]
          
          # B_nuggVar_Binom_JittTemp3[[l]] = ((fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_JittAcc"]][3,1])^2 - 0.1)/0.1
          # CI_nuggVar_Binom_JittTemp3[[l]] = fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_JittAcc"]][3,4]
          
          B_margVar_Binom_JittTemp3[[l]] = (fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_JittAcc"]][2,1] - 1)/1
          CI_margVar_Binom_JittTemp3[[l]] = fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_JittAcc"]][2,4]
          
          B_alpha_Binom_JittTemp3[[l]] = fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_JittAcc"]][4,1] - 0
          CI_alpha_Binom_JittTemp3[[l]] = fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_JittAcc"]][4,4]
          
        } else {
          #Standard model with Binomial Simulations
          CI_range_Binom_standardTemp3[[l]] = NA
          
          # B_nuggVar_Binom_standardTemp3[[l]] = ((fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_standard"]][3,1])^2 - 0.1)/0.1
          # CI_nuggVar_Binom_standardTemp3[[l]] = fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_standard"]][3,4]
          
          B_margVar_Binom_standardTemp3[[l]] = NA
          CI_margVar_Binom_standardTemp3[[l]] = NA
          
          B_alpha_Binom_standardTemp3[[l]] = NA
          CI_alpha_Binom_standardTemp3[[l]] = NA
          
          #Jittering Accounted model with Binomial Simulations
          
          CI_range_Binom_JittTemp3[[l]] = NA
          
          # B_nuggVar_Binom_JittTemp3[[l]] = ((fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_JittAcc"]][3,1])^2 - 0.1)/0.1
          # CI_nuggVar_Binom_JittTemp3[[l]] = fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["summParam_JittAcc"]][3,4]
          
          B_margVar_Binom_JittTemp3[[l]] = NA
          CI_margVar_Binom_JittTemp3[[l]] = NA
          
          B_alpha_Binom_JittTemp3[[l]] = NA
          CI_alpha_Binom_JittTemp3[[l]] = NA
          
        } 
      }
      
      #Standard model with Gaussian Simulations
      
      B_range_Gauss_standardTemp2[[h]] = B_range_Gauss_standardTemp3
      CI_range_Gauss_standardTemp2[[h]] = CI_range_Gauss_standardTemp3

      B_nuggVar_Gauss_standardTemp2[[h]] = B_nuggVar_Gauss_standardTemp3
      CI_nuggVar_Gauss_standardTemp2[[h]] = CI_nuggVar_Gauss_standardTemp3

      B_margVar_Gauss_standardTemp2[[h]] = B_margVar_Gauss_standardTemp3
      CI_margVar_Gauss_standardTemp2[[h]] = CI_margVar_Gauss_standardTemp3

      B_alpha_Gauss_standardTemp2[[h]] = B_alpha_Gauss_standardTemp3
      CI_alpha_Gauss_standardTemp2[[h]] = CI_alpha_Gauss_standardTemp3
      
      # #Jittering Accounted model with Gaussian Simulations
      # 
      B_range_Gauss_JittTemp2[[h]] = B_range_Gauss_JittTemp3
      CI_range_Gauss_JittTemp2[[h]] = CI_range_Gauss_JittTemp3

      B_nuggVar_Gauss_JittTemp2[[h]] = B_nuggVar_Gauss_JittTemp3
      CI_nuggVar_Gauss_JittTemp2[[h]] = CI_nuggVar_Gauss_JittTemp3

      B_margVar_Gauss_JittTemp2[[h]] = B_margVar_Gauss_JittTemp3
      CI_margVar_Gauss_JittTemp2[[h]] = CI_margVar_Gauss_JittTemp3

      B_alpha_Gauss_JittTemp2[[h]] = B_alpha_Gauss_JittTemp3
      CI_alpha_Gauss_JittTemp2[[h]] = CI_alpha_Gauss_JittTemp3

      #Standard model with Binomial Simulations
      
      B_range_Binom_standardTemp2[[h]] = B_range_Binom_standardTemp3
      CI_range_Binom_standardTemp2[[h]] = CI_range_Binom_standardTemp3

      # B_nuggVar_Binom_standardTemp2[[h]] = B_nuggVar_Binom_standardTemp3
      # CI_nuggVar_Binom_standardTemp2[[h]] = CI_nuggVar_Binom_standardTemp3

      B_margVar_Binom_standardTemp2[[h]] = B_margVar_Binom_standardTemp3
      CI_margVar_Binom_standardTemp2[[h]] = CI_margVar_Binom_standardTemp3

      B_alpha_Binom_standardTemp2[[h]] = B_alpha_Binom_standardTemp3
      CI_alpha_Binom_standardTemp2[[h]] = CI_alpha_Binom_standardTemp3
      
      #Jittering Accounted model with Binomial Simulations
      
      B_range_Binom_JittTemp2[[h]] = B_range_Binom_JittTemp3
      CI_range_Binom_JittTemp2[[h]] = CI_range_Binom_JittTemp3

      # B_nuggVar_Binom_JittTemp2[[h]] = B_nuggVar_Binom_JittTemp3
      # CI_nuggVar_Binom_JittTemp2[[h]] = CI_nuggVar_Binom_JittTemp3

      B_margVar_Binom_JittTemp2[[h]] = B_margVar_Binom_JittTemp3
      CI_margVar_Binom_JittTemp2[[h]] = CI_margVar_Binom_JittTemp3

      B_alpha_Binom_JittTemp2[[h]] = B_alpha_Binom_JittTemp3
      CI_alpha_Binom_JittTemp2[[h]] = CI_alpha_Binom_JittTemp3

    }
    
    #Standard model with Gaussian Simulations
    
    B_range_Gauss_standardTemp1[[j]] = B_range_Gauss_standardTemp2
    CI_range_Gauss_standardTemp1[[j]] = CI_range_Gauss_standardTemp2

    B_nuggVar_Gauss_standardTemp1[[j]] = B_nuggVar_Gauss_standardTemp2
    CI_nuggVar_Gauss_standardTemp1[[j]] = CI_nuggVar_Gauss_standardTemp2

    B_margVar_Gauss_standardTemp1[[j]] = B_margVar_Gauss_standardTemp2
    CI_margVar_Gauss_standardTemp1[[j]] = CI_margVar_Gauss_standardTemp2

    B_alpha_Gauss_standardTemp1[[j]] = B_alpha_Gauss_standardTemp2
    CI_alpha_Gauss_standardTemp1[[j]] = CI_alpha_Gauss_standardTemp2
    # 
    # #Jittering Accounted model with Gaussian Simulations
    # 
    B_range_Gauss_JittTemp1[[j]] = B_range_Gauss_JittTemp2
    CI_range_Gauss_JittTemp1[[j]] = CI_range_Gauss_JittTemp2

    B_nuggVar_Gauss_JittTemp1[[j]] = B_nuggVar_Gauss_JittTemp2
    CI_nuggVar_Gauss_JittTemp1[[j]] = CI_nuggVar_Gauss_JittTemp2

    B_margVar_Gauss_JittTemp1[[j]] = B_margVar_Gauss_JittTemp2
    CI_margVar_Gauss_JittTemp1[[j]] = CI_margVar_Gauss_JittTemp2

    B_alpha_Gauss_JittTemp1[[j]] = B_alpha_Gauss_JittTemp2
    CI_alpha_Gauss_JittTemp1[[j]] = CI_alpha_Gauss_JittTemp2
    
    #Standard model with Binomial Simulations
    
    B_range_Binom_standardTemp1[[j]] = B_range_Binom_standardTemp2
    CI_range_Binom_standardTemp1[[j]] = CI_range_Binom_standardTemp2

    # B_nuggVar_Binom_standardTemp1[[j]] = B_nuggVar_Binom_standardTemp2
    # CI_nuggVar_Binom_standardTemp1[[j]] = CI_nuggVar_Binom_standardTemp2

    B_margVar_Binom_standardTemp1[[j]] = B_margVar_Binom_standardTemp2
    CI_margVar_Binom_standardTemp1[[j]] = CI_margVar_Binom_standardTemp2

    B_alpha_Binom_standardTemp1[[j]] = B_alpha_Binom_standardTemp2
    CI_alpha_Binom_standardTemp1[[j]] = CI_alpha_Binom_standardTemp2
    
    #Jittering Accounted model with Binomial Simulations
    
    B_range_Binom_JittTemp1[[j]] = B_range_Binom_JittTemp2
    CI_range_Binom_JittTemp1[[j]] = CI_range_Binom_JittTemp2

    # B_nuggVar_Binom_JittTemp1[[j]] = B_nuggVar_Binom_JittTemp2
    # CI_nuggVar_Binom_JittTemp1[[j]] = CI_nuggVar_Binom_JittTemp2

    B_margVar_Binom_JittTemp1[[j]] = B_margVar_Binom_JittTemp2
    CI_margVar_Binom_JittTemp1[[j]] = CI_margVar_Binom_JittTemp2

    B_alpha_Binom_JittTemp1[[j]] = B_alpha_Binom_JittTemp2
    CI_alpha_Binom_JittTemp1[[j]] = CI_alpha_Binom_JittTemp2
    
  }
  #Standard model with Gaussian Simulations
  
  B_range_Gauss_standard[[i]] = B_range_Gauss_standardTemp1
  CI_range_Gauss_standard[[i]] = CI_range_Gauss_standardTemp1

  B_nuggVar_Gauss_standard[[i]] = B_nuggVar_Gauss_standardTemp1
  CI_nuggVar_Gauss_standard[[i]] = CI_nuggVar_Gauss_standardTemp1

  B_margVar_Gauss_standard[[i]] = B_margVar_Gauss_standardTemp1
  CI_margVar_Gauss_standard[[i]] = CI_margVar_Gauss_standardTemp1

  B_alpha_Gauss_standard[[i]] = B_alpha_Gauss_standardTemp1
  CI_alpha_Gauss_standard[[i]] = CI_alpha_Gauss_standardTemp1
  
  # #Jittering Accounted model with Gaussian Simulations
  # 
  B_range_Gauss_Jitt[[i]] = B_range_Gauss_JittTemp1
  CI_range_Gauss_Jitt[[i]] = CI_range_Gauss_JittTemp1

  B_nuggVar_Gauss_Jitt[[i]] = B_nuggVar_Gauss_JittTemp1
  CI_nuggVar_Gauss_Jitt[[i]] = CI_nuggVar_Gauss_JittTemp1

  B_margVar_Gauss_Jitt[[i]] = B_margVar_Gauss_JittTemp1
  CI_margVar_Gauss_Jitt[[i]] = CI_margVar_Gauss_JittTemp1

  B_alpha_Gauss_Jitt[[i]] = B_alpha_Gauss_JittTemp1
  CI_alpha_Gauss_Jitt[[i]] = CI_alpha_Gauss_JittTemp1
  # 
  
  #Standard model with Binomial Simulations
  
  B_range_Binom_standard[[i]] = B_range_Binom_standardTemp1
  CI_range_Binom_standard[[i]] = CI_range_Binom_standardTemp1

  # B_nuggVar_Binom_standard[[i]] = B_nuggVar_Binom_standardTemp1
  # CI_nuggVar_Binom_standard[[i]] = CI_nuggVar_Binom_standardTemp1

  B_margVar_Binom_standard[[i]] = B_margVar_Binom_standardTemp1
  CI_margVar_Binom_standard[[i]] = CI_margVar_Binom_standardTemp1

  B_alpha_Binom_standard[[i]] = B_alpha_Binom_standardTemp1
  CI_alpha_Binom_standard[[i]] = CI_alpha_Binom_standardTemp1
  
  #Jittering Accounted model with Binomial Simulations
  
  B_range_Binom_Jitt[[i]] = B_range_Binom_JittTemp1
  CI_range_Binom_Jitt[[i]] = CI_range_Binom_JittTemp1

  # B_nuggVar_Binom_Jitt[[i]] = B_nuggVar_Binom_JittTemp1
  # CI_nuggVar_Binom_Jitt[[i]] = CI_nuggVar_Binom_JittTemp1

  B_margVar_Binom_Jitt[[i]] = B_margVar_Binom_JittTemp1
  CI_margVar_Binom_Jitt[[i]] = CI_margVar_Binom_JittTemp1

  B_alpha_Binom_Jitt[[i]] = B_alpha_Binom_JittTemp1
  CI_alpha_Binom_Jitt[[i]] = CI_alpha_Binom_JittTemp1
}


# 

################################################################################


#Table 1 : Gaussian Standard Model with 160 km of range
column1 = c("Jittering Factor", "Parameters", "$\\alpha$", "$\\tau$", "$\\sigma^{2}_{N}$", "$\\sigma^{2}_{SF}$")
column2 = c("&", "True", 1, 160, 0.1, 1)
column3 = c(1, "Bias", mean(unlist(B_alpha_Gauss_standard[[1]][[1]][[1]]), na.rm=TRUE), mean(unlist(B_range_Gauss_standard[[1]][[1]][[1]]), na.rm=TRUE),
            mean(unlist(B_nuggVar_Gauss_standard[[1]][[1]][[1]]), na.rm=TRUE), mean(unlist(B_margVar_Gauss_standard[[1]][[1]][[1]]), na.rm=TRUE))

column4 = c(1, "CI Length", mean(unlist(CI_alpha_Gauss_standard[[1]][[1]][[1]]), na.rm=TRUE), mean(unlist(CI_range_Gauss_standard[[1]][[1]][[1]]), na.rm=TRUE), 
            mean(unlist(CI_nuggVar_Gauss_standard[[1]][[1]][[1]]), na.rm=TRUE), mean(unlist(CI_margVar_Gauss_standard[[1]][[1]][[1]]), na.rm=TRUE))

column5 = c(4, "Bias", mean(unlist(B_alpha_Gauss_standard[[1]][[1]][[2]]), na.rm=TRUE), mean(unlist(B_range_Gauss_standard[[1]][[1]][[2]]), na.rm=TRUE), 
            mean(unlist(B_nuggVar_Gauss_standard[[1]][[1]][[2]]), na.rm=TRUE), mean(unlist(B_margVar_Gauss_standard[[1]][[1]][[2]]), na.rm=TRUE))

column6 = c(4, "CI Length", mean(unlist(CI_alpha_Gauss_standard[[1]][[1]][[2]]), na.rm=TRUE), mean(unlist(CI_range_Gauss_standard[[1]][[1]][[2]]), na.rm=TRUE), 
            mean(unlist(CI_nuggVar_Gauss_standard[[1]][[1]][[2]]), na.rm=TRUE), mean(unlist(CI_margVar_Gauss_standard[[1]][[1]][[2]]), na.rm=TRUE))


#Table 2 : Gaussian Jittering Accounted Model with 160 km of range
column7 = c("Jittering Factor", "Parameters", "$\\alpha$", "$\\tau$", "$\\sigma^{2}_{N}$", "$\\sigma^{2}_{SF}$")
column8 = c("&", "True", 1, 160, 0.1, 1)

column9 = c(1, "Bias", mean(unlist(B_alpha_Gauss_Jitt[[1]][[1]][[1]]), na.rm=TRUE), mean(unlist(B_range_Gauss_Jitt[[1]][[1]][[1]]), na.rm=TRUE),
            mean(unlist(B_nuggVar_Gauss_Jitt[[1]][[1]][[1]]), na.rm=TRUE), mean(unlist(B_margVar_Gauss_Jitt[[1]][[1]][[1]]), na.rm=TRUE))

column10 = c(1, "CI Length", mean(unlist(CI_alpha_Gauss_Jitt[[1]][[1]][[1]]), na.rm=TRUE), mean(unlist(CI_range_Gauss_Jitt[[1]][[1]][[1]]), na.rm=TRUE), 
            mean(unlist(CI_nuggVar_Gauss_Jitt[[1]][[1]][[1]]), na.rm=TRUE), mean(unlist(CI_margVar_Gauss_Jitt[[1]][[1]][[1]]), na.rm=TRUE))

column11 = c(4, "Bias", mean(unlist(B_alpha_Gauss_Jitt[[1]][[1]][[2]]), na.rm=TRUE), mean(unlist(B_range_Gauss_Jitt[[1]][[1]][[2]]), na.rm=TRUE), 
            mean(unlist(B_nuggVar_Gauss_Jitt[[1]][[1]][[2]]), na.rm=TRUE), mean(unlist(B_margVar_Gauss_Jitt[[1]][[1]][[2]]), na.rm=TRUE))

column12 = c(4, "CI Length", mean(unlist(CI_alpha_Gauss_Jitt[[1]][[1]][[2]]), na.rm=TRUE), mean(unlist(CI_range_Gauss_Jitt[[1]][[1]][[2]]), na.rm=TRUE), 
            mean(unlist(CI_nuggVar_Gauss_Jitt[[1]][[1]][[2]]), na.rm=TRUE), mean(unlist(CI_margVar_Gauss_Jitt[[1]][[1]][[2]]), na.rm=TRUE))




#Table 3 : Gaussian Standard Model with 340 km of range

column13 = c("Jittering Factor", "Parameters", "$\\alpha$", "$\\tau$", "$\\sigma^{2}_{N}$", "$\\sigma^{2}_{SF}$")
column14 = c("&", "True", 1, 340, 0.1, 1)
column15 = c(1, "Bias", mean(unlist(B_alpha_Gauss_standard[[1]][[2]][[1]]), na.rm=TRUE), mean(unlist(B_range_Gauss_standard[[1]][[2]][[1]]), na.rm=TRUE),
            mean(unlist(B_nuggVar_Gauss_standard[[1]][[2]][[1]]), na.rm=TRUE), mean(unlist(B_margVar_Gauss_standard[[1]][[2]][[1]]), na.rm=TRUE))

column16 = c(1, "CI Length", mean(unlist(CI_alpha_Gauss_standard[[1]][[2]][[1]]), na.rm=TRUE), mean(unlist(CI_range_Gauss_standard[[1]][[2]][[1]]), na.rm=TRUE), 
            mean(unlist(CI_nuggVar_Gauss_standard[[1]][[2]][[1]]), na.rm=TRUE), mean(unlist(CI_margVar_Gauss_standard[[1]][[2]][[1]]), na.rm=TRUE))

column17 = c(4, "Bias", mean(unlist(B_alpha_Gauss_standard[[1]][[2]][[2]]), na.rm=TRUE), mean(unlist(B_range_Gauss_standard[[1]][[2]][[2]]), na.rm=TRUE), 
            mean(unlist(B_nuggVar_Gauss_standard[[1]][[2]][[2]]), na.rm=TRUE), mean(unlist(B_margVar_Gauss_standard[[1]][[2]][[2]]), na.rm=TRUE))

column18 = c(4, "CI Length", mean(unlist(CI_alpha_Gauss_standard[[1]][[2]][[2]]), na.rm=TRUE), mean(unlist(CI_range_Gauss_standard[[1]][[2]][[2]]), na.rm=TRUE), 
            mean(unlist(CI_nuggVar_Gauss_standard[[1]][[2]][[2]]), na.rm=TRUE), mean(unlist(CI_margVar_Gauss_standard[[1]][[2]][[2]]), na.rm=TRUE))


# #Table 4 : Gaussian Jittering Accounted Model with 340 km of range
column19 = c("Jittering Factor", "Parameters", "$\\alpha$", "$\\tau$", "$\\sigma^{2}_{N}$", "$\\sigma^{2}_{SF}$")
column20 = c("&", "True", 1, 340, 0.1, 1)

column21 = c(1, "Bias", mean(unlist(B_alpha_Gauss_Jitt[[1]][[2]][[1]]), na.rm=TRUE), mean(unlist(B_range_Gauss_Jitt[[1]][[2]][[1]]), na.rm=TRUE),
            mean(unlist(B_nuggVar_Gauss_Jitt[[1]][[2]][[1]]), na.rm=TRUE), mean(unlist(B_margVar_Gauss_Jitt[[1]][[2]][[1]]), na.rm=TRUE))

column22 = c(1, "CI Length", mean(unlist(CI_alpha_Gauss_Jitt[[1]][[2]][[1]]), na.rm=TRUE), mean(unlist(CI_range_Gauss_Jitt[[1]][[2]][[1]]), na.rm=TRUE),
            mean(unlist(CI_nuggVar_Gauss_Jitt[[1]][[2]][[1]]), na.rm=TRUE), mean(unlist(CI_margVar_Gauss_Jitt[[1]][[2]][[1]]), na.rm=TRUE))

column23 = c(4, "Bias", mean(unlist(B_alpha_Gauss_Jitt[[1]][[2]][[2]]), na.rm=TRUE), mean(unlist(B_range_Gauss_Jitt[[1]][[2]][[2]]), na.rm=TRUE),
            mean(unlist(B_nuggVar_Gauss_Jitt[[1]][[2]][[2]]), na.rm=TRUE), mean(unlist(B_margVar_Gauss_Jitt[[1]][[2]][[2]]), na.rm=TRUE))

column24 = c(4, "CI Length", mean(unlist(CI_alpha_Gauss_Jitt[[1]][[2]][[2]]), na.rm=TRUE), mean(unlist(CI_range_Gauss_Jitt[[1]][[2]][[2]]), na.rm=TRUE),
            mean(unlist(CI_nuggVar_Gauss_Jitt[[1]][[2]][[2]]), na.rm=TRUE), mean(unlist(CI_margVar_Gauss_Jitt[[1]][[2]][[2]]), na.rm=TRUE))




# Tables for Binomial Simulations

#Table 5 : Standard Model with Binomial simulations and 160 km of range
column25 = c("Jittering Factor", "Parameters", "$\\alpha$", "$\\tau$", "$\\sigma^{2}_{SF}$")
column26 = c("&", "True", 0, 160, 1)
column27 = c(1, "Bias", mean(unlist(B_alpha_Binom_standard[[1]][[1]][[1]]), na.rm=TRUE), mean(unlist(B_range_Binom_standard[[1]][[1]][[1]]), na.rm=TRUE),
            mean(unlist(B_margVar_Binom_standard[[1]][[1]][[1]]), na.rm=TRUE))

column28 = c(1, "CI Length", mean(unlist(CI_alpha_Binom_standard[[1]][[1]][[1]]), na.rm=TRUE), mean(unlist(CI_range_Binom_standard[[1]][[1]][[1]]), na.rm=TRUE), 
            mean(unlist(CI_margVar_Binom_standard[[1]][[1]][[1]]), na.rm=TRUE))

column29 = c(4, "Bias", mean(unlist(B_alpha_Binom_standard[[1]][[1]][[2]]), na.rm=TRUE), mean(unlist(B_range_Binom_standard[[1]][[1]][[2]]), na.rm=TRUE), 
            mean(unlist(B_margVar_Binom_standard[[1]][[1]][[2]]), na.rm=TRUE))

column30 = c(4, "CI Length", mean(unlist(CI_alpha_Binom_standard[[1]][[1]][[2]]), na.rm=TRUE), mean(unlist(CI_range_Binom_standard[[1]][[1]][[2]]), na.rm=TRUE), 
            mean(unlist(CI_margVar_Binom_standard[[1]][[1]][[2]]), na.rm=TRUE))


#Table 6 : Jittering Accounted Model with Binomial simulations and 160 km of range
column31 = c("Jittering Factor", "Parameters", "$\\alpha$", "$\\tau$", "$\\sigma^{2}_{SF}$")
column32 = c("&", "True", 0, 160, 1)

column33 = c(1, "Bias", mean(unlist(B_alpha_Binom_Jitt[[1]][[1]][[1]]), na.rm=TRUE), mean(unlist(B_range_Binom_Jitt[[1]][[1]][[1]]), na.rm=TRUE),
            mean(unlist(B_margVar_Binom_Jitt[[1]][[1]][[1]]), na.rm=TRUE))

column34 = c(1, "CI Length", mean(unlist(CI_alpha_Binom_Jitt[[1]][[1]][[1]]), na.rm=TRUE), mean(unlist(CI_range_Binom_Jitt[[1]][[1]][[1]]), na.rm=TRUE), 
            mean(unlist(CI_margVar_Binom_Jitt[[1]][[1]][[1]]), na.rm=TRUE))

column35 = c(4, "Bias", mean(unlist(B_alpha_Binom_Jitt[[1]][[1]][[2]]), na.rm=TRUE), mean(unlist(B_range_Binom_Jitt[[1]][[1]][[2]]), na.rm=TRUE), 
            mean(unlist(B_margVar_Binom_Jitt[[1]][[1]][[2]]), na.rm=TRUE))

column36 = c(4, "CI Length", mean(unlist(CI_alpha_Binom_Jitt[[1]][[1]][[2]]), na.rm=TRUE), mean(unlist(CI_range_Binom_Jitt[[1]][[1]][[2]]), na.rm=TRUE), 
           mean(unlist(CI_margVar_Binom_Jitt[[1]][[1]][[2]]), na.rm=TRUE))




#Table 7 :Standard Model with Binomial Simulations and 340 km of range

column37 = c("Jittering Factor", "Parameters", "$\\alpha$", "$\\tau$", "$\\sigma^{2}_{SF}$")
column38 = c("&", "True", 0, 340, 1)
column39 = c(1, "Bias", mean(unlist(B_alpha_Binom_standard[[1]][[2]][[1]]), na.rm=TRUE), mean(unlist(B_range_Binom_standard[[1]][[2]][[1]]), na.rm=TRUE),
           mean(unlist(B_margVar_Binom_standard[[1]][[2]][[1]]), na.rm=TRUE))

column40 = c(1, "CI Length", mean(unlist(CI_alpha_Binom_standard[[1]][[2]][[1]]), na.rm=TRUE), mean(unlist(CI_range_Binom_standard[[1]][[2]][[1]]), na.rm=TRUE), 
            mean(unlist(CI_margVar_Binom_standard[[1]][[2]][[1]]), na.rm=TRUE))

column41 = c(4, "Bias", mean(unlist(B_alpha_Binom_standard[[1]][[2]][[2]]), na.rm=TRUE), mean(unlist(B_range_Binom_standard[[1]][[2]][[2]]), na.rm=TRUE), 
           mean(unlist(B_margVar_Binom_standard[[1]][[2]][[2]]), na.rm=TRUE))

column42 = c(4, "CI Length", mean(unlist(CI_alpha_Binom_standard[[1]][[2]][[2]]), na.rm=TRUE), mean(unlist(CI_range_Binom_standard[[1]][[2]][[2]]), na.rm=TRUE), 
            mean(unlist(CI_margVar_Binom_standard[[1]][[2]][[2]]), na.rm=TRUE))


# #Table 8 : Jittering Accounted Model with Binomial Simulations and 340 km of range
column43 = c("Jittering Factor", "Parameters", "$\\alpha$", "$\\tau$", "$\\sigma^{2}_{SF}$")
column44 = c("&", "True", 0, 340, 1)

column45 = c(1, "Bias", mean(unlist(B_alpha_Binom_Jitt[[1]][[2]][[1]]), na.rm=TRUE), mean(unlist(B_range_Binom_Jitt[[1]][[2]][[1]]), na.rm=TRUE),
            mean(unlist(B_margVar_Binom_Jitt[[1]][[2]][[1]]), na.rm=TRUE))

column46 = c(1, "CI Length", mean(unlist(CI_alpha_Binom_Jitt[[1]][[2]][[1]]), na.rm=TRUE), mean(unlist(CI_range_Binom_Jitt[[1]][[2]][[1]]), na.rm=TRUE),
            mean(unlist(CI_margVar_Binom_Jitt[[1]][[2]][[1]]), na.rm=TRUE))

column47 = c(4, "Bias", mean(unlist(B_alpha_Binom_Jitt[[1]][[2]][[2]]), na.rm=TRUE), mean(unlist(B_range_Binom_Jitt[[1]][[2]][[2]]), na.rm=TRUE),
            mean(unlist(B_margVar_Binom_Jitt[[1]][[2]][[2]]), na.rm=TRUE))

column48 = c(4, "CI Length", mean(unlist(CI_alpha_Binom_Jitt[[1]][[2]][[2]]), na.rm=TRUE), mean(unlist(CI_range_Binom_Jitt[[1]][[2]][[2]]), na.rm=TRUE),
           mean(unlist(CI_margVar_Binom_Jitt[[1]][[2]][[2]]), na.rm=TRUE))





data1 = data.frame(column1, column2, column3, column4, column5, column6)
data2 = data.frame(column7, column8, column9, column10, column11, column12)
data3 = data.frame(column13, column14, column15, column16, column17, column18)
data4 = data.frame(column19, column20, column21, column22, column23, column24)
data5 = data.frame(column25, column26, column27, column28, column29, column30)
data6 = data.frame(column31, column32, column33, column34, column35, column36)
data7 = data.frame(column37, column38, column39, column40, column41, column42)
data8 = data.frame(column43, column44, column45, column46, column47, column48)

xtable(data1, caption="Average biases and average CI lengths for Gaussian simulations with range = 160 km (standard model)")
xtable(data2, caption="Average biases and average CI lengths for Gaussian simulations with range = 160 km (jittering accounted model)")
xtable(data3, caption="Average biases and average CI lengths for Gaussian simulations with range = 340 km (standard model)")
xtable(data4, caption="Average biases and average CI lengths for Gaussian simulations with range = 340 km (jittering accounted model)")

xtable(data5, caption="Average biases and average CI lengths for Binomial simulations with range = 160 km (standard model)")
xtable(data6, caption="Average biases and average CI lengths for Binomial simulations with range = 160 km (jittering accounted model)")
xtable(data7, caption="Average biases and average CI lengths for Binomial simulations with range = 340 km (standard model)")
xtable(data8, caption="Average biases and average CI lengths for Binomial simulations with range = 340 km (jittering accounted model)")



#Measuring the runtime

FitTimeStandardBinom = list()
FitTimeJittBinom = list()
FitTimeStandardGauss = list()
FitTimeJittGauss = list()

for (i in 1:nLikelihood){
  FitTimeStandardBinom1 = list()
  FitTimeJittBinom1 = list()
  FitTimeStandardGauss1 = list()
  FitTimeJittGauss1 = list()
  for (j in 1:nRange){
    FitTimeStandardBinom2 = list()
    FitTimeJittBinom2 = list()
    FitTimeStandardGauss2 = list()
    FitTimeJittGauss2 = list()
   
    for (h in 1:nScale){
      
      FitTimeStandardBinom3 = list()
      FitTimeJittBinom3 = list()
      FitTimeStandardGauss3 = list()
      FitTimeJittGauss3 = list()
      
      for (l in 1:nSim){
        if(class(fittedGaussian[[1]][[1]][[i]][[j]][[h]][[l]])  != "try-error"){
          FitTimeStandardGauss3[[l]] = fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["fitStandard_end"]] - fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["fitStandard_start"]]
          FitTimeJittGauss3[[l]] = fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["fitJittAcc_end"]] - fittedGaussian[[1]][[1]][[1]][[j]][[h]][[l]][["fitJittAcc_start"]]
        }else {
          FitTimeStandardGauss3[[l]] = NA
          FitTimeJittGauss3[[l]] = NA
        }
        
        if(class(fittedBinomial[[1]][[1]][[i]][[j]][[h]][[l]])  != "try-error"){
          FitTimeStandardBinom3[[l]] = fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["fitStandard_end"]] - fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["fitStandard_start"]]
          FitTimeJittBinom3[[l]] = fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["fitJittAcc_end"]] - fittedBinomial[[1]][[1]][[1]][[j]][[h]][[l]][["fitJittAcc_start"]]

        } else {
          #Standard model with Binomial Simulations
          FitTimeStandardBinom3[[l]]  = NA
          FitTimeJittBinom3[[l]] = NA
        } 
      }
      
      
      FitTimeStandardBinom2[[h]]  = FitTimeStandardBinom3
      FitTimeJittBinom2[[h]] = FitTimeJittBinom3
      
      FitTimeStandardGauss2[[h]] = FitTimeStandardGauss3
      FitTimeJittGauss2[[h]] = FitTimeJittGauss3
    }
    
    FitTimeStandardBinom1[[j]]  = FitTimeStandardBinom2
    FitTimeJittBinom1[[j]] = FitTimeJittBinom2
    
    FitTimeStandardGauss1[[j]] = FitTimeStandardGauss2
    FitTimeJittGauss1[[j]] = FitTimeJittGauss2

  }
  FitTimeStandardBinom[[i]]  = FitTimeStandardBinom1
  FitTimeJittBinom[[i]] = FitTimeJittBinom1
  
  FitTimeStandardGauss[[i]] = FitTimeStandardGauss1
  FitTimeJittGauss[[i]] = FitTimeJittGauss1
}




tGaussJitt_range160DHS = mean(unlist(FitTimeJittGauss[[1]][[1]][[1]]), na.rm=TRUE)
tGaussJitt_range160Extra = mean(unlist(FitTimeJittGauss[[1]][[1]][[2]]), na.rm=TRUE)
tGaussJitt_range340DHS = mean(unlist(FitTimeJittGauss[[1]][[2]][[1]]), na.rm=TRUE)
tGaussJitt_range340Extra = mean(unlist(FitTimeJittGauss[[1]][[2]][[2]]), na.rm=TRUE)

tBinomJitt_range160DHS = mean(unlist(FitTimeJittBinom[[1]][[1]][[1]]), na.rm=TRUE)
tBinomJitt_range160Extra = mean(unlist(FitTimeJittBinom[[1]][[1]][[2]]), na.rm=TRUE)
tBinomJitt_range340DHS = mean(unlist(FitTimeJittBinom[[1]][[2]][[1]]), na.rm=TRUE)
tBinomJitt_range340Extra =mean(unlist(FitTimeJittBinom[[1]][[2]][[2]]), na.rm=TRUE)


save(tGaussJitt_range160DHS, tGaussJitt_range160Extra, tGaussJitt_range340DHS, tGaussJitt_range340Extra,
     tBinomJitt_range160DHS, tBinomJitt_range160Extra, tBinomJitt_range340DHS, tBinomJitt_range340Extra, file = "~/directory")




library(xtable)
xtable(data1, caption="Average model fitting time for Gaussian simulations")

xtable(data2, caption="Average model fitting time for Binomial simulations")



##############Plots of predictions on the map of Kenya##########################


proj = "+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"

adm1Kenya_trnsfrmd = spTransform(adm1Kenya,proj)
dfKenya <- fortify(adm1Kenya_trnsfrmd, region = "NAME_1")

locObs= data.frame(East= kenya.data$east, North=kenya.data$north)
# for Range = 160, factor = 1 (DHS jittering)

# #Standard  (Mean)
d = data.frame(East=fitted[[1]][[1]][[1]][[1]][[1]][[1]][["predCoords"]][, 1], North=fitted[[1]][[1]][[1]][[1]][[1]][[1]][["predCoords"]][, 2],
               Predictions=fitted[[1]][[1]][[1]][[1]][[1]][[1]][["summStandard"]][,1])

summary(d$Predictions)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00933 0.05729 0.18788 0.28031 0.51207 0.86191     # do this to decide which values to use for the "breaks" argument

ggplot(d, aes(East,North)) + 
  geom_tile(aes(fill=Predictions), width=25,height=25) + #theme_bw() 
  geom_path(data = dfKenya, aes(long,lat, group = group),colour = "white", inherit.aes = FALSE)+
  theme(axis.text.y = element_text(size = 35), axis.text.x = element_text(size = 35)) +
  theme(axis.title.x=element_text(size = rel(2))) + theme(axis.title.y=element_text(size = rel(2)))+
  theme(legend.title = element_text(size = rel(3))) + coord_fixed() + 
  xlab("Easting (km)") + ylab("Northing (km)")  + theme(legend.text=element_text(size=30))+
  scale_fill_viridis_c(option = "viridis", begin = 0.2, end = 1, limits = c(0,0.9), breaks = c(0, 0.5, 0.9)) +geom_point(data = locObs, color = "red", size=1, shape="plus")+
  guides(fill = guide_colourbar(barwidth = 2.5, barheight = 25, title = labs("pred."), title.vjust=3) ) 


ggsave("meanStandard.png", path = "~/directory")

#Jitt accounted

d = data.frame(East=fitted[[1]][[1]][[1]][[1]][[1]][[1]][["predCoords"]][, 1], 
               North=fitted[[1]][[1]][[1]][[1]][[1]][[1]][["predCoords"]][, 2], 
               Predictions=fitted[[1]][[1]][[1]][[1]][[1]][[1]][["summJittAcc"]][,1])

summary(d$Predictions)

# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.009453 0.056821 0.186520 0.279762 0.510756 0.860494       # do this to decide which values to use for the "breaks" argument


ggplot(d, aes(East,North)) + 
  geom_tile(aes(fill=Predictions), width=25,height=25) + #theme_bw() 
  geom_path(data = dfKenya, aes(long,lat, group = group),colour = "white", inherit.aes = FALSE)+
  theme(axis.text.y = element_text(size = 35), axis.text.x = element_text(size = 35)) +
  theme(axis.title.x=element_text(size = rel(2))) + theme(axis.title.y=element_text(size = rel(2)))+
  theme(legend.title = element_text(size = rel(3))) + coord_fixed() + 
  xlab("Easting (km)") + ylab("Northing (km)")  + theme(legend.text=element_text(size=30))+
  scale_fill_viridis_c(option = "viridis", begin = 0.2, end = 1, limits = c(0,0.9), breaks = c(0, 0.5, 0.9)) +geom_point(data = locObs, color = "red", size=1, shape="plus")+
  guides(fill = guide_colourbar(barwidth = 2.5, barheight = 25, title = labs("pred."), title.vjust=3) ) 


ggsave("meanJitt.png", path = "~/directory")

# Coefficient of variation

#Range = 160, factor = 1 (DHS jittering)

# #Standard 
d = data.frame(East=fitted[[1]][[1]][[1]][[1]][[1]][[1]][["predCoords"]][, 1], North=fitted[[1]][[1]][[1]][[1]][[1]][[1]][["predCoords"]][, 2],
               uncertainty=(fitted[[1]][[1]][[1]][[1]][[1]][[1]][["summStandard"]][,3]/fitted[[1]][[1]][[1]][[1]][[1]][[1]][["summStandard"]][,1])*100)

summary(d$uncertainty)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.052  18.453  54.561  57.449  91.451 149.626 

ggplot(d, aes(East,North)) + 
  geom_tile(aes(fill=uncertainty), width=25,height=25) + #theme_bw() 
  geom_path(data = dfKenya, aes(long,lat, group = group),colour = "white", inherit.aes = FALSE)+
  theme(axis.text.y = element_text(size = 35), axis.text.x = element_text(size = 35)) +
  theme(axis.title.x=element_text(size = rel(2))) + theme(axis.title.y=element_text(size = rel(2)))+
  theme(legend.title = element_text(size = rel(3))) + coord_fixed() + 
  xlab("Easting (km)") + ylab("Northing (km)")  + theme(legend.text=element_text(size=30))+
  scale_fill_viridis_c(option = "viridis", begin = 0.2, end = 1, limits = c(3,150), breaks = c(3, 50, 100, 150)) +geom_point(data = locObs, color = "red", size=1, shape="plus")+
  guides(fill = guide_colourbar(barwidth = 2.5, barheight = 25, title = labs("CV (%)"), title.vjust=3) ) 

ggsave("cvStandard.png", path = "~/directory")

#Jitt accounted

d = data.frame(East=fitted[[1]][[1]][[1]][[1]][[1]][[1]][["predCoords"]][, 1], North=fitted[[1]][[1]][[1]][[1]][[1]][[1]][["predCoords"]][, 2],
               uncertainty=(fitted[[1]][[1]][[1]][[1]][[1]][[1]][["summJittAcc"]][,3]/fitted[[1]][[1]][[1]][[1]][[1]][[1]][["summJittAcc"]][,1])*100)
summary(d$uncertainty)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.142  19.076  55.809  58.973  93.859 158.930

ggplot(d, aes(East,North)) + 
  geom_tile(aes(fill=uncertainty), width=25,height=25) + #theme_bw() 
  geom_path(data = dfKenya, aes(long,lat, group = group),colour = "white", inherit.aes = FALSE)+
  theme(axis.text.y = element_text(size = 35), axis.text.x = element_text(size = 35)) +
  theme(axis.title.x=element_text(size = rel(2))) + theme(axis.title.y=element_text(size = rel(2)))+
  theme(legend.title = element_text(size = rel(3))) + coord_fixed() + 
  xlab("Easting (km)") + ylab("Northing (km)")  + theme(legend.text=element_text(size=30))+
  scale_fill_viridis_c(option = "viridis", begin = 0.2, end = 1, limits = c(3,160), breaks = c(3, 50, 100, 160)) +geom_point(data = locObs, color = "red", size=1, shape="plus")+
  guides(fill = guide_colourbar(barwidth = 2.5, barheight = 25, title = labs("CV (%)"), title.vjust=3) ) 


ggsave("cvJitt.png", path = "~/directory")


# Plot observation and prediction locations together on Kenya map 

d = data.frame(East=predCoords[,1], North=predCoords[, 2])

#map the points
ggplot(d, aes(East,North)) + coord_fixed()+
  geom_path(data = dfKenya, aes(long,lat, group = group),colour = "white", inherit.aes = FALSE)+
  theme(axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20)) +
  theme(axis.title.x=element_text(size = rel(2))) + theme(axis.title.y=element_text(size = rel(2)))+
  xlab("Easting (km)") + ylab("Northing (km)")  + theme(legend.text=element_text(size=30))+
  geom_point(data = locObs, color = "red", size=1, shape="plus")+
  geom_point(data = predCoords, color = "blue", size=1, shape="plus")+
  guides(fill = guide_colourbar(barwidth = 2.5, barheight = 25, title = labs("CV (%)"), title.vjust=3) ) 

ggsave("locs.png", path = "~/directory")



# Scatter plots of predictions and coefficient of variation from Model-j versus Model-S

predModelS = fitted[[1]][[1]][[1]][[1]][[1]][[1]][["summStandard"]][,1]
predModelJ = fitted[[1]][[1]][[1]][[1]][[1]][[1]][["summJittAcc"]][,1]

scatterPredictions = data.frame(x = predModelS, y = predModelJ)

ggplot(scatterPredictions, aes(x,y)) +coord_fixed() +#+theme_bw() +
  geom_point(color = "red") +
  geom_abline(slope=1, intercept = 0, color = "black", size=1)+ scale_y_log10() +scale_x_log10() +
  xlab("Model-S") + ylab("Model-J")+
  theme(axis.text.y = element_text(size = 25), axis.text.x = element_text(size = 25)) +
  theme(axis.title.x=element_text(size = rel(3))) + theme(axis.title.y=element_text(size = rel(3)))

ggsave("predictionsCrossPlot.png", path = "~/directory")


uncertaintyModelS=(fitted[[1]][[1]][[1]][[1]][[1]][[1]][["summStandard"]][,3]/fitted[[1]][[1]][[1]][[1]][[1]][[1]][["summStandard"]][,1])*100
uncertaintyModelJ=(fitted[[1]][[1]][[1]][[1]][[1]][[1]][["summJittAcc"]][,3]/fitted[[1]][[1]][[1]][[1]][[1]][[1]][["summJittAcc"]][,1])*100

scatterCV = data.frame(x = uncertaintyModelS, y = uncertaintyModelJ)

ggplot(scatterCV, aes(x,y)) +coord_fixed()+ #theme_bw() +
  geom_point(color = "red") +
  #geom_line(color = "white") +
  geom_abline(slope=1, intercept = 0, color = "black", size =1)+
  xlab("Model-S") + ylab("Model-J")+
  theme(axis.text.y = element_text(size = 25), axis.text.x = element_text(size = 25)) +
  theme(axis.title.x=element_text(size = rel(3))) + theme(axis.title.y=element_text(size = rel(3)))

ggsave("CVCrossPlot.png", path = "~/directory")































