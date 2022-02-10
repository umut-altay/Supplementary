################################################################################
# Simulation.R                                                                 #
#    This is the main script for the code used to generate results for our     #
#    manuscript.                                                               #
################################################################################



  ## Initialization ----

  
# Libraries
library(SUMMER)
library(rgdal)
library(foreign)
library(INLA)
library(TMB)
library(scoringRules)
library(foreach)
library(geosphere)


# Functions 
source("modSPDEJitter.R")
source("makeIntegrationPoints.R")
source('functions.R')



  ## Setup ----



# Spatial Range (in kilometers) 
rangeSc <- c(160, 340)

# Sub-integration points

nSubRPerPoint <- 10
nSubAPerPoint <- 10

# Jittering Scheme  ( 1--> DHS jittering, 4--> Extra jittering) 
scaleSc <- c(1, 4)

# Likelihood : 1st element  : 0/1 stands for Gaussian/Binomial

#Gaussian
likelihoodSc <- as.matrix(rbind(likelihood = c(0), nuggetVar = c(0.1), p = c(0)))
#Binomial
likelihoodSc <- as.matrix(rbind(likelihood = c(1), nuggetVar = c(0), p = c(0.5)))

# Provincial Boundaries 
boundarySc <- TRUE   #: jittering is done by respecting admin1 borders
boundarySc <- FALSE  #: jittering is done by not respecting admin1 borders

# Number of simulations per scenario
nSim <- 50  


 
## Simulate Data ----



set.seed(2345)

# Geography data 

# # Geography and demography data obtained from SUMMER package

#data(kenyaMaps, package = "SUMMER")             # Maps
#data(kenyaPopulationData, package = "SUMMER")   # Prediction grid

# Geography of survey data
kenya.obsLoc <-readOGR("KEGE71FL.shp")
kenya.obsLoc <- subset(kenya.obsLoc, kenya.obsLoc$LONGNUM != 0)
  
# Responses of survey data
kenya.resp <- read.dta("KEIR72FL.DTA")
  
# Relevant geographical information for clusters
kenya.data <- data.frame(clusterID = kenya.obsLoc@data[["DHSCLUST"]],
                         lon = kenya.obsLoc@coords[,1],
                         lat = kenya.obsLoc@coords[,2], 
                         east = rep(NA, length(kenya.obsLoc@coords[,1])),
                         north = rep(NA, length(kenya.obsLoc@coords[,1])), 
                         urban = kenya.obsLoc@data[["URBAN_RURA"]],
                         adm2 = kenya.obsLoc$ADM1DHS)
  
# Assign UTM37 coordinates to clusters
kenya.data[,c("east", "north")] <- convertDegToKM(kenya.data[,c("lon", "lat")])
  
  
# Extract contraception numbers
contr.data <- getContraception(kenya.resp)
  
# Merge geography and contraception responses
kenya.data <- merge(kenya.data, contr.data, by="clusterID")
  
# NB: row 773 is outside the map in adm2Kenya. Take it away to make tings easier.
kenya.data <- kenya.data[-773,]
 

# Choose a subset of the prediction locations
predCoords <- cbind(popMatKenya$east, popMatKenya$north)
nPred <- 1000
predCoords <- predCoords[1+floor(runif(nPred, max = dim(predCoords)[1])),]


# Create admin1 point in polygon table (object "check1") to check against new locations
true_latLon <- cbind(kenya.data[,"lon"], kenya.data[,"lat"])
true_latLon <- SpatialPoints(true_latLon, proj4string=adm1Kenya@proj4string, bbox = NULL)
check1 <- over(true_latLon, adm1Kenya, returnList = FALSE)
locKM <- kenya.data[,c("east", "north")]


## Simulation 
#Settings
intercept <- 1
sigma.sim <- 1
simLoc <- rbind(as.matrix(kenya.data[, c("east", "north")]), predCoords)
ns <- rep(100, dim(simLoc)[1])


#Simulate the responses and the spatial field

simulatedData <- list()

for (g in 1:length(boundarySc)){
  temp1 <- list()
  for (i in 1:length(likelihoodSc[1,])){
    temp2 <- list()
    for (j in 1:length(rangeSc)){
      temp3 <- list()
      for (h in 1:length(scaleSc)){
        temp4 <- list()
        for(l in 1:nSim){
          
          sigma.cluster.sim <- sqrt(likelihoodSc[[2, i]])
          range.sim <- rangeSc[[j]]
          sim.data <- simulateResponses(loc = simLoc,
                                        ns = ns,
                                        intercept = intercept,
                                        space.range = range.sim,
                                        space.sigma = sigma.sim,
                                        gauss.sim = sigma.cluster.sim)
          
          temp4[[l]] <- sim.data
        }
        temp3[[h]] <- temp4
      }
      temp2[[j]] <- temp3
    }
    temp1[[i]] <- temp2
  }
  simulatedData[[g]] <- temp1
}

#Construct the mesh
mesh.s <- inla.mesh.2d(loc.domain = cbind(kenya.data$east, kenya.data$north),
                       n=5000, 
                       max.n=10000,
                       offset=-.08,
                       cutoff=4, 
                       max.edge=c(25, 50))



# Compile .cpp files

# Model-S
compile( "standard.cpp")
dyn.load( dynlib("standard") )

# Model-J
compile( "jittAccounted.cpp")
dyn.load( dynlib("jittAccounted"))

save.image("tempImage.RData")
load("tempImage.RData")
# 

# set seeds for main and parallel processes
set.seed(123)
totalIter <- length(boundarySc)*length(likelihoodSc[1,])*length(rangeSc)*length(scaleSc)*nSim
# set seeds
allSeeds <- sample(1:1000000, totalIter, replace=FALSE)

#install.packages("parallel")
# set up parallel processes
library(parallel)
cl <- makeCluster(16)
clusterEvalQ(cl, {
  setwd(" ")           #set working directory
  library(SUMMER)
  library(rgdal)
  library(foreign)
  library(INLA)
  library(TMB)
  library(scoringRules)
  
  # Load functions
  source("modSPDEJitter.R")
  source("makeIntegrationPoints.R")
  source('functions.R')
  
  load("tempImage.RData")
  })
clusterExport(cl, c("allSeeds", "totalIter"))

# Prepare inputs for model fitting with TMB

inputs <- list()
locations <- list()
startTime <- proc.time()[3]
for (g in 1:length(boundarySc)){
  tempInput1 <- list()
  tempLoc1 <- list()
  for (i in 1:length(likelihoodSc[1,])){
    tempInput2 <- list()
    tempLoc2 <- list()
    for (j in 1:length(rangeSc)){
      tempInput3 <- list()
      tempLoc3 <- list()
      for (h in 1:length(scaleSc)){
        tempInput4 <- list()
        tempLoc4 <- list()
        
        theseIters <- (g-1)*length(likelihoodSc[1,])*length(rangeSc)*length(scaleSc)*nSim + 
            (i-1) * length(rangeSc)*length(scaleSc)*nSim + 
          (j-1) * length(scaleSc) + 
          (h-1) * nSim + 1:nSim
        clusterExport(cl, "theseIters")
        # set seeds

                # for(l in 1:nSim){
        innerForFun <- function(l) {
          # print(paste0("g: ", g, "/", length(boundarySc), 
          #              ", i: ", i, "/", length(likelihoodSc[1,]), ", ", 
          #              ", j: ", j, "/", length(rangeSc), 
          #              ", h: ", h, "/", length(scaleSc), 
          #              ", l: ", l, "/", nSim))
          
          # set seed
          thisIter = theseIters[l]
          set.seed(allSeeds[thisIter])
          
          boundary = boundarySc[[g]]
          jScale = scaleSc[[h]]
          flag2 = likelihoodSc[[1, i]]
          sigma.cluster.sim = sqrt(likelihoodSc[[2, i]])
          range.sim = rangeSc[[j]]
          
          
          modelParams <- list(intercept = intercept,
                             range.sim = range.sim,  # in kilometers
                             sigma.sim = 1,
                             sigma.cluster.sim = sigma.cluster.sim)
          
          otherValues <- list(USpatial = 1,
                             alphaSpatial = 0.05,
                             flag2 = flag2,
                             jScale = jScale)
          
          
          # TMB Input
          locObs <- Displace(scale = jScale, 
                            locKM = locKM, 
                            urbanRural = kenya.data[,"urban"], 
                            KenyaShapeFile = adm1Kenya, 
                            check1 = check1, 
                            boundary = boundary)
          
          
          locObs <- cbind(locObs[[1]][["east"]], locObs[[1]][["north"]])
          
          if(boundarySc) {
            thisMap <- adm1Kenya
          } else {
            thisMap <- NULL
          }
          
          # data
          tmbInput <- prepare_input(sim.data = simulatedData[[g]][[i]][[j]][[h]][[l]], 
                                   locObs = locObs, 
                                   modelParams = modelParams, 
                                   otherValues = otherValues, 
                                   jScale = jScale, 
                                   urban = kenya.data$urban,
                                   mesh.s = mesh.s, 
                                   adminMap = thisMap, 
                                   nSubRPerPoint=nSubRPerPoint, 
                                   nSubAPerPoint=nSubAPerPoint)
          
          # tempInput4[[l]] = tmbInput
          # tempLoc4[[l]] = locObs
          
          list(tmbInput=tmbInput, locObs=locObs)
        }
        
        clusterExport(cl, c("g", "i", "j", "h"))
        # nOverSim = 60 # simulate more than is
        res = parLapply(cl, 1:nSim, innerForFun)
        
        # extract results
        tempInput4 <- lapply(res, function(x) {x$tmbInput})
        tempLoc4 <- lapply(res, function(x) {x$locObs})
          
        tempInput3[[h]] <- tempInput4
        tempLoc3[[h]] <- tempLoc4
        
        # estimate time remaining
        currTime <- proc.time()[3]
        timeDiff <- currTime - startTime
        numIter <- (g-1)*length(likelihoodSc[1,])*length(rangeSc)*length(scaleSc)*nSim + 
          (i-1) * length(rangeSc)*length(scaleSc)*nSim + 
          (j-1) * length(scaleSc) + 
          (h-1) * nSim + 
          l
        propDone <- numIter/totalIter
        timeRemainingEstimate <- timeDiff * (1 / propDone) - timeDiff
        print(paste0("Estimated time remaining: ", round(timeRemainingEstimate/60/60, digits=2), " hours..."))
      }
      tempInput2[[j]] <- tempInput3
      tempLoc2[[j]] <- tempLoc3
    }
    
    tempInput1[[i]] <- tempInput2
    tempLoc1[[i]] <- tempLoc2
  }
  inputs[[g]] <- tempInput1
  locations[[g]] <- tempLoc1
}
stopCluster(cl)


# parameters
tmb_params <- list(alpha = 0.0, # intercept
                   log_tau = 5, # Log tau (i.e. log spatial precision, Epsilon)
                   log_kappa = -4, # SPDE parameter related to the range
                   Epsilon_s = rep(0, mesh.s[['n']]), # RE on mesh vertices
                   log_nug_std = log(sqrt(0.1))
)

# random effects
rand_effs <- c('Epsilon_s')


# Fit the models with TMB
nLoc <- length(kenya.data$east)


AdministrativeBorders <- list()

for (g in 1:length(boundarySc)){
  Likelihoods <- list()
  for (i in 1:length(likelihoodSc[1,])){
    SpatialRanges <- list()
    for (j in 1:length(rangeSc)){
      JitteringFactors <- list()
      for (h in 1:length(scaleSc)){
        
        data1 <- list()
        data2 <- list()
        u.sim <- list()
        flag2 <- list()
        for (l in 1:nSim){
          data1[[l]] = inputs[[g]][[i]][[j]][[h]][[l]][["data_standard"]]
          data2[[l]] = inputs[[g]][[i]][[j]][[h]][[l]][["data_jittAccounted"]] # data for jittering accounted
          u.sim[[l]] = simulatedData[[g]][[i]][[j]][[h]][[l]][["u.sim"]][-(1:nLoc)]
          flag2[[l]] = likelihoodSc[[1, i]]
        }
        
        Simulations <- list()
        for(l in 1:nSim){ 
          Simulations[[l]] <- try(FitSamplePredict(nLoc =nLoc, 
                                                  intercept = 1, 
                                                  data1 = data1[[l]], 
                                                  data2 = data2[[l]], 
                                                  parameters = tmb_params, 
                                                  random = rand_effs, 
                                                  flag2 = flag2[[l]],
                                                  predCoords = predCoords, 
                                                  mesh.s = mesh.s,
                                                  u.sim = u.sim[[l]]), TRUE)
        }
        
        save(Simulations, g, i, j, h, file = "tempResults.RData")
        
        JitteringFactors[[h]] <- Simulations
      }
      SpatialRanges[[j]] <- JitteringFactors
    }
    
    Likelihoods[[i]] <- SpatialRanges
  }
  
  AdministrativeBorders[[g]] <- Likelihoods
}

fitted <- list(AdministrativeBorders)

save(AdministrativeBorders, fitted, file = "fullResults.RData")

quit()

