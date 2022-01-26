################################################################################
# functions.R                                                                  #
#    These are the helper functions for our main script.                       #
################################################################################

# Convert from lat,lon to UTM37
convertDegToKM = function(loc){
  locLatLon = SpatialPoints(loc, 
                            proj4string = CRS("+proj=longlat +datum=WGS84"))
  locKM = spTransform(locLatLon,
                      CRS("+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"))
  return(locKM@coords[,c(1,2)])
}

convertKMToDeg = function(loc) {
  locSP = SpatialPoints(loc, proj4string=CRS("+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"))
  lonLatCoords = spTransform(locSP, CRS("+proj=longlat +datum=WGS84"))
  attr(lonLatCoords, "coords")
}

################################################################################

# Read contraception
getContraception = function(dat){
  # Convert the answers given to dat$v302a from yes/no into 1s and 0s :
  yesNo = rep(NA, length(dat[,1]))
  for (i in 1:length(dat[,1])){
    if (dat$v302a[i] =="no"){
      yesNo[i] = 0
    } else {
      yesNo[i] = 1
    }
  }
  
  # Associate with clusters
  contraceptUse = data.frame(yesNo = yesNo,
                             clustID = dat$v001)
  
  #sum 1s with respect to cluster IDs
  answers_x = aggregate(contraceptUse$yesNo,
                        by = list(clusterID = contraceptUse$clustID),
                        FUN = sum)
  answers_n= aggregate(contraceptUse$yesNo,
                       by = list(clusterID = contraceptUse$clustID),
                       FUN = length)
  answers_joint = merge(answers_x, answers_n, by="clusterID")
  colnames(answers_joint) = c("clusterID", "ys", "ns")
  
  return(answers_joint)
}

################################################################################

# Matérn covariance function
covMatern = function(dMat, range, stdDev){
  Sig = inla.matern.cov(nu = 1,
                        kappa = sqrt(8*1)/range,
                        x = dMat,
                        corr = TRUE)
  Sig = stdDev^2*Sig
}

################################################################################

# Simulate from model
simulateResponses = function(loc, ns, intercept, space.range, space.sigma, gauss.sim){
  # Spatial covariance matrix
  covMat = covMatern(dMat = as.matrix(dist(loc)),
                     range = space.range,
                     stdDev = space.sigma)
  
  # Add small noise
  covMat = covMat + space.sigma^2*diag(x = 1e-6^2, nrow = dim(covMat)[1], ncol = dim(covMat)[2])
  
  # Simulate spatial effect
  L = t(chol(covMat))
  u.sim = L%*%rnorm(dim(L)[1])
  
  # Add intercept
  lin.pred = intercept + u.sim
  
  # Simulate Gaussian responses
  y.gauss = rnorm(n = length(lin.pred),
                  mean = lin.pred,
                  sd = gauss.sim)
  
  # Simulate Binomial responses
  ps = expit(lin.pred)
  y.binom = rbinom(n = length(lin.pred),
                   size = ns,
                   prob = ps)
  
  return(data.frame(u.sim = u.sim,
                    lin.pred = lin.pred,
                    y.gauss = y.gauss,
                    y.binom = y.binom))
}

################################################################################

# # Simulate jittered coordinates without respecting admin2 areas
# jitterCoord = function(loc, urb, scale){
#   # Number of locations
#   nLoc = dim(loc)[1]
#   
#   # Distance and angle
#   unitDist = runif(n = nLoc)
#   angle = runif(n = nLoc, min = 0, max = 2*pi)
#   
#   # Distance scaling
#   dSc = rep(2, nLoc)
#   dSc[!urb] = 5
#   
#   
#   for (i in 1:length(type)){
#     if (type[[i]]=="U"){
#       distance[[i]]=runif(1, min = 0, max = 2*s)}
#     else {
#       if (runif(1) < 0.01){
#         distance[[i]]=runif(1, min = 0, max = 10*s)
#       } else{
#         distance[[i]]=runif(1, min = 0, max = 5*s)
#       }
#       list(rand.dist=distance)
#     }}
# }

################################################################################
# Generate random distances with respect to the selected jittering scale
# type : a vector of location types : U for urban, R for rural
# s = scaling factor (1, 3 ,5, 10)
random.distance<- function(type, s){      
  distance<- rep(0, length(type))
  for (i in 1:length(type)){
    if (type[[i]]=="U"){
      distance[[i]]=runif(1, min = 0, max = 2*s)}
    else {
      if (runif(1) < 0.01){
        distance[[i]]=runif(1, min = 0, max = 10*s)
      } else{
        distance[[i]]=runif(1, min = 0, max = 5*s)
      }
      list(rand.dist=distance)
    }}
  return(distance)
}

################################################################################
#Generate random displacement angles (between 0 - 2pi) 

#Argument: length: Number of observations.
random.angle<- function(length){
  angle<- (runif(length, min = 0, max = 2*pi))
  return(angle)
}

################################################################################
# Displace locations w.r.t distance and angle generated by two functions above
displace <- function(east, north, angle, distance){
  locx=rep(0, length(north))
  locy=rep(0, length(north))	
  for (i in 1:length(north)){
    (locy[i]=north[i]+((distance[i])*sin(angle[i])))&
      (locx[i]=east[i]+((distance[i])*cos(angle[i])))}
  results=data.frame(locx=locx, locy=locy)
  return(results)
}

################################################################################
# Bring 3 functions above together. Jitter coordinate sets by respecting/not respecting admin1 areas.
# scale-->jittering scale,  locKM--> true locations as east/north (km)
# urbanRural--> should be a vector of U/R, boundary--> TRUE if we respect admin1 boundaries

Displace = function(scale, locKM, urbanRural, KenyaShapeFile, check1, boundary){
  # Jitter each true coordinate one by one and then check the administrative areas they are in now
  #If they landed into a different area then their previous one, jitter that location again
  #Stop when it is jittered and stayed in the same area as befor. Continue jittering with the next location.
  eastOriginal = locKM[, "east"]
  northOriginal = locKM[,"north"]
  nLoc = length(eastOriginal)
  jitteredCoords = list()
  for (j in 1:length(scale)){
    newLocationSet=data.frame(east = rep(NA, nLoc), north = rep(NA, nLoc))
    for (i in 1:nLoc){
      repeat{
        #arguments to be used in jittering:
        east = eastOriginal[i]; north = northOriginal[i]; angle = random.angle(1)
        distance = random.distance(type = urbanRural[i], s = scale[j])
        #jitter location i with respect to scale j  (in east/north format)
        newPoint_eastNorth = displace(east = east, north = north, angle = angle, distance = distance)
        # If we respect admin1 boundaries, check the new point against initial point in polygon table (check1)
        if (boundary == "TRUE"){
          #convert jittered location i to a spatialPoints object
          newPoint_spatialPointsObject = SpatialPoints(cbind(newPoint_eastNorth[,1], newPoint_eastNorth[,2]), proj4string = CRS("+units=km +proj=utm +zone=37 +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +no_defs"), bbox = NULL)
          # Transform it also into longitude/latitude format
          newPoint_longLat <- sp::spTransform(newPoint_spatialPointsObject, KenyaShapeFile@proj4string)
          #see which admin1 area it landed in:
          check2 <- over(newPoint_longLat, KenyaShapeFile, returnList = FALSE)
          #compare it with its previous admin1 area, keep jittering the same point until two areas match
          if ((is.na(check2[,"NAME_1"][[1]]) == FALSE) & (check2[,"NAME_1"][[1]] == check1[,"NAME_1"][[i]])){
            break
          }else{next}
        }else{break}
      } #fill in the jittered location i
      newLocationSet[[i,1]] = newPoint_eastNorth[[1,1]]
      newLocationSet[[i,2]] = newPoint_eastNorth[[1,2]]
    }
    jitteredCoords[[j]] = newLocationSet
  }
  return(jitteredCoords)
}

################################################################################

# prepare TMB input w.r.t Gaussian/Binomial response
# flag2-->0! for Gaussian/Binomial response   
# sim.data--> contains spatial field, simulated Gaussian and Binomial responses

prepare_input = function(sim.data, locObs, modelParams, otherValues, jScale, urban, mesh.s, 
                         adminMap=NULL, nSubAPerPoint=10, nSubRPerPoint=10, testMode=FALSE){
  
  # extract arguments
  USpatial = otherValues[["USpatial"]]
  alphaSpatial = otherValues[["alphaSpatial"]]
  flag2 = otherValues[["flag2"]]
  #flag1 = otherValues[["flag1"]]
  
  range.sim = modelParams[["range.sim"]]
  
  # number of observed locations
  nLoc = length(locObs[,1])
  
  #response variable Gaussian/Binomial
  if (flag2 == 0){
    ys = sim.data$y.gauss[1:nLoc]
    ns = rep(1, nLoc)
  } else {
    ys = sim.data$y.binom[1:nLoc]
    ns = rep(100, nLoc)
  }
  
  
  # spde components
  
  #jittering the points a bit just to make a mesh
  spde = getSPDEPrior(mesh.s, U=USpatial, alpha=alphaSpatial)
  A.proj = inla.spde.make.A(mesh = mesh.s, loc = cbind(locObs[,1], locObs[,2]))
  
  # TMB input for standard (jittering is not accounted for) model
  data_standard = list(num_i = nrow(locObs),  # Total number of observations
                       num_s = mesh.s$n,
                       y_i   =  ys, # num. of pos. obs in the cluster
                       n_i   = ns,  # num. of exposures in the cluster
                       X_alpha  = matrix(1, nrow = nrow(locObs), ncol = 1),# des.mat
                       M0    = spde[['param.inla']][['M0']], # SPDE sparse matrix
                       M1    = spde[['param.inla']][['M1']], # SPDE sparse matrix
                       M2    = spde[['param.inla']][['M2']], # SPDE sparse matrix
                       Aproj = A.proj,             # Projection matrix
                       options = c(1, ## if 1, use normalization trick
                                   1), ## if 1, run adreport
                       flag1 = 1,# normalization flag.
                       flag2 = flag2, #(0/1 for Gaussian/Binomial)
                       alpha_pri = c(0, sqrt(1000)), ## normal
                       matern_pri = c(range.sim, 0.5, USpatial =USpatial , alphaSpatial = alphaSpatial))
  
  # TMB input for the model that accounts for jittering
  
  # convert urban U/R into TRUE/FALSE
  for (i in 1:length(urban)){
    if (urban[i]=='U'){
      urban[i]='TRUE'
    }else{
      urban[i]='FALSE'
    }
  }
  urbanVals=as.logical(urban)
  
  intPointInfo = makeAllIntegrationPoints(coords = cbind(locObs[,1], locObs[,2]), urbanVals, 
                                          numPointsUrban=1+15*4, numPointsRural=1+15*9, 
                                          scalingFactor = jScale, 
                                          JInnerUrban=5, JOuterUrban=0, 
                                          JInnerRural=5, JOuterRural=5, 
                                          adminMap=adminMap, 
                                          nSubAPerPoint=nSubAPerPoint, 
                                          nSubRPerPoint=nSubRPerPoint, 
                                          testMode=testMode)
  if(testMode) {
    return(intPointInfo)
  }
  
  wUrban = intPointInfo$wUrban
  wRural = intPointInfo$wRural
  n_integrationPointsUrban = ncol(wUrban)
  n_integrationPointsRural = ncol(wRural)
  
  
  # Construct projection matrices, and get other relevant info for TMB
  out = makeJitterDataForTMB(intPointInfo, ys , urbanVals, ns, spdeMesh=mesh.s)
  ysUrban = out$ysUrban
  ysRural = out$ysRural
  nsUrban = out$nsUrban
  nsRural = out$nsRural
  AUrban = out$AUrban
  ARural = out$ARural
  
  # Compile inputs for TMB
  data_jittAccounted <- list(num_iUrban = length(ysUrban),  # Total number of urban observations
                             num_iRural = length(ysRural),  # Total number of rural observations
                             num_s = mesh.s[['n']], # num. of vertices in SPDE mesh
                             y_iUrban   = ysUrban, # num. of pos. urban obs in the cluster
                             y_iRural   = ysRural, # num. of pos. rural obs in the cluster
                             n_iUrban   = nsUrban,  # num. of urban exposures in the cluster
                             n_iRural   = nsRural,  # num. of rural exposures in the cluster
                             n_integrationPointsUrban = n_integrationPointsUrban, 
                             n_integrationPointsRural = n_integrationPointsRural, 
                             wUrban = wUrban, 
                             wRural = wRural, 
                             X_alphaUrban  = matrix(1, nrow = length(ysUrban), ncol = 1),# des.mat for urban observations
                             X_alphaRural  = matrix(1, nrow = length(ysRural), ncol = 1),# des.mat for rural observations
                             M0    = spde[['param.inla']][['M0']], # SPDE sparse matrix
                             M1    = spde[['param.inla']][['M1']], # SPDE sparse matrix
                             M2    = spde[['param.inla']][['M2']], # SPDE sparse matrix
                             AprojUrban = AUrban,             # Projection matrix (urban)
                             AprojRural = ARural,             # Projection matrix (rural)
                             options = c(1, ## if 1, use normalization trick
                                         1), ## if 1, run adreport
                             # normalization flag.
                             flag1 = 1,
                             flag2 = flag2, #(0/1 for Gaussian/Binomial)
                             alpha_pri = c(0, sqrt(1000)), ## normal
                             matern_pri = c(range.sim, 0.5, USpatial = USpatial , alphaSpatial = alphaSpatial) 
  )
  
  return(list(data_standard = data_standard, data_jittAccounted = data_jittAccounted, mesh.s = mesh.s))
}

################################################################################

# Model fitting, posterior sampling and prediction with standard model
FitSamplePredict = function(nLoc, intercept, data1, data2, parameters, random, flag2, predCoords, mesh.s, u.sim){
  
  overallSim_start = Sys.time()
  
  map=list(log_nug_std= as.factor(NA))
  
  fitStandard_start = Sys.time()
  # 1.STANDARD MODEL
  # make the autodiff generated likelihood func & gradient
  if (flag2 != 0){
    objSimple <- MakeADFun(data=data1,
                           parameters=parameters,
                           map = map,
                           random=random,
                           hessian=TRUE,
                           DLL='standard')
  } else {
    objSimple <- MakeADFun(data=data1,
                           parameters=parameters,
                           random=random,
                           hessian=TRUE,
                           DLL='standard')
    
  }
  
  # We can normalize the GMRF outside of the nested optimization,
  # avoiding unnecessary and expensive cholesky operations.
  objSimple <- normalize(objSimple, flag="flag1", value = 0)
  
  #newtonOption(objSimple, smartsearch=TRUE)
  
  # Run TMB ----
  opt0 <- nlminb(start  = objSimple[['par']],
                 objective = objSimple[['fn']],
                 gradient = objSimple[['gr']],
                 lower = rep(-10, length(objSimple[['par']])),
                 upper = rep( 10, length(objSimple[['par']])),
                 control = list(trace=0))   
  
  # Get standard errors
  SD0 <- TMB::sdreport(objSimple,
                       par.fixed = opt0$par,
                       getJointPrecision = TRUE,
                       bias.correct = TRUE,
                       bias.correct.control = list(sd = TRUE))
  
  fitStandard_end = Sys.time()
  
  predictionStandard_start = Sys.time()
  # Take samples from fitted model
  mu <- c(SD0$par.fixed,SD0$par.random)
  
  ## Simulate draws
  rmvnorm_prec <- function(mu, chol_prec, n.sims) {
    z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    L <- chol_prec #Cholesky(prec, super=TRUE)
    z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
    z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
    z <- as.matrix(z)
    mu + z
  }
  
  L <- Cholesky(SD0[['jointPrecision']], super = T)
  
  A.pred = inla.spde.make.A(mesh = mesh.s, loc = predCoords)
  
  t.draws <- rmvnorm_prec(mu = mu , chol_prec = L, n.sims = 10000)
  
  # Summarize the draws
  parnames <- c(names(SD0[['par.fixed']]), names(SD0[['par.random']]))
  epsilon_simple_draws  <- t.draws[parnames == 'Epsilon_s',]
  alpha_simple_draws    <- matrix(t.draws[parnames == 'alpha',], nrow = 1)
  
  log_tau_simple_draws  <- t.draws[parnames == 'log_tau',]
  log_kappa_simple_draws  <- t.draws[parnames == 'log_kappa',]
  # 
  sp_range_simple_draws = sqrt(8.0)/exp(log_kappa_simple_draws)
  sp_sigma_simple_draws = 1.0 / sqrt(4.0 * 3.14159265359 *
                                       exp(2.0 * log_tau_simple_draws) * exp(2.0 * log_kappa_simple_draws))
  log_nug_std_simple_draws  <- t.draws[parnames == 'log_nug_std',]
  nugVar_simple_draws = exp(log_nug_std_simple_draws)^2 
  
  
  # Project from mesh to raster, add intercept
  pred_simple <- as.matrix(A.pred %*% epsilon_simple_draws)
  pred_simple <- sweep(pred_simple, 2, alpha_simple_draws, '+')
  
  predictionStandard_end = Sys.time()
  
  #Calculate prediction scores
  fieldSimulated = u.sim
  trueIntercept = rep(intercept, length(fieldSimulated))
  simulatedField_trueIntercept = fieldSimulated + trueIntercept
  
  #Coverage
  nPred = length(pred_simple[,1])
  coverage1 = 0
  for (m in 1:nPred){
    qnt = quantile(pred_simple[m,], c(0.025, 0.975))
    if ((qnt[[1]]<simulatedField_trueIntercept[[m]])&(qnt[[2]]>simulatedField_trueIntercept[[m]])){
      coverage1 = coverage1 + 1
    }
  }
  
  coverage1 = coverage1/nPred
  
  #scores when jittering is not accounted for
  Logscores1 = logs_sample(y = simulatedField_trueIntercept, dat = pred_simple)
  Logscores1 = mean(Logscores1)
  
  CRPSscores1 = crps_sample(y = simulatedField_trueIntercept, dat = pred_simple, method = "edf")
  CRPSscores1 = mean(CRPSscores1)
  
  if (flag2 != 0){
    # Convert to probability scale
    pred_simple = expit(pred_simple)
  }
  # Find the median and sd across draws, as well as 95% intervals
  summStandard <- cbind(mean = (apply(pred_simple, 1, mean)),
                        median = (apply(pred_simple, 1, median)),
                        sd     = (apply(pred_simple, 1, sd)),
                        lower = (apply(pred_simple, 1, quantile, .025)),
                        upper = (apply(pred_simple, 1, quantile, .975)))
  
  summParam_standard = cbind(median = c(range = median(sp_range_simple_draws), margVariance = median(sp_sigma_simple_draws), nugVAR = median(nugVar_simple_draws), alpha = median(alpha_simple_draws)),
                             lower = c(quantile(sp_range_simple_draws, probs = 0.025), quantile(sp_sigma_simple_draws, probs = 0.025), quantile(nugVar_simple_draws, probs = 0.025), quantile(alpha_simple_draws, probs = 0.025)),
                             upper = c(quantile(sp_range_simple_draws, probs = 0.975), quantile(sp_sigma_simple_draws, probs = 0.975), quantile(nugVar_simple_draws, probs = 0.975), quantile(alpha_simple_draws, probs = 0.975)),
                             length = c((quantile(sp_range_simple_draws, probs = 0.975)- c(quantile(sp_range_simple_draws, probs = 0.025))), (quantile(sp_sigma_simple_draws, probs = 0.975)- quantile(sp_sigma_simple_draws, probs = 0.025)), (quantile(nugVar_simple_draws, probs = 0.975) - quantile(nugVar_simple_draws, probs = 0.025)), (quantile(alpha_simple_draws, probs = 0.975) - quantile(alpha_simple_draws, probs = 0.025))))
  # ################################################################################
  # 2.JITTERING IS ACCOUNTED FOR
  
  fitJittAcc_start = Sys.time()
  
  if (flag2 != 0){
    paramsNew <- list(alpha = opt0[["par"]][["alpha"]], # intercept
                      log_tau = opt0[["par"]][["log_tau"]], # Log tau (i.e. log spatial precision, Epsilon)
                      log_kappa = opt0[["par"]][["log_kappa"]], # SPDE parameter related to the range
                      Epsilon_s = rep(0, mesh.s[['n']]), # RE on mesh vertices
                      log_nug_std = log(sqrt(0.1))
    )
    
  } else {
    paramsNew <- list(alpha = opt0[["par"]][["alpha"]], # intercept
                      log_tau = opt0[["par"]][["log_tau"]], # Log tau (i.e. log spatial precision, Epsilon)
                      log_kappa = opt0[["par"]][["log_kappa"]], # SPDE parameter related to the range
                      Epsilon_s = rep(0, mesh.s[['n']]), # RE on mesh vertices
                      log_nug_std = opt0[["par"]][["log_nug_std"]])
    
  }
  
  
  
  ## Make the autodiff generated likelihood func & gradient
  if (flag2 != 0){
    obj <- MakeADFun(data=data2,
                     parameters=paramsNew,
                     map = map,
                     random=random,
                     hessian=TRUE,
                     DLL='jittAccounted')
  } else {
    
    obj <- MakeADFun(data=data2,
                     parameters=paramsNew,
                     random=random,
                     hessian=TRUE,
                     DLL='jittAccounted')
    
  }
  
  ## We can normalize the GMRF outside of the nested optimization,
  ## avoiding unnecessary and expensive cholesky operations.
  obj <- normalize(obj, flag="flag1", value = 0)
  
  #newtonOption(obj, smartsearch=TRUE)
  
  #if (flag2 != 0){
  #save(opt0, data1, data2, i, j, g, l, h, flag2, file = "crashingInputBinomial2ranges2scales.RData") # save the recent values that cause nlminb to crash
  #} else {
  #save(opt0, data1, data2, i, j, g, l, h, flag2, file = "crashingInputGaussian2ranges2scales.RData") # save the recent values that cause nlminb to crash
  #}
  
  # * Run TMB ----
  opt0 <- nlminb(start       =    obj[['par']],
                 objective   =    obj[['fn']],
                 gradient    =    obj[['gr']],
                 lower = rep(-10, length(obj[['par']])),
                 upper = rep( 10, length(obj[['par']])),
                 control     =    list(trace=1))    
  
  # par.fixed = c(0.904220,2.98627,-4.30992,-1.04693)
  # par.fixed2 = c(0.864414,  3.14235, -4.40464, -0.687948)
  # TMB Posterior Sampling ----
  # Get standard errors
  SD0 <- TMB::sdreport(obj,
                       par.fixed = opt0$par,
                       getJointPrecision=TRUE,
                       bias.correct = TRUE,
                       bias.correct.control = list(sd = TRUE))
  
  fitJittAcc_end = Sys.time()
  
  predictionJittAcc_start = Sys.time()
  
  # Take samples from fitted model
  mu <- c(SD0$par.fixed,SD0$par.random)
  
  # Simulate draws
  rmvnorm_prec <- function(mu, chol_prec, n.sims) {
    z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    L <- chol_prec #Cholesky(prec, super=TRUE)
    z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
    z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
    z <- as.matrix(z)
    mu + z
  }
  L <- Cholesky(SD0[['jointPrecision']], super = T)
  
  t.draws <- rmvnorm_prec(mu = mu , chol_prec = L, n.sims = 10000)
  
  # Summarize the draws
  parnames <- c(names(SD0[['par.fixed']]), names(SD0[['par.random']]))
  epsilon_tmb_draws  <- t.draws[parnames == 'Epsilon_s',]
  alpha_tmb_draws    <- matrix(t.draws[parnames == 'alpha',], nrow = 1)
  
  log_tau_tmb_draws  <- t.draws[parnames == 'log_tau',]
  log_kappa_tmb_draws  <- t.draws[parnames == 'log_kappa',]
  
  sp_range_tmb_draws = sqrt(8.0)/exp(log_kappa_tmb_draws)
  sp_sigma_tmb_draws = 1.0 / sqrt(4.0 * 3.14159265359 *
                                    exp(2.0 * log_tau_tmb_draws) * exp(2.0 * log_kappa_tmb_draws))
  log_nug_std_tmb_draws  <- t.draws[parnames == 'log_nug_std',]
  nugVar_tmb_draws = exp(log_nug_std_tmb_draws)^2
  
  # Project from mesh to raster, add intercept
  pred_tmb <- as.matrix(A.pred %*% epsilon_tmb_draws)
  pred_tmb <- sweep(pred_tmb, 2, alpha_tmb_draws, '+')
  
  predictionJittAcc_end = Sys.time()
  
  #Coverage
  nPred = length(pred_tmb[,1])
  coverage2 = 0
  for (m in 1:nPred){
    qnt = quantile(pred_tmb[m,], c(0.025, 0.975))
    if ((qnt[[1]]<simulatedField_trueIntercept[[m]])&(qnt[[2]]>simulatedField_trueIntercept[[m]])){
      coverage2 = coverage2 + 1
    }
  }
  
  coverage2 = coverage2/nPred
  
  # scores when jittering is accounted for
  Logscores2 = logs_sample(y = simulatedField_trueIntercept, dat = pred_tmb)
  Logscores2 = mean(Logscores2)
  
  CRPSscores2 = crps_sample(y = simulatedField_trueIntercept, dat = pred_tmb, method = "edf")
  CRPSscores2 = mean(CRPSscores2)
  
  if (flag2 != 0){
    # Convert to probability scale
    pred_tmb = expit(pred_tmb)
  }
  ## Find the median and sd across draws, as well as 90% intervals
  summJittAcc <- cbind(mean  = (apply(pred_tmb, 1, mean)),
                       median = (apply(pred_tmb, 1, median)),
                       sd     = (apply(pred_tmb, 1, sd)),
                       lower = (apply(pred_tmb, 1, quantile, .025)),
                       upper = (apply(pred_tmb, 1, quantile, .975)))
  
  summParam_JittAcc = cbind(median = c(range = median(sp_range_tmb_draws), margVariance = median(sp_sigma_tmb_draws), nugVAR = median(nugVar_tmb_draws), alpha = median(alpha_tmb_draws)),
                            lower = c(quantile(sp_range_tmb_draws, probs = 0.025), quantile(sp_sigma_tmb_draws, probs = 0.025), quantile(nugVar_tmb_draws, probs = 0.025), quantile(alpha_tmb_draws, probs = 0.025)),
                            upper = c(quantile(sp_range_tmb_draws, probs = 0.975), quantile(sp_sigma_tmb_draws, probs = 0.975), quantile(nugVar_tmb_draws, probs = 0.975), quantile(alpha_tmb_draws, probs = 0.975)),
                            length = c((quantile(sp_range_tmb_draws, probs = 0.975)- c(quantile(sp_range_tmb_draws, probs = 0.025))), (quantile(sp_sigma_tmb_draws, probs = 0.975)- quantile(sp_sigma_tmb_draws, probs = 0.025)), (quantile(nugVar_tmb_draws, probs = 0.975) - quantile(nugVar_tmb_draws, probs = 0.025)), (quantile(alpha_tmb_draws, probs = 0.975) - quantile(alpha_tmb_draws, probs = 0.025))))
  
  
  overallSim_end = Sys.time()
  
  return(list(fitStandard_start = fitStandard_start,
              fitStandard_end = fitStandard_end,
              predictionStandard_start = predictionStandard_start,
              predictionStandard_end = predictionStandard_end,
              fitJittAcc_start = fitJittAcc_start,
              fitJittAcc_end = fitJittAcc_end,
              predictionJittAcc_start = predictionJittAcc_start,
              predictionJittAcc_end = predictionJittAcc_end,
              overallSim_start = overallSim_start,
              overallSim_end = overallSim_end,
              Logscores1 = Logscores1, 
              CRPSscores1 = CRPSscores1,
              coverage1 = coverage1,
              Logscores2 = Logscores2, 
              CRPSscores2 = CRPSscores2,
              coverage2 = coverage2,
              summParam_JittAcc = summParam_JittAcc,
              summParam_standard = summParam_standard,
              summStandard = summStandard,
              summJittAcc = summJittAcc,
              predCoords = predCoords))
}

# functions for calculating point in polygon. Unfortunately, these all seem to be 
# slower than `over'
# # function adapted from SDMTools and 
# # https://stackoverflow.com/questions/49828692/determine-if-a-given-lat-lon-belong-to-a-polygon
# pnt.in.poly2 <- function(pnts, poly.pnts){
#   if (poly.pnts[1, 1] == poly.pnts[nrow(poly.pnts), 1] && 
#       poly.pnts[1, 2] == poly.pnts[nrow(poly.pnts), 2]){ 
#     poly.pnts = poly.pnts[-1, ]
#   }
#   out = pip(pnts[, 1], pnts[, 2], nrow(pnts), poly.pnts[,1], poly.pnts[, 2], nrow(poly.pnts))
#   return(out)
# }

# cppFunction('SEXP pip(SEXP pntx, SEXP pnty, SEXP pntn, SEXP polyx, SEXP polyy, SEXP polyn)
# {
#   // define constants
#   double TWOPI = 2 * PI;
#   double epsilon = 0.000000000001; // threshold value
# 
# 	//define the pointers to the variables
# 	PROTECT(pntx = Rf_coerceVector(pntx, REALSXP)); double *ptx = REAL(pntx); // pnts x values
# 	PROTECT(pnty = Rf_coerceVector(pnty, REALSXP)); double *pty = REAL(pnty); // pnts y values
# 	PROTECT(pntn = Rf_coerceVector(pntn, INTSXP)); int npt = INTEGER(pntn)[0]; // number of points
# 	PROTECT(polyx = Rf_coerceVector(polyx, REALSXP)); double *plx = REAL(polyx); // polygon x values
# 	PROTECT(polyy = Rf_coerceVector(polyy, REALSXP)); double *ply = REAL(polyy); // polygon y values
# 	PROTECT(polyn = Rf_coerceVector(polyn, INTSXP)); int npl = INTEGER(polyn)[0]; // number of polygon points
# 
# 	//define the output variables
# 	SEXP ans; int *out;
# 	PROTECT(ans = Rf_allocVector(INTSXP, npt)); out = INTEGER(ans); //pointer to output dataset
# 
# 	//define some other variables
# 	int ii, jj;
# 	double x, x1, x2, y, y1, y2, dy, dx, dd;
# 
# 	//cycle through the points
# 	for (ii=0;ii<npt;ii++) {
# 		//cycle through the polygon vertices and sum the angles
# 		double angle = 0.0;
# 		for (jj=0;jj<npl;jj++) {
# 			//define the points
# 			x1 = plx[jj]; x2 = plx[(jj+1) % npl]; x = ptx[ii];
# 			y1 = ply[jj]; y2 = ply[(jj+1) % npl]; y = pty[ii];
# 			//check if point are vertix
# 			if (x == x1 && y == y1) { angle = PI+1; break; }
# 			//check if point is on border line between 2 points
# 			if (x == x1 && x == x2) { if ((y1 <= y && y <= y2) || (y1 >= y && y >= y2)) { angle = PI+1; break; } } // check point between two horizontal points
# 			if (y == y1 && y == y2) { if ((x1 <= x && x <= x2) || (x1 >= x && x >= x2)) { angle = PI+1; break; } } // check point between two verticle points
# 			dy = (y1==y2) ? -9999:(y1-y)/(y1-y2); //check if the relative change in x == relative change in y
# 			dx = (x1==x2) ? -9999:(x1-x)/(x1-x2); //check if the relative change in x == relative change in y
# 			dd = dy-dx; dd = (dd<0) ? -dd:dd;
# 			if (dd < epsilon && dy>0 && dy<1) { angle = PI+1; break; } // if dx == dy and dy is between 0 & 1 ... point is on the border line
# 			// && dy > 0 && dy < 1
# 			//if not a vertex or on border lines... sum the angles
# 			double dtheta = atan2(y2 - y, x2 - x) - atan2(y1 - y, x1 - x);
# 			while (dtheta > PI) dtheta -= TWOPI;
# 			while (dtheta < -PI) dtheta += TWOPI;
# 			angle += dtheta;
# 		}
# 		//write out if point is in polygon
# 		if (fabs(angle) < PI) { out[ii] = 0; } else { out[ii] = 1; }
# 	}
# 
# 	//return the output data
# 	UNPROTECT(7);
#     return(ans);
# 
# }')

# function adapted from SDMTools and 
# https://stackoverflow.com/questions/49828692/determine-if-a-given-lat-lon-belong-to-a-polygon
# pnt.in.poly2 <- function(pnts, poly.pnts){
#   if (poly.pnts[1, 1] == poly.pnts[nrow(poly.pnts), 1] && 
#       poly.pnts[1, 2] == poly.pnts[nrow(poly.pnts), 2]){ 
#     poly.pnts = poly.pnts[-1, ]
#   }
#   out = pip(pnts[, 1], pnts[, 2], nrow(pnts), poly.pnts[,1], poly.pnts[, 2], nrow(poly.pnts))
#   return(out)
# }

# function from 
# https://stackoverflow.com/questions/36683825/how-to-check-if-a-point-is-in-a-polygon-effectively-using-r-for-large-data-set
# cppFunction('bool pnpoly(const arma::rowvec& point, const arma::mat& bp) {
#     // Implementation of the ray-casting algorithm is based on
#     // 
#     unsigned int i, j;
# 
#     double x = point(0), y = point(1);
# 
#     bool inside = false;
#     for (i = 0, j = bp.n_rows - 1; i < bp.n_rows; j = i++) {
#       double xi = bp(i,0), yi = bp(i,1);
#       double xj = bp(j,0), yj = bp(j,1);
# 
#       // See if point is inside polygon
#       inside ^= (((yi >= y) != (yj >= y)) && (x <= (xj - xi) * (y - yi) / (yj - yi) + xi));
#     }
# 
#     // Is the cat alive or dead?
#     return inside;
# }', depends="RcppArmadillo")



makeGreenSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  # sequential_hcl(n, h1=260, c1=80, l1=30, l2=90, p1=1.5, rev=TRUE)
  if(!ggplot)
    sequential_hcl(n, h1=128, c1=100, l1=72, l2=95, p1=1.0, rev=TRUE)
  else
    scale_colour_continuous_sequential(h1=128, c1=100, l1=72, l2=95, p1=1.0, rev=TRUE, n_interp=n)
}

# given continuous color scale and range, chooses colors based on a set of values
getColorsFromScale = function(vals, valRange=NULL, cols, scaleFun=function(x) {x}, 
                              forceValuesInRange=FALSE) {
  
  if(is.null(valRange)) {
    nas = !is.finite(scaleFun(vals))
    valRange = range(vals[!nas])
  }
  
  if(forceValuesInRange) {
    vals[vals < valRange[1]] = valRange[1]
    vals[vals > valRange[2]] = valRange[2]
  }
  
  valRange = scaleFun(valRange)
  vals = scaleFun(vals)
  vals = vals - valRange[1]
  vals = vals/(valRange[2] - valRange[1])
  col = cols[round(vals*(length(cols)-1))+1]
  
  col
}

plotWithColor = function(x, y, z, zlim=NULL, colScale=tim.colors(), 
                         legend.mar=7, new=TRUE, scaleFun = function(x) {x}, scaleFunInverse = function(x) {x}, 
                         n.ticks=5, min.n=5, ticks=NULL, tickLabels=NULL, legend.width=1.2, addColorBar=TRUE, 
                         legendArgs=list(), leaveRoomForLegend=TRUE, forceColorsInRange=FALSE, orderI=NULL, 
                         ordering=c("none", "increasing", "decreasing"), colorName = c("col", "bg"), ...) {
  ordering = match.arg(ordering)
  colorName = match.arg(colorName)
  
  # remove NA points
  nas = is.na(x) | is.na(y) | is.na(z)
  if(any(nas)) {
    warning("Removing NAs")
    x = x[!nas]
    y = y[!nas]
    z = z[!nas]
  }
  
  # do setup for ploting data if necessary
  if(is.null(zlim)) {
    nas = !is.finite(scaleFun(z))
    zlim = range(z[!nas])
  }
  
  # order the plotting of the points
  if(is.null(orderI)) {
    if(ordering == "increasing") {
      orderI = sort(z, index.return=TRUE)$ix
    } else if(ordering == "decreasing") {
      orderI = sort(z, decreasing=TRUE, index.return=TRUE)$ix
    } else {
      orderI = 1:length(z)
    }
  }
  x = x[orderI]
  y = y[orderI]
  z = z[orderI]
  
  # if(forceColorsInRange) {
  #   z[z > zlim[2]] = zlim[2]
  #   z[z < zlim[1]] = zlim[1]
  # }
  
  # get colors of points
  cols = getColorsFromScale(z, zlim, colScale, scaleFun, forceColorsInRange)
  
  # generate new plot if necessary
  # browser()
  if(new) {
    # set graphical parameters so the legend won't overlap with plot
    currPar = par()
    newPar = currPar
    newMar = newPar$mar
    if(leaveRoomForLegend) {
      newMar[4] = max(newMar[4], legend.mar)
      newPar$mar = newMar
    }
    if(currPar$mar[4] != newMar[4])
      suppressWarnings({par(newPar)})
    
    # par( oma=c( 0,0,0,6)) # leave room for the legend
    if(colorName == "col") {
      do.call("plot", c(list(x=x, y=y, col=cols), list(...)))
    } else {
      do.call("plot", c(list(x=x, y=y, bg=cols), list(...)))
    }
  } else {
    if(colorName == "col") {
      do.call("points", c(list(x=x, y=y, col=cols), list(...)))
    } else {
      do.call("points", c(list(x=x, y=y, bg=cols), list(...)))
    }
  }
  
  if(addColorBar) {
    # add legend
    # par( oma=c(0,0,0,2))
    if(is.null(ticks))
      ticks = scaleFun(pretty(zlim, n=n.ticks, min.n=min.n))
    else
      ticks = scaleFun(ticks)
    if(is.null(tickLabels))
      tickLabels = scaleFunInverse(ticks)
    
    # par( oma=c( 0,0,0,3))
    
    # set list of arguments to image.plot
    legendArgs$zlim=scaleFun(zlim)
    legendArgs$nlevel=length(colScale)
    legendArgs$legend.only=TRUE
    legendArgs$horizontal=FALSE
    legendArgs$col=colScale
    legendArgs$add = TRUE
    if(is.null(legendArgs$axis.args))
      legendArgs$axis.args=list(at=ticks, labels=tickLabels)
    else {
      legendArgs$axis.args$at=ticks
      legendArgs$axis.args$labels=tickLabels
    }
    legendArgs$legend.mar=legend.mar
    legendArgs$legend.width=legend.width
    do.call("image.plot", legendArgs)
    
    # image.plot(zlim=zlim, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE, 
    #            col=cols, add = TRUE)
  }
  invisible(NULL)
}
