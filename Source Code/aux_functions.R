# Aunxillary function for paper weighted ML approach in supply chain
# Author: Erfanul Hoque
# Date: Jan 2021


######################################################################
# Generate AR(1) model # D_t + \mu = \phi (D_{t-1} - \mu) + a_t
######################################################################
# Simulate AR(p) model

AR_sim <- function(t=100, mean = 0, ar.coef = 0.75, ma.coef = NULL, order = NULL, error = c("norm", "t"), 
                   sd = 0, df=3.5){
  error <- match.arg(error)
  # simulating data from AR(p) model
  D.t <- if (error == "norm") {
    # Note: arima.sim() will assume normally distributed errors with SD = 1
    arima.sim(n = t, list(order = order, ar = ar.coef, ma=ma.coef), sd = sd) + mean
  } else { #long-tailed
    arima.sim(n = t, list(order = order, ar = ar.coef, ma = ma.coef), rand.gen = function(t, ...) rt(t, df = df)) + mean
  }
  return(D.t)
}

# data generation from AR(1)
#set.seed(2021)
#D.t <- AR_sim(t=400, mean = 0, ar.coef = 0.8, order = c(1, 0, 0), error = "norm", sd= 5) # AR(1)
#D.t <- AR_sim(t=400, mean = 500, ar.coef = 0.8, order = c(1, 0, 0), error = "t", df=3.5) # AR(1)
#head(D.t)
# sd(D.t)
# sqrt(5^2/(1-0.8^2))

######################################################################
# Estimate AR(p) model : D_t + \mu = \phi (D_{t-1} - \mu) + a_t
######################################################################
# fit the AR model to simulated data using "arima" function 
# arima: ARIMA Modelling of Time Series

AR_est <- function(data, order = c(1, 0, 0), include.mean = TRUE, # for "arima" function
                   #include.constant = FALSE, # for "Arima" function
                   method = c("CSS-ML", "ML", "CSS")){
  # method: fitting method- maximum likelihood or minimize conditional sum-of-squares. 
  # The default (unless there are missing values) is to use conditional-sum-of-squares 
  # to find starting values, then maximum likelihood.
  #res <- arima(x = data, order = order, include.mean = include.mean)
  res <- forecast::Arima(data, order = order, include.constant = TRUE)
  return(res)
}
# Implementation
#res <- AR_est(data = D.t, order = c(1, 0, 0))
#Mu <- as.numeric(res$coef)[2]
#se.Mu <- as.numeric( sqrt(diag(res$var.coef))[2])
#phi <- as.numeric(res$coef)[1]
#se.phi <- as.numeric( sqrt(diag(res$var.coef))[1])

#######################################################################
# Function to calculate lag L (L= 1, 2, 3, 4, 5,...) values
#######################################################################
lag_data <- function(data, L){
  data <- as.vector(data)
  lag.data <- sapply(1:L, function(l) dplyr::lag(data, l))
  return(lag.data)
}
# Calculate lag values
# lag.data <- lag_data(data = data, L = L)


#######################################################################
# Function to calculate the optimal weights over number of simulations
######################################################################

ML_approach <- function(data, L = 5, replication = 500){
  L <- L
  # Calculate lag values
  lag.data <- lag_data(data = data, L = L)
  # Do the simulation for number of replication times
  # Creating a matrix to store the weights
  all_wts <- matrix(nrow = replication, ncol = ncol(lag.data))
  # Creating an empty vector to store demand forecast
  dem_fore <- matrix(0, nrow = nrow(lag.data) - L, ncol = replication)
  # Creating an empty vector to store Standard deviation
  meanFore <- vector('numeric', length = replication)
  riskFore <- vector('numeric', length = replication)
  # Creating an empty vector to store Sharpe Ratio
  sharpe_ratio <- vector('numeric', length = replication)
  
  error <- matrix(0, nrow = nrow(lag.data) - L, ncol = replication)
  # Creating an empty vector to store results
  SSE <- vector('numeric', length = replication)
  #all_wts <- vector('numeric', length = replication)
  
  # Do the simulation
  for (i in 1:replication) {
    set.seed(2021 + i)
    # Create random weights first.
    wts <- runif(n = L)
    wts <- wts/sum(wts) # standardized the weight to get sum equal 1
    all_wts[i,] <- wts
    # demand forecast (forecast is valid only when t >= L+1)
    demand_fore <- sapply((L+1):nrow(lag.data), function(l) sum(wts * lag.data[l,]))
    # Storing demand forecast values
    dem_fore[,i] <- demand_fore
    # Calculate the risk
    meanForeF <- mean(demand_fore)
    meanFore[i] <- meanForeF 
    dem_risk <- sd(demand_fore)
    riskFore[i] <- dem_risk
    sharpe_ratio[i] <- meanForeF / dem_risk
    
    error[,i] <- (data[-c(1:L)] - dem_fore[,i])
    SSE[i] <- sum(error[,i]^2) 
    
  }
  # Storing the values in the table
  res1 <- data.frame(Weight = all_wts, SSE = SSE)
  res_values <- list(results = res1, demandForecast = dem_fore, error= error, 
                     MeanDemand = meanFore, Risk = riskFore, SRatio = sharpe_ratio,
                     Lagvalues = lag.data) 
  #res1 <- data.frame(Weight = all_wts, MeanDemand = meanFore, Risk = riskFore, SRatio = sharpe_ratio)
  #res_values <- list(results = res1, demandForecast = dem_fore, Lagvalues = lag.data)
  return(res_values)
}

#######################################################################
# Function to calculate the optimal order quantity based on MSE
#######################################################################
ML_approachQt <- function(data, L = 5, replication = 500){
  L <- L
  # Calculate lag values
  lag.data <- lag_data(data = data, L = L)
  # Do the simulation for number of replication times
  # Creating a matrix to store the weights
  all_wts <- matrix(nrow = replication, ncol = ncol(lag.data))
  # Creating an empty vector to store Qt
  dem_fore <- matrix(0, nrow = nrow(lag.data) - L, ncol = replication)
  # Creating an empty vector to store Standard deviation
  meanFore <- vector('numeric', length = replication)
  riskFore <- vector('numeric', length = replication)
  # Creating an empty vector to store Sharpe Ratio
  sharpe_ratio <- vector('numeric', length = replication)
  error <- matrix(0, nrow = nrow(lag.data) - L, ncol = replication)
  # Creating an empty vector to store results
  SSE <- vector('numeric', length = replication)
  
  # Do the simulation
  for (i in 1:replication) {
    set.seed(2020 + i)
    # Create random weights first.
    wts <- runif(n = L)
    wts <- wts/sum(wts) # standardized the weight to get sum equal 1
    all_wts[i,] <- wts
    # demand forecast (forecast is valid only when t >= L+1)
    d.t <- data[length(data)]
    demand_fore <- sapply((L+1):nrow(lag.data), function(l) sum(wts * lag.data[l,]) + d.t)
    # Storing Qt values
    dem_fore[,i] <- demand_fore
    # Calculate the risk
    meanForeF <- mean(demand_fore)
    meanFore[i] <- meanForeF 
    dem_risk <- sd(demand_fore)
    riskFore[i] <- dem_risk
    sharpe_ratio[i] <- meanForeF / dem_risk
    
    error[,i] <- (data[-c(1:L)] - dem_fore[,i])
    SSE[i] <- sum(error[,i]^2) 
    
  }
  # Storing the values in the table
  res1 <- data.frame(Weight = all_wts, SSE = SSE)
  res_values <- list(results = res1, demandForecast = dem_fore, error= error, 
                     MeanDemand = meanFore, Risk = riskFore, SRatio = sharpe_ratio,
                     Lagvalues = lag.data)
  # Storing the values in the table
  #res1 <- data.frame(Weight = all_wts, MeanDemand = meanFore, Risk = riskFore, SRatio = sharpe_ratio)
  #res_values <- list(results = res1, demandForecast = dem_fore, Lagvalues = lag.data)
  return(res_values)
}


#######################################################################################
# Function to calculate the order qantitly and variance of sum of L-step ahead forecasts
#######################################################################################
q.t_minvar <- function(data, phi, Optweights, Optlagvalues, minvar, minMSE, 
                       error.var = NULL, model=c("AR1", "AR1seasonal", "RW")){
  model <- match.arg(model)
  # calculate sum of L-step ahead forecast
  d.t <- data[length(data)]
  d.t2 <- data[length(data)-1]
  L <- length(Optweights)
  if(model == "AR1"){
    Sfore_naive <- (((phi/(1-phi))*(1-phi^L))*d.t )
    Sfore_naiveQt <- (((phi/(1-phi))*(1-phi^L))*d.t ) - (((phi/(1-phi))*(1-phi^L))*d.t2)
    var_naive <- ((phi^2/(1-phi)^2)*(1-phi^L)^2)* var(data) # var(data): sample var of orignal demand
    MSE_naive <- error.var*(1 + (1+phi^2) + (1+phi^2+phi^4) + (1+phi^2+phi^4+phi^6)+
                              (1+phi^2+phi^4+phi^6+phi^8) )        
  }
  else {
    Sfore_naive <- L * d.t
    Sfore_naiveQt <- (L * d.t) - (L * d.t2)
    var_naive <- L^2 * (length(data)*error.var) 
  }
  # calculate sum of L-step ahead forecast  
  Sfore_ML <- sum(Optweights * Optlagvalues)
  # calcualte variance of sum of L-step ahead forecasts
  var_ML <- minvar
  MSE_ML <- minMSE
  
  # Calcualte the order quantity based based on ML approach and using estimated AR parameter
  q.t_naive <- d.t + Sfore_naiveQt
  q.t_naive
  q.t_ML <- d.t + Sfore_ML
  q.t_ML
  
  res <- data.frame(Sfore_naive, Sfore_ML, var_naive, MSE_naive, 
                    var_ML, MSE_ML, q.t_naive, q.t_ML)
  return(res)
}


######################################################################
# Calculate Bullwhip Measure
#####################################################################
# Function to calculate Bullwhip measures

# phi: estimated sutoregressive paramter
# L: lead time

bullwhip <- function(phi, errorvar, L=10, minvar, minMSE, method = c("ML", "TS")) {
  method <- match.arg(method)
  
  var_Dt <- errorvar/(1-phi^2) 
  
  if(method == "TS"){ # AR(1)
    numer <- (1 + phi) * (1 - 2*phi^(L+1)) + (2*phi^(2*(L+1)))
    denom <- (1 - phi)
    ratio <- numer / denom
    var_Qt <- (errorvar/(1-phi^2)) * ratio 
    BW_value1 <- var_Qt / var_Dt # definition 1 from 2015 paper (Chitra)
    BW_value2 <- var_Qt - var_Dt # definition 2 from 2018 paper (Chitra)
  } else
  { # ML method
    var_Qt <- minvar # minvar is the variance of weighted sum
    var_Qt1 <- minMSE # minvar is the variance of weighted sum
    BW_value1 <- var_Qt / var_Dt # definition 1 from 2015 paper (Chitra)
    BW_value2 <- var_Qt - var_Dt # definition 2 from 2018 paper (Chitra)
    BW.value_MSE <- var_Qt1 / var_Dt
  }
  #return(BW_value1)
  return(list(BW = c(BW_value1, BW_value2), VarQt = var_Qt))
}

#Implementation
#source("AR_sim.R")
#source("AR_est.R")
#set.seed(2021)
#D.t <- AR_sim(t=500, mean = 2, ar.coef=0.75, order = c(1, 0, 0), sd = 1, func = "arima.sim") # AR(1)
#res <- AR_est(data = D.t, order = c(1, 0, 0))
#phi <- as.numeric(par$coef)[1]
#errorvar <- res$sigma2
#bullwhip(phi = phi, errorvar =errorvar, L=1)

################################################################
# SMA approach
################################################################
#library(TTR)

sma_method <- function(data, L = 5, method= c("sma","DDsma"), dist = c("norm", "t"), alpha=0.05){
  if(method == "sma"){
    sma.values <- SMA(data, n=L)
    ss <- sma.values
    sd.sma <- rollapply(data, width = L, FUN = sd, by.column = TRUE, fill = NA, align = "right")
  }else{
    # DDSMA
    all_res <- ML_approach(data = data, L = 5, replication = 5000)
    sma.values <- all_res$demandForecast[,which.min(all_res$results$SSE)]
    # sd estimate for DDWMA
    ss <- c(rep(NA,L),sma.values)
    sd.sma <- sapply((L+1):length(data), function(t) sqrt(sum((data[(t-L+1):t] - ss[t])^2) / (L-1)))
  }
  
  #sapply(L:399, function(t) mean(Pt[(t-5+1):t]))
  #Dt <- data[-c(1:L)]
  #sapply(1:394, function(t) sqrt(sum((Dt[t:(t+5-1)] - sma.values[t])^2) / L-1))
  
  # Bolinger band using 2-sigma
  bb.lower.sd <- sma.values - 2*sd.sma
  bb.upper.sd <- sma.values + 2*sd.sma
  
  #calculate residual of demand
  #ss <- c(rep(NA,5), sma.values) 
  res <- na.omit(data - ss)
  # Sign correlation
  rho.cal <- function(y) cor(y-mean(y), sign(y-mean(y)))
  rho <- rho.cal(res)
  
  if(dist == "t"){
    nu.fun <- function (x) rho*(x-1)*beta(x/2,1/2)-2*sqrt(x-2)
    nu <- uniroot(nu.fun, c(2, 15))$root
  }else {nu <- 4}
  
  if(method == "sma"){
    #data driven volatility estimate (DDVE)
    vol.cal <-function(y, rho){return(mean(abs(y-mean(y)))/rho)}
    vol.sma <- rollapply(data, width = L, FUN = vol.cal, rho = rho, by.column = TRUE, fill = NA,
                         align = "right")
  }else{
    # vol.sma  estimate for DDSMA
    #ss <- c(rep(NA,5),sma.values)   
    #vol.sma <- sapply((L+1):length(data), function(t) sum(abs(data[(t-5+1):t]- ss[t]))/rho)
    vol.sma <- sapply((L+1):length(data), function(t) mean(abs(data[(t-5+1):t]- ss[t]))/rho)
  }
  
  if(dist == "norm"){
    bbt.lower.vol1 <- sma.values - qnorm((1 - (alpha/2)), mean = 0, sd = 1)*vol.sma
    bbt.upper.vol1 <- sma.values + qnorm((1 - (alpha/2)), mean = 0, sd = 1)*vol.sma
  } else{
    bbt.lower.vol1 <- sma.values - qstd((1 - (alpha/2)), mean = 0, sd = 1, nu=nu)*vol.sma
    bbt.upper.vol1 <- sma.values + qstd((1 - (alpha/2)), mean = 0, sd = 1, nu=nu)*vol.sma
  }
  
  if (method == "sma"){
    BB.res <- na.omit(as.data.frame(cbind(data[-c(1:L)], sma.values[-c(1:L)], bb.lower.sd[-c(1:L)], bb.upper.sd[-c(1:L)],
                                          bbt.lower.vol1[-c(1:L)], bbt.upper.vol1[-c(1:L)])))
    colnames(BB.res)<-c("Demand", "SMA", "BBL_SD", "BBU_SD","BBL_sigma", "BBU_signam")
  }else{
    BB.res <- na.omit(as.data.frame(cbind(data[-c(1:L)], sma.values, bb.lower.sd, bb.upper.sd,
                                          bbt.lower.vol1, bbt.upper.vol1)))
    colnames(BB.res)<-c("Demand", "DDWMA", "BBL_SD", "BBU_SD","BBL_sigma", "BBU_signam")
  }
  
  return(BB.res)
}  


#######################################################################
# Function to calculate the optimal weights for volatility forecast 
######################################################################
# data: the original data
ML_approach_Vol <- function(data, L = 5, replication = 500){
  L <- L
  # First calculate the optimal weighted demand forecast based on DDWMA
  all_res <- ML_approach(data = data, L = 5, replication = 5000)
  DDWMA.values <- all_res$demandForecast[,which.min(all_res$results$SSE)]
  opt.weight.demand <- all_res$results[which.min(all_res$results$SSE),]
  opt.error.demand <- all_res$error[,which.min(all_res$results$SSE)]
  # Extract the residuals of demand for volatility estimate
  res <- all_res$error[,which.min(all_res$results$SSE)]
  # calculate residual of demand for volatility estimate in alternative way
  # ss <- c(rep(NA,L),DDWMA.values) # first L values are missing
  # res <- na.omit(data - ss)
  
  # Sign correlation and d.f
  rho.cal <- function(y) cor(y-mean(y), sign(y-mean(y)))
  rho <- rho.cal(res)
  
  nu.fun <- function (x) rho*(x-1)*beta(x/2,1/2)-2*sqrt(x-2)
  nu <- uniroot(nu.fun, c(2, 15))$root
  
  # convert residual to standardized using rho to get original volatility
  #res.original <- res/rho 
  res.original.vol <- abs(res)/rho 
  
  # Calculate L-lag values for volatility forecast (L window size)
  lag.data.res <- lag_data(data = res.original.vol, L = L)
  
  # Creating a matrix to store the weights fro volatility
  all_wts <- matrix(nrow = replication, ncol = L)
  # Creating an empty vector to store volatility forecast
  vol_fore <- matrix(0, nrow = nrow(lag.data.res) - L, ncol = replication)
  # Creating an empty vector to store Standard deviation
  mean.vol <- vector('numeric', length = replication)
  sd.vol <- vector('numeric', length = replication)
  # Creating an empty vector to store FESS
  error.vol <- matrix(0, nrow = nrow(lag.data.res) - L, ncol = replication)
  #Creating an empty vector to store SSE
  SSE.vol <- vector('numeric', length = replication)
  
  # Do the simulation to get optimal weight for volatility based on standardized residual data
  for (i in 1:replication) {
    set.seed(2021 + i)
    # Create random weights first.
    wts <- runif(n = L)
    wts <- wts/sum(wts) # standardized the weight to get sum equal 1
    all_wts[i,] <- wts
    # calculate volatility forecast based on L window size (note: first L observations will be missing)
    vol.fore <- sapply((L+1):nrow(lag.data.res), function(t) sum(wts * lag.data.res[t,]))
    # Storing demand forecast values
    vol_fore[,i] <- vol.fore
    
    # Calculate the risk
    meanVolF <- mean(vol.fore)
    mean.vol[i] <- meanVolF 
    vol.risk <- sd(vol.fore)
    sd.vol[i] <- vol.risk
    
    # Calculate one-step ahead forecast error (FESS) to choose optimal weight
    error.vol[,i] <- (res.original.vol[-c(1:L)] - vol_fore[,i])
    SSE.vol[i] <- sum(error.vol[,i]^2) 
  }
  
  # Storing the values in the data table
  res1 <- data.frame(Weight = all_wts, SSE = SSE.vol)
  res_values <- list(results = res1, VolForecast = vol_fore, error.vol = error.vol, 
                     demandForecast=DDWMA.values, opt.weight.demand = opt.weight.demand,
                     error.demand = opt.error.demand,
                     meanVolatility = mean.vol, sdVolatility = sd.vol, est.df=nu, rho=rho) 
  return(res_values)
}

# set.seed(2021)
# data <- AR_sim(t=500, mean = 60, ar.coef = 0.80, order = c(1, 0, 0), error = "t", df=3)
# all_res <- ML_approach_Vol(data = data, L = 5, replication = 500)


################################################################
# Extended Bollinger band using weighted algorithm
################################################################

sma_methodnew <- function(data, L = 5, dist = c("norm", "t"), alpha=0.05){
  
  all_res <- ML_approach_Vol(data = data, L = 5, replication = 5000)
  # Obtain the optimal weight based on minimum FESS
  Optm.res <- all_res$results
  min_vol <- Optm.res[which.min(Optm.res$SSE),]  # subtract the results to calculate the optimal weights
  min_vol
  
  # optimal volatility and weigts
  vol.for <- all_res$VolForecast[,which.min(Optm.res$SSE)]
  opt.weight.vol <- min_vol
  error.vol <- all_res$error.vol[,which.min(Optm.res$SSE)]
  # optimal demand forecast and optimal weight 
  DDWMA.values <- all_res$demandForecast[-c(1:L)]
  opt.weight.demand <- all_res$opt.weight.demand
  error.demand <- all_res$error.demand[-c(1:L)]
  # Sign correlation
  rho <- all_res$rho
  nu <- all_res$est.df
  
  
  if(dist == "norm"){
    bbt.lower.vol1 <- DDWMA.values - qnorm((1 - (alpha/2)), mean = 0, sd = 1)*vol.for
    bbt.upper.vol1 <- DDWMA.values + qnorm((1 - (alpha/2)), mean = 0, sd = 1)*vol.for
  } else{
    bbt.lower.vol1 <- DDWMA.values - qstd((1 - (alpha/2)), mean = 0, sd = 1, nu=nu)*vol.for
    bbt.upper.vol1 <- DDWMA.values + qstd((1 - (alpha/2)), mean = 0, sd = 1, nu=nu)*vol.for
  }
  
  BB.res <- na.omit(as.data.frame(cbind(data[-c(1:(2*L))],DDWMA.values,bbt.lower.vol1, bbt.upper.vol1,vol.for)))
  colnames(BB.res)<-c("Demand", "DDWMA", "BBL_Weightedvol", "BBU_Weightedvol", "weighted_volatility")
  
  opt.weights <- rbind(opt.weight.vol=opt.weight.vol, opt.weight.demand=opt.weight.demand)
  rho.df <- c(rho, df)
  errors <- data.frame(error.demand = error.demand, error.vol=error.vol)
  return(list(BB.res, opt.weights, rho.df, errors))
}  



#######################################################################  
  # Functions for simulation
#######################################################################

  #######################################################################
  # Function to calculate the optimal weights over number of simulations
  ######################################################################
  
  ML_approachTest <- function(data, L = 5, opt.wt, replication=1){
    L <- L
    # Calculate lag values
    lag.data <- lag_data(data = data, L = L)
    # Creating an empty vector to store demand forecast
    dem_fore <- matrix(0, nrow = nrow(lag.data) - L, ncol = replication)
    # Creating an empty vector to store Standard deviation
    meanFore <- vector('numeric', length = replication)
    riskFore <- vector('numeric', length = replication)
    # Creating an empty vector to store Sharpe Ratio
    sharpe_ratio <- vector('numeric', length = replication)
    
    error <- matrix(0, nrow = nrow(lag.data) - L, ncol = replication)
    # Creating an empty vector to store results
    SSE <- vector('numeric', length = replication)
    #all_wts <- vector('numeric', length = replication)
    
    # Do the simulation
    for (i in 1:replication) {
      set.seed(2021 + i)
      wts <- opt.wt
      # demand forecast (forecast is valid only when t >= L+1)
      demand_fore <- sapply((L+1):nrow(lag.data), function(l) sum(wts * lag.data[l,]))
      # Storing demand forecast values
      dem_fore[,i] <- demand_fore
      # Calculate the risk
      meanForeF <- mean(demand_fore)
      meanFore[i] <- meanForeF 
      dem_risk <- sd(demand_fore)
      riskFore[i] <- dem_risk
      sharpe_ratio[i] <- meanForeF / dem_risk
      
      error[,i] <- (data[-c(1:L)] - dem_fore[,i])
      SSE[i] <- sum(error[,i]^2) 
      
    }
    # Storing the values in the table
    res1 <- data.frame(Weight = wts, SSE = SSE)
    res_values <- list(results = res1, demandForecast = dem_fore, error= error, 
                       MeanDemand = meanFore, Risk = riskFore, SRatio = sharpe_ratio,
                       Lagvalues = lag.data) 
    #res1 <- data.frame(Weight = all_wts, MeanDemand = meanFore, Risk = riskFore, SRatio = sharpe_ratio)
    #res_values <- list(results = res1, demandForecast = dem_fore, Lagvalues = lag.data)
    return(res_values)
  }
  
  #######################################################################
  # Function to calculate the optimal order quantity based on MSE
  #######################################################################
  ML_approachQtTest <- function(data, L = 5, opt.wt, replication = 1){
    L <- L
    # Calculate lag values
    lag.data <- lag_data(data = data, L = L)
    # Creating an empty vector to store Qt
    dem_fore <- matrix(0, nrow = nrow(lag.data) - L, ncol = replication)
    # Creating an empty vector to store Standard deviation
    meanFore <- vector('numeric', length = replication)
    riskFore <- vector('numeric', length = replication)
    # Creating an empty vector to store Sharpe Ratio
    sharpe_ratio <- vector('numeric', length = replication)
    error <- matrix(0, nrow = nrow(lag.data) - L, ncol = replication)
    # Creating an empty vector to store results
    SSE <- vector('numeric', length = replication)
    
    # Do the simulation
    for (i in 1:replication) {
      set.seed(2020 + i)
      # Create random weights first.
      wts <- opt.wt
      # demand forecast (forecast is valid only when t >= L+1)
      d.t <- data[length(data)]
      demand_fore <- sapply((L+1):nrow(lag.data), function(l) sum(wts * lag.data[l,]) + d.t)
      # Storing Qt values
      dem_fore[,i] <- demand_fore
      # Calculate the risk
      meanForeF <- mean(demand_fore)
      meanFore[i] <- meanForeF 
      dem_risk <- sd(demand_fore)
      riskFore[i] <- dem_risk
      sharpe_ratio[i] <- meanForeF / dem_risk
      
      error[,i] <- (data[-c(1:L)] - dem_fore[,i])
      SSE[i] <- sum(error[,i]^2) 
      
    }
    # Storing the values in the table
    res1 <- data.frame(Weight = wts, SSE = SSE)
    res_values <- list(results = res1, demandForecast = dem_fore, error= error, 
                       MeanDemand = meanFore, Risk = riskFore, SRatio = sharpe_ratio,
                       Lagvalues = lag.data)
    # Storing the values in the table
    #res1 <- data.frame(Weight = all_wts, MeanDemand = meanFore, Risk = riskFore, SRatio = sharpe_ratio)
    #res_values <- list(results = res1, demandForecast = dem_fore, Lagvalues = lag.data)
    return(res_values)
  }
  
  #####################################################################
  ##  Calculate order quantity and bullwhip effects (L= L+1)
  #####################################################################
  Qt_BWsimTest <- function(data, L = 5, opt.wt, phi.fit, error.var){
    D.t <- data
    all_res <- ML_approachQtTest(data = D.t, L = L+1, opt.wt = opt.wt) # for Q_t calculation
    # Obtain the optimal weight based on minimum variance
    min_sd <-all_res$results
    min_sd
    
    # Estimated AR parameters using training data
    phi.fit <- phi.fit # AR parameter
    error.var <- error.var
    
    # extract the arguments to implement the "q.t_minvar"" function
    L <- L + 1
    Optweights <- opt.wt
    Optweights
    Optlagvalues <- all_res$Lagvalues[length(D.t),]
    Optlagvalues
    minvar <- all_res$Risk^2
    minMSE <- min_sd$SSE/length(D.t)
    
    Final_Qt <- q.t_minvar(data = D.t, phi = phi.fit, Optweights = Optweights, 
                           Optlagvalues = Optlagvalues, minvar = minvar, 
                           error.var = error.var, minMSE=minMSE, model = "AR1")
    
    #Final_Qt
    ## Calculate Bullwhip measure
    BW_ML <- bullwhip(phi = phi.fit, errorvar = error.var, minvar = minvar, 
                      minMSE=minMSE, L= L, method = "ML")
    BW_TS <- bullwhip(phi = phi.fit, errorvar = error.var, minvar = minvar, 
                      minMSE=minMSE, L= L, method = "TS")
    BW_ML; BW_TS
    
    res <- list(Qt = Final_Qt[7:8], BW_ML = BW_ML, BW_TS = BW_TS)
    return(res)
  }
  
  
  #######################################################################
  # Function to calculate the results based on test sample for simulation
  #######################################################################
  
  train_test <- function(data){
    #######################################
    # Training and testing
    ######################################
    # Divide the data fro training and testing
    smp_siz = floor(0.75*length(data))
    train = data[1:smp_siz] 
    test = data[-(1:smp_siz)]
    #tail (train)
    #head (test)
    
    # Implement to training data to get parameters
    
    # Implement DDWMA
    all_res <- ML_approach(data = train , L = 5, replication = 5000)
    
    # Obtain the optimal weight based on minimum variance
    Optm.res <- all_res$results # subtract the results to calculate the optimal weights
    min_sd <- Optm.res[which.min(Optm.res$SSE),]
    min_sd
    
    # Estimate the AR parameters using "arima" function
    par_est <- AR_est(data = train , order = c(1, 0, 0), include.mean = TRUE)
    phi.fit <- as.numeric(par_est$coef)[1] # AR parameter
    error.var <- par_est$sigma2
    
    # Implement to testing data
    
    # AR(1) model
    fit2 <- forecast::Arima(test, model = par_est)
    onestep <- fitted(fit2)
    MSE_ar <- forecast::accuracy(onestep,test)[2]^2
    
    # forecast error
    Onefore_naive <- onestep[length(test)]
    Onefore_naive
    #SSE_ar <- sum(par_est$residuals[-1]^2)
    #MSE_ar <- mean(par_est$residuals[-1]^2)
    #MSE_ar
    
    # DDWMA model
    # extract the arguments to implement the "q.t_minvar"" function
    L <- 5
    Optweights <- min_sd[1:L]
    Optweights
    
    
    all_res <- ML_approachTest(data=test, L = 5, opt.wt=Optweights, replication=1)
    min_sd <- all_res$results
    minMSE <- min_sd$SSE/length(test)
    
    # calculate sum of L-step ahead forecast  
    Sfore_ML <- all_res$demandForecast[length(test)-L,]
    MSE_ML <- minMSE
    
    Simfinal_res <- data.frame(Onefore_naive, Sfore_ML, MSE_ar, MSE_ML) 
    Simfinal_res
    
    # xtable::xtable(Simfinal_res[c(1:2,4,6)], digits=2)
    
    #####################################################################
    ##  Calculate order quantity and bullwhip effects (L= L+1)
    #####################################################################
    
    QtBW <- Qt_BWsimTest(data = test, L=5, opt.wt = Optweights, phi.fit = phi.fit, error.var = error.var)
    QtBW <- data.frame(t(c(unlist(QtBW$Qt), QtBW$BW_TS$BW[1], Bullwhip_ML=QtBW$BW_ML$BW[1])))
    QtBW <- QtBW[c(2,4)]
    QtBW
    #xtable::xtable(QtBW[c(2,4)],digits=2)
    final_res <- list(Simfinal_res, QtBW, all_res)
    return(final_res)
  }
  
  ################################################################
  # SMA approach for test data
  ################################################################
  #library(TTR)
sma_methodTest <- function(data, L = 5, method= c("sma","DDsma"), dist = c("norm", "t"), alpha=0.05){
    
    #######################################
    # Training and testing
    ######################################
    # Divide the data fro training and testing
    smp_siz = floor(0.75*length(data))
    train.data = data[1:smp_siz]  
    test.data = data[-(1:smp_siz)]
    
    if(method == "sma"){
      sma.values <- SMA(test.data, n=L)
      ss <- sma.values
      sd.sma <- rollapply(test.data, width = L, FUN = sd, by.column = TRUE, fill = NA, align = "right")
    }else{
      # DDSMA
      # Implement DDWMA
      all_res1 <- ML_approach(data = train.data , L = 5, replication = 5000)
      # Obtain the optimal weight based on minimum variance
      Optm.res <- all_res1$results # subtract the results to calculate the optimal weights
      min_sd <- Optm.res[which.min(Optm.res$SSE),]
      L <- L
      Optweights <- min_sd[1:L]
      Optweights
      # use optimal weight for test data
      all_res <- ML_approachTest(data = test.data, L = 5, opt.wt = Optweights)
      sma.values <- all_res$demandForecast
      # sd estimate for DDSMA
      ss <- c(rep(NA,L),sma.values)
      sd.sma <- sapply((L+1):length(test.data), function(t) sqrt(sum((test.data[(t-L+1):t] - ss[t])^2) / (L-1)))
    }
    
    # Bolinger band using 2-sigma
    bb.lower.sd <- sma.values - 2*sd.sma
    bb.upper.sd <- sma.values + 2*sd.sma
    
    #calculate residual of demand
    #ss <- c(rep(NA,5), sma.values) 
    res <- na.omit(test.data- ss)
    # Sign correlation
    rho.cal <- function(y) cor(y-mean(y), sign(y-mean(y)))
    rho <- rho.cal(res)
    
    if(dist == "t"){
      nu.fun <- function (x) rho*(x-1)*beta(x/2,1/2)-2*sqrt(x-2)
      nu <- uniroot(nu.fun, c(2, 15))$root
    }else {nu <- 4}
    
    if(method == "sma"){
      #data driven volatility estimate (DDVE)
      vol.cal <-function(y, rho){return(mean(abs(y-mean(y)))/rho)}
      vol.sma <- rollapply(test.data, width = L, FUN = vol.cal, rho = rho, by.column = TRUE, fill = NA,
                           align = "right")
    }else{
      vol.sma <- sapply((L+1):length(test.data), function(t) mean(abs(test.data[(t-5+1):t]- ss[t]))/rho)
    }
    
    if(dist == "norm"){
      bbt.lower.vol1 <- sma.values - qnorm((1 - (alpha/2)), mean = 0, sd = 1)*vol.sma
      bbt.upper.vol1 <- sma.values + qnorm((1 - (alpha/2)), mean = 0, sd = 1)*vol.sma
    } else{
      bbt.lower.vol1 <- sma.values - qstd((1 - (alpha/2)), mean = 0, sd = 1, nu=nu)*vol.sma
      bbt.upper.vol1 <- sma.values + qstd((1 - (alpha/2)), mean = 0, sd = 1, nu=nu)*vol.sma
    }
    
    if (method == "sma"){
      BB.res <- na.omit(as.data.frame(cbind(test.data[-c(1:L)], sma.values[-c(1:L)], bb.lower.sd[-c(1:L)], bb.upper.sd[-c(1:L)],
                                            bbt.lower.vol1[-c(1:L)], bbt.upper.vol1[-c(1:L)])))
      colnames(BB.res)<-c("Demand", "SMA", "BBL_SD", "BBU_SD","BBL_sigma", "BBU_signam")
    }else{
      BB.res <- na.omit(as.data.frame(cbind(test.data[-c(1:L)], sma.values, bb.lower.sd, bb.upper.sd,
                                            bbt.lower.vol1, bbt.upper.vol1)))
      colnames(BB.res)<-c("Demand", "DDSMA", "BBL_SD", "BBU_SD","BBL_sigma", "BBU_signam")
    }
    return(BB.res)
}
 
# Function to get Dynamic RAF plots
sma_methodTestRFS <- function(data, L = 5, method= c("sma","DDsma"), dist = c("norm", "t"), alpha=0.05){
  
  #######################################
  # Training and testing
  ######################################
  # Divide the data fro training and testing
  smp_siz = floor(0.75*length(data))
  train.data = data[1:smp_siz]  
  test.data = data[-(1:smp_siz)]
  
  if(method == "sma"){
    sma.values <- SMA(test.data, n=L)
    ss <- sma.values
    sd.sma <- rollapply(test.data, width = L, FUN = sd, by.column = TRUE, fill = NA, align = "right")
  }else{
    # DDSMA
    # Implement DDWMA
    all_res1 <- ML_approach(data = train.data , L = 5, replication = 5000)
    # Obtain the optimal weight based on minimum variance
    Optm.res <- all_res1$results # subtract the results to calculate the optimal weights
    min_sd <- Optm.res[which.min(Optm.res$SSE),]
    L <- L
    Optweights <- min_sd[1:L]
    Optweights
    # use optimal weight for test data
    all_res <- ML_approachTest(data = test.data, L = 5, opt.wt = Optweights)
    sma.values <- all_res$demandForecast
    # sd estimate for DDSMA
    ss <- c(rep(NA,L),sma.values)
    sd.sma <- sapply((L+1):length(test.data), function(t) sqrt(sum((test.data[(t-L+1):t] - ss[t])^2) / (L-1)))
  }
  
  # Bolinger band using 2-sigma
  bb.lower.sd <- sma.values - 2*sd.sma
  bb.upper.sd <- sma.values + 2*sd.sma
  
  #calculate residual of demand
  #ss <- c(rep(NA,5), sma.values) 
  res <- na.omit(test.data- ss)
  # Sign correlation
  rho.cal <- function(y) cor(y-mean(y), sign(y-mean(y)))
  rho <- rho.cal(res)
  
  if(dist == "t"){
    nu.fun <- function (x) rho*(x-1)*beta(x/2,1/2)-2*sqrt(x-2)
    nu <- uniroot(nu.fun, c(2, 15))$root
  }else {nu <- 4}
  
  if(method == "sma"){
    #data driven volatility estimate (DDVE)
    vol.cal <-function(y, rho){return(mean(abs(y-mean(y)))/rho)}
    vol.sma <- rollapply(test.data, width = L, FUN = vol.cal, rho = rho, by.column = TRUE, fill = NA,
                         align = "right")
  }else{
    vol.sma <- sapply((L+1):length(test.data), function(t) mean(abs(test.data[(t-5+1):t]- ss[t]))/rho)
  }
  
  if(dist == "norm"){
    bbt.lower.vol1 <- sma.values - qnorm((1 - (alpha/2)), mean = 0, sd = 1)*vol.sma
    bbt.upper.vol1 <- sma.values + qnorm((1 - (alpha/2)), mean = 0, sd = 1)*vol.sma
  } else{
    bbt.lower.vol1 <- sma.values - qstd((1 - (alpha/2)), mean = 0, sd = 1, nu=nu)*vol.sma
    bbt.upper.vol1 <- sma.values + qstd((1 - (alpha/2)), mean = 0, sd = 1, nu=nu)*vol.sma
  }
  
  if (method == "sma"){
    BB.res <- na.omit(as.data.frame(cbind(test.data[-c(1:L)], sma.values[-c(1:L)], bb.lower.sd[-c(1:L)], bb.upper.sd[-c(1:L)],
                                          bbt.lower.vol1[-c(1:L)], bbt.upper.vol1[-c(1:L)])))
    colnames(BB.res)<-c("Demand", "SMA", "BBL_SD", "BBU_SD","BBL_sigma", "BBU_signam")
  }else{
    BB.res <- na.omit(as.data.frame(cbind(test.data[-c(1:L)], sma.values, bb.lower.sd, bb.upper.sd,
                                          bbt.lower.vol1, bbt.upper.vol1)))
    colnames(BB.res)<-c("Demand", "DDSMA", "BBL_SD", "BBU_SD","BBL_sigma", "BBU_signam")
  }
  return(list(BB.res, sma.values, vol.sma))
}

 
#######################################################################
# Function to calculate the percentile quantiles
######################################################################

  Qt_percentileS <- function(data, alpha=0.05, dist = c("norm", "t")){
    PT <- sma_methodTest(data = data, L=5, method="sma", dist = dist, alpha=alpha)
    # DDSMA
    PT_ddsma <- sma_methodTest(data = data, L=5, method="DDsma", dist = dist, alpha = alpha)
    # Percentile order quantity
    Qt1 <- PT$BBU_SD[nrow(PT)] + PT$Demand[nrow(PT)]
    Qt2 <- PT_ddsma$BBU_signam[nrow(PT_ddsma)] + PT_ddsma$Demand[nrow(PT_ddsma)]
    return(c(Qt1, Qt2))
  }
  