############################################################
####################### Fitting vines ######################
############################################################

# Packages
install.packages("R.matlab")
install.packages("CDvine")
install.packages("VineCopula")
install.packages("network")

# Pad and libraries
rm(list=ls(all=TRUE))
setwd("E:\\Users\\jpvdveld\\Onderzoek\\Data\\2_detrended+vines")
library(CDVine)
library(R.matlab)
library(network)
library(VineCopula)

# Vine fitting TPEpE

names <- c('MPI-rcp45corr_1_TPEpE_') #Change this when fitting different timeseries

files <- character(1) #Empty vector
for (i in 1:1){files[i] <- paste(names[i],'vines.mat',sep="")}

for (k in 1:1){
  
  # Select and load dataset
  name <- names[k]
  file <- files[k]
  data <- list()
  for (i in 1:12){ #For each month
    iname <- paste(name,sprintf('%d',i),sep="") #Makes the name
    data[i] <- readMat(paste(iname,'.mat',sep=""))} #Uses the name to get a file
  
  # Fit vine 
  family <- array(rep(NaN,12*4*4), c(12,4,4))
  pars <- array(rep(NaN,12*4*4), c(12,4,4))
  White1 <- array(rep(NaN, 12*2), c(12,2))
  copMatrix <- array(rep(NaN,12*4*4), c(12,4,4))
  copMatrixbasis <-matrix(data=c(4,3,2,1, 0,3,2,1,0,0,2,1,0,0,0,1), nrow=4,ncol=4)
  for (i in 1:12){ #for each month
    utest <- c()
    xtest <- data[[i]]
    for (j in 1:4){ # 4 variables that are used in the vine copula
      Ftest <- ecdf(xtest[,j])
      tmp <- Ftest(xtest[,j])
      utest <- cbind(utest,tmp)}
    res <- RVineCopSelect(utest, familyset=c(1,3,4,5),copMatrixbasis, selectioncrit="AIC", indeptest=FALSE, rotations=FALSE, level=0.05, method="mle")
    family[i, , ] <- res$family
    pars[i, , ] <- res$par
    copMatrix[i, , ] <- res$Matrix
    Whitestor <- RVineGofTest(utest, res, method="White", B=0)
    White1[i,1] <- Whitestor$White
    White1[i,2] <- Whitestor$p.value
    }
  
  # Save results
  writeMat(file, fams=as.matrix(family), pars=as.matrix(pars), copmatrix=as.matrix(copMatrix))
  
}

# Vine fitting TpPT

names <- c('MPI-rcp45corr_1_TpPT_')

files <- character(1)
for (i in 1:1){files[i] <- paste(names[i],'vines.mat',sep="")}

for (k in 1:1){
  
  # Select and load dataset
  name <- names[k]
  file <- files[k]
  data <- list()
  for (i in 1:12){
    iname <- paste(name,sprintf('%d',i),sep="")
    data[i] <- readMat(paste(iname,'.mat',sep=""))}
  
  # Fit vine 
  family <- array(rep(NaN,12*3*3), c(12,3,3))
  pars <- array(rep(NaN,12*3*3), c(12,3,3))
  copMatrix <- array(rep(NaN,12*3*3), c(12,3,3))
  copMatrixbasis <-matrix(data=c(3,2,1,0,2,1,0,0,1), nrow=3,ncol=3)
  White2 <- array(rep(NaN, 12*2), c(12,2))
  for (i in 1:12){ 
    utest <- c()
    xtest <- data[[i]]
    for (j in 1:3){ # 3 variables that are used in the vine copula 
      Ftest <- ecdf(xtest[,j])
      tmp <- Ftest(xtest[,j])
      utest <- cbind(utest,tmp)}
    res <- RVineCopSelect(utest, familyset=c(1,3,4,5), copMatrixbasis, selectioncrit="AIC", indeptest=FALSE, rotations = FALSE, level=0.05, method="mle")
    family[i, , ] <- res$family
    pars[i, , ] <- res$par
    copMatrix[i, , ] <- res$Matrix
    Whitestor <- RVineGofTest(utest, res, method="White", B=0)
    White2[i,1] <- Whitestor$White
    White2[i,2] <- Whitestor$p.value
    }
  
  # Save results
  writeMat(file, fams=as.matrix(family), pars=as.matrix(pars), copmatrix=as.matrix(copMatrix))
  
}

