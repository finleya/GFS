rm(list=ls())
#  
# WAIC and LOO for exponential and Weibul models using
#  fire scar data
#  Be sure to set working directory!!
#
setwd(" ")
#
#  The package loo is used to compute LOO and WAIC
#
library("loo")

# fire scar data
y<-c(  
  2, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8, 8, 
  8, 9, 9, 9, 9, 9, 9, 9, 10, 11, 11, 12, 12, 13, 13,
  13, 13, 13, 14, 14, 14, 14, 15, 16, 16, 17, 19, 20,
  21, 24, 25, 25, 30, 30, 31, 31, 31, 31, 31, 31, 33,
  33, 34, 36, 37, 39, 41, 44, 45, 47, 48, 51, 52, 52,
  53, 53, 53, 53, 53, 57, 60, 62, 76, 77, 164)
Nobs <- length(y)
# read joint posterior samples
expon <- read.table("exp.out",header=FALSE,row.names = NULL)
weib <- read.table("Weibull.out",header=FALSE,row.names = NULL)
k <- 10000  # posterior sample size
# intitialize parameter matrices and vectors

lambda <- rep(0,k)
nu <- rep(0,k)
gamma <- rep(0,k)
# copy posterior samples into appropriate matrices and vectors
lambda[1:k] <- expon[1:k,2]
gamma[1:k] <- weib[1:k,2]
nu[1:k] <- weib[((k+1):(2*k)),2] 
rm(expon,weib)
#
# intilalize log-likelihood matrices
#
log_lik_e <- matrix(0,nrow=k,ncol=Nobs)
log_lik_w <- matrix(0,nrow=k,ncol=Nobs)
# Compute loglikelihood for each tree and each set of
# parameters in joint poterior sample
for (i in 1:k){
  # parameters have index i
  for (j in 1:Nobs){
    # observations have index j
    log_lik_e[i,j] <- log(lambda[i])-(y[j]*lambda[i])
    log_lik_w[i,j] <- log(nu[i])+log(gamma[i])+((nu[i]-1)*log(y[j])) - 
      (gamma[i]*(y[j]^nu[i]))
  }
}
#
# Get LOO and WAIC from loo
#
LOO_e <- loo(log_lik_e)
WAIC_e <- waic(log_lik_e)
LOO_w <- loo(log_lik_w)
WAIC_w <- waic(log_lik_w)


LOO_e; LOO_w
WAIC_e; WAIC_w


