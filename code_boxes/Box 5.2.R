rm(list=ls())
#  
# WAIC and LOO for models a - d, doe weight example
#
#  Be sure to set working directory!!
#
setwd(" ")
#
#  The package loo is used to compute LOO and WAIC
#
library("loo")

# input doe weights
y <- c(37, 35, 46, 48, 37, 41, 44, 51,
      41, 48, 41, 43, 43, 31, 35, 42, 
      43, 43, 42, 50, 50, 55, 33, 39)
yr <- c(1,1,1,1,1,1,1,1,
       2,2,2,2,2,2,2,2,
       2,2,2,2,2,2,2,2)
Nobs <- length(y)
# read joint posterior samples

theta_a <- read.table("model a.out",header=FALSE,row.names = NULL)
theta_b <- read.table("model b.out",header=FALSE,row.names = NULL)
theta_c <- read.table("model c.out",header=FALSE,row.names = NULL)
theta_d <- read.table("model d.out",header=FALSE,row.names = NULL)

k <- 10000  # posterior sample size
# intitialize parameter matrices and vectors

mu_a <- matrix(0,nrow=2,ncol=k)
mu_b <- matrix(0,nrow=2,ncol=k)
mu_c <- rep(0,k)
mu_d <- rep(0,k)

sig_a <- matrix(0,nrow=2,ncol=k)
sig_b <- rep(0,k)
sig_c <- matrix(0,nrow=2,ncol=k)
sig_d <- rep(0,k)

# copy posterior samples into appropriate matrices and vectors

mu_a[1,1:k] <- theta_a[1:k,2]
mu_a[2,1:k] <- theta_a[(k+1):(2*k),2]
sig_a[1,1:k] <- theta_a[((2*k)+1):(3*k),2] 
sig_a[2,1:k] <- theta_a[((3*k)+1):(4*k),2] 

mu_b[1,1:k] <- theta_b[1:k,2]
mu_b[2,1:k] <- theta_b[(k+1):(2*k),2]
sig_b[1:k] <- theta_b[((2*k)+1):(3*k),2] 

mu_c[1:k] <- theta_c[1:k,2]
sig_c[1,1:k] <- theta_c[(k+1):(2*k),2] 
sig_c[2,1:k] <- theta_c[((2*k)+1):(3*k),2] 

mu_d[1:k] <- theta_d[1:k,2]
sig_d[1:k] <- theta_d[(k+1):(2*k),2] 

rm(theta_a, theta_b, theta_c, theta_d)

#
# intilalize log-likelihood matrices
#
log_lik_a <- matrix(0,nrow=k,ncol=Nobs)
log_lik_b <- matrix(0,nrow=k,ncol=Nobs)
log_lik_c <- matrix(0,nrow=k,ncol=Nobs)
log_lik_d <- matrix(0,nrow=k,ncol=Nobs)
#
# Compute loglikelihood for each tree and each set of
# parameters in joint poterior sample
#
for (i in 1:k){
  # parameters have index i
  for (j in 1:Nobs){
    # observations have index j
    log_lik_a[i,j] <- dnorm(y[j],mean = mu_a[yr[j],i],sd = sig_a[yr[j],i],log=TRUE)
    log_lik_b[i,j] <- dnorm(y[j],mean = mu_b[yr[j],i],sd = sig_b[i],log=TRUE)
    log_lik_c[i,j] <- dnorm(y[j],mean = mu_c[i],sd = sig_c[yr[j],i],log=TRUE)
    log_lik_d[i,j] <- dnorm(y[j],mean = mu_d[i],sd = sig_d[i],log=TRUE)
  }
}
#
# Get LOO and WAIC from loo
#
LOO_a <- loo(log_lik_a)
WAIC_a <- waic(log_lik_a)
LOO_b <- loo(log_lik_b)
WAIC_b <- waic(log_lik_b)
LOO_c <- loo(log_lik_c)
WAIC_c <- waic(log_lik_c)
LOO_d <- loo(log_lik_d)
WAIC_d <- waic(log_lik_d)

WAIC_a; WAIC_b; WAIC_c; WAIC_d
LOO_a; LOO_b;LOO_c; LOO_d




