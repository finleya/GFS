rm(list=ls())
# Analysis of Beetle data set: WAIC and LOO for hierachical 
# and non-heirarchical models
#
#  Be sure to set working directory!!
setwd(" ")
library("nlme")
library("loo")
library("boot")
#  Data:
rep = c(rep(1,8),rep(2,8))
conc = c(49.06, 52.99, 56.91, 60.84, 64.76, 68.69, 72.61, 76.54)
conc = c(conc,conc)
logcon = log(conc)
corr_conc = logcon - mean(logcon)
pctkill = c( 6.9, 23.3, 32.9, 51.9, 76.7, 93.6,  96.7, 100.0, 
             13.3, 20.0, 26.5, 48.3, 87.9, 85.7, 100.0, 100.0)
total = c(29, 30, 28, 27, 30, 31, 30, 29,
          30, 30, 34, 29, 33, 28, 32, 31)
dead = round(pctkill*total/100)
beetle = cbind(rep,conc,total,dead)
Nobs = length(dead)
#
# read posterior samples
#
k = 25000  # k = posterior sample size

#  Hierarchical model

# Initialize matrices to hold BUGS output

alpha = matrix(0,nrow=k,ncol=2)
beta = matrix(0,nrow=k,ncol=2)

# read BUGS output from hierarchical model

y = read.table("beetle_hier.out",header=FALSE,row.names = NULL)
alpha[,1] = y[1:k,2]
alpha[,2] = y[(k+1):(2*k),2]
beta[,1] = y[((4*k)+1):(5*k),2]
beta[,2] = y[((5*k+1)):(6*k),2]
#
# intilalize log-likelihood matrices
#
log_lik = matrix(0,nrow=k,ncol=Nobs)
#
# Compute loglikelihood for each tree and each set of
# parameters in joint poterior sample
#
for (i in 1:k){
        # parameters have index i
        for (j in 1:Nobs){
                # observations have index j
                p = inv.logit(alpha[i,rep[j]]+
                    beta[i,rep[j]]*corr_conc[j])
                
                log_lik[i,j] = dbinom(dead[j],size=total[j],
                                      prob=p,log=TRUE)
        }
}
#
# Get LOO and WAIC from loo
#
LOOh = loo(log_lik)
WAICh = waic(log_lik)

LOOh
WAICh
##################################################################
# Non-hierarchical model

# Initialize matrices to hold BUGS output

alpha = matrix(0,nrow=k,ncol=1)
beta = matrix(0,nrow=k,ncol=1)

# read BUGS output from non-hierarchical model

y = read.table("beetle_non_hier.out",header=FALSE,row.names = NULL)
alpha[,1] = y[1:k,2]
beta[,1] = y[(k+1):(2*k),2]

#
# intilalize log-likelihood matrices
#
log_lik = matrix(0,nrow=k,ncol=Nobs)
#
# Compute loglikelihood for each tree and each set of
# parameters in joint poterior sample
#
for (i in 1:k){
  # parameters have index i
  for (j in 1:Nobs){
    # observations have index j
    p = inv.logit(alpha[i]+beta[i]*corr_conc[j])
    log_lik[i,j] = dbinom(dead[j],size=total[j],
                          prob=p,log=TRUE)
  }
}
#
# Get LOO and WAIC from loo
#
LOOnh = loo(log_lik)
WAICnh = waic(log_lik)
LOOnh
WAICnh
