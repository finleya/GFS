rm(list=ls())
#  Be sure to set directory
setwd(" ")
library("loo")
library("boot")
x = c(
0.245, 0.247,	0.285,	0.299,	0.327,  0.347,	0.356,	0.36,	0.363,	0.364,	
0.398,	0.4,	0.409,	0.421,	0.432,	0.473,	0.509,	0.529,	0.561,	0.569,	
0.594,	0.638,	0.656,	0.816,	0.853,	0.938,	1.036,1.045)
y = c(0,0,1,1,1,1,0,1,0,1,0,1,0,1,0,1,1,1,0,0,1,1,1,1,1,1,1,1)
xbar = mean(x)
corr_x = x-xbar
Nobs=length(y)
#
# read posterior samples
#
k = 25000  # k = posterior sample size
# Initialize matrices to hold BUGS output
alpha = matrix(0,nrow=k,ncol=1)
beta = matrix(0,nrow=k,ncol=1)
# read BUGS output
z = read.table("spider.out",header=FALSE,row.names = NULL)
alpha[,1] = z[1:k,2]
beta[,1] = z[(k+1):(2*k),2]
# intilalize log-likelihood matrices
log_lik = matrix(0,nrow=k,ncol=Nobs)
# Compute loglikelihood for each tree and each set of
# parameters in joint poterior sample
for (i in 1:k){
  # parameters have index i
  for (j in 1:Nobs){
    # observations have index j
    p = inv.logit(alpha[i]+beta[i]*corr_x[j])
    log_lik[i,j] = dbinom(y[j],size=1,
                          prob=p,log=TRUE)
  }
}
# Get LOO and WAIC from loo
LOO = loo(log_lik) ; LOO
WAIC = waic(log_lik); WAIC 


