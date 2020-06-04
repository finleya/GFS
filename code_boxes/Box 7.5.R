rm(list=ls())
#  Be sure to set directory
setwd(" ")
library("loo")
library("boot")
y = c(0,0,1,1,1,1,0,1,0,1,0,1,0,1,0,1,1,1,0,0,1,1,1,1,1,1,1,1)
Nobs=length(y)
# read posterior samples
k = 25000  # k = posterior sample size
# Initialize matrices to hold BUGS output
alpha = matrix(0,nrow=k,ncol=1)
# read BUGS output
z = read.table("spider_no_covriate.out",header=FALSE,row.names = NULL)
alpha[,1] = z[1:k,2]
# intilalize log-likelihood matrices
log_lik = matrix(0,nrow=k,ncol=Nobs)
# Compute loglikelihood for each tree and each set of
# parameters in joint poterior sample
for (i in 1:k){
  # parameters have index i
  for (j in 1:Nobs){
    # observations have index j
    p = inv.logit(alpha[i])
    log_lik[i,j] = dbinom(y[j],size=1,
                          prob=p,log=TRUE)
  }
}
# Get LOO and WAIC from loo
LOO = loo(log_lik) ; LOO
WAIC = waic(log_lik); WAIC 


