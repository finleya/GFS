rm(list=ls())
#
#  Be sure to set working directory!!
#
setwd("//Users//edwingreen//Documents//Bayes text")
v = trees$Volume
d = trees$Girth
h = trees$Height
d2h = d*d*h
Nobs = length(d) # number of observations
library("bayesplot")
library("coda")
library("loo")
k =15000  # posterior sample size
# read posterior samples
# set to directory where posterior samples are stored
#setwd("  ") 
# read in posterior samples (change file name as appropriate)
y = read.table("trees.out",header=FALSE,row.names = NULL)
alpha = y[1:k,2]
beta = y[(k+1):(2*k),2]
rsq = y[((2*k)+1):(3*k),2]
sigma = y[((3*k)+1):(4*k),2]
rm(y)
par(mfrow=c(1,3))
a = density(alpha)
b = density(beta)
c = density(sigma)
plot(a,main=" ",xlab=expression(alpha),
     ylab="density",lwd=1)
plot(b,main=" ",xlab=expression(beta),
     ylab="density",lwd=1)
plot(c,main=" ",xlab=expression(sigma),
     ylab="density",lwd=1)
vrep=matrix(0,nrow=k,ncol=Nobs)
for (i in 1:k){
  mu = alpha[i] + beta[i]*(d2h)    
  vrep[i,] = rnorm(n=Nobs,mean=mu,sd=sigma[i])
}
color_scheme_set("brightblue")
ind = sample(1:15000, 50, replace=FALSE)
ppc_dens_overlay(v, vrep[ind, ])
ind2 = sample(1:15000, 10, replace=FALSE)
ppc_boxplot(v, vrep[ind2, ],notch=FALSE)