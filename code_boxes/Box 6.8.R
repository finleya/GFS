rm(list=ls())
#
#  Be sure to set working directory!!
#
setwd("//Users//edwingreen//Documents//Bayes text")
library(nlme)
library("bayesplot")

head(BodyWeight)
D = BodyWeight[BodyWeight[,4]==2,]
y = D$weight
x = D$Time

Nobs = length(y)
Nrats = 4
Nperrat = Nobs/Nrats
k = 25000  # k = posterior sample size
alpha = matrix(0,nrow=k,ncol=Nrats)
beta = matrix(0,nrow=k,ncol=Nrats)
rsq = matrix(0,nrow=k,ncol=1)
sigma = matrix(0,nrow=k,ncol=1)
alpha_mu = rep(0,k)
alpha_sig = rep(0,k)
beta_mu = rep(0,k)
beta_sig = rep(0,k)

# read posterior samples

z = read.table("rat diet 2.out",header=FALSE,row.names = NULL)
alpha[,1] = z[1:k,2]
alpha[,2] = z[(k+1):(2*k),2]
alpha[,3] = z[((2*k)+1):(3*k),2]
alpha[,4] = z[((3*k)+1):(4*k),2]
alpha_mu[] = z[((4*k)+1):(5*k),2]
alpha_sig[] = z[((5*k)+1):(6*k),2]
beta[,1] = z[((6*k)+1):(7*k),2]
beta[,2] = z[((7*k)+1):(8*k),2]
beta[,3] = z[((8*k)+1):(9*k),2]
beta[,4] = z[((9*k)+1):(10*k),2]
beta_mu[] = z[((10*k)+1):(11*k),2]
beta_sig[]= z[((11*k)+1):(12*k),2]
sigma = z[((12*k)+1):(13*k),2]
rm(z)

yrep1=matrix(0,nrow=k,ncol=Nperrat)
yrep2=matrix(0,nrow=k,ncol=Nperrat)
yrep3=matrix(0,nrow=k,ncol=Nperrat)
yrep4=matrix(0,nrow=k,ncol=Nperrat)

for (i in 1:k){
  mu1 = alpha[i,1] + beta[i,1]*x[1:11]  
  yrep1[i,] = rnorm(n=Nperrat,mean=mu1,sd=sigma[i])
  mu2 = alpha[i,2] + beta[i,2]*x[12:22]  
  yrep2[i,] = rnorm(n=Nperrat,mean=mu2,sd=sigma[i])
  mu3 = alpha[i,3] + beta[i,3]*x[23:33]  
  yrep3[i,] = rnorm(n=Nperrat,mean=mu3,sd=sigma[i])
  mu4 = alpha[i,4] + beta[i,4]*x[34:44]  
  yrep4[i,] = rnorm(n=Nperrat,mean=mu4,sd=sigma[i])
}
  yrep = cbind(yrep1,yrep2,yrep3,yrep4)

color_scheme_set("brightblue")
ind = sample(1:k, 50, replace=FALSE)
ppc_dens_overlay(y, yrep[ind, ])
ind2 = sample(1:k, 10, replace=FALSE)
ppc_boxplot(y, yrep[ind2, ],notch=FALSE)
