rm(list=ls())
#
#  Be sure to set working directory!!
#
setwd("  ")

library(nlme)
library("bayesplot")

head(BodyWeight)
D1 <- BodyWeight[BodyWeight[,4]==1,]
wt1 <- as.matrix(D1$weight)
t1 <- as.matrix(D1$Time)
r1 <- as.matrix(as.numeric(as.character(D1$Rat)))
n1 <- length(wt1)
diet1 <- rep(1,n1)

D2 <- BodyWeight[BodyWeight[,4]==2,]
wt2 <- as.matrix(D2$weight)
t2 <- as.matrix(D2$Time)
r2 <- as.matrix(as.numeric(as.character(D2$Rat)) - 8)
n2 <- length(wt2)
diet2 <- rep(2,n2)

D3 <- BodyWeight[BodyWeight[,4]==3,]
wt3 <- as.matrix(D3$weight)
t3 <- as.matrix(D3$Time)
r3 <- as.matrix(as.numeric(as.character(D3$Rat)) - 12)
n3 <- length(wt3)
diet3 <- rep(3,n3)

wt <- rbind(wt1,wt2,wt3)
wt <- as.vector(wt)
t <- rbind(t1,t2,t3)
r <- rbind(r1,r2,r3)
diet <- c(diet1,diet2,diet3)

rm(D1,D2,D3,diet1,diet2,diet3)

Nobs <- length(wt)
Nrats1 <- 8
Nrats2 <- 4
Nrats3 <- 4
Ndiets <- 3
Nperrat <- 11

k <- 25000  # k = posterior sample size

# Initialize matrices to hold BUGS output

alpha.mean <- matrix(0,nrow=k,ncol=Ndiets)
alpha.sig <- matrix(0,nrow=k,ncol=Ndiets)
alpha1 <- matrix(0,nrow=k,ncol=Nrats1)
alpha2 <- matrix(0,nrow=k,ncol=Nrats2)
alpha3 <- matrix(0,nrow=k,ncol=Nrats3)
beta.mean <- matrix(0,nrow=k,ncol=Ndiets)
beta.sig <- matrix(0,nrow=k,ncol=Ndiets)
beta1 <- matrix(0,nrow=k,ncol=Nrats1)
beta2 <- matrix(0,nrow=k,ncol=Nrats2)
beta3 <- matrix(0,nrow=k,ncol=Nrats3)

mu.alpha <- rep(0,k)
mu.beta <- rep(0,k)
sig.alpha <- rep(0,k)
sig.beta <- rep(0,k)
sigma <- matrix(0,nrow=k,ncol=Ndiets)

# read BUGS output

y <- read.table("full rat growth model.out",header=FALSE,row.names = NULL)
alpha.mean[,1] <- y[1:k,2]
alpha.mean[,2] <- y[(k+1):(2*k),2]
alpha.mean[,3] <- y[((2*k)+1):(3*k),2]
alpha.sig[,1] <- y[((3*k)+1):(4*k),2]
alpha.sig[,2] <- y[((4*k)+1):(5*k),2]
alpha.sig[,3] <- y[((5*k)+1):(6*k),2]
alpha1[,1] <- y[((6*k)+1):(7*k),2]
alpha1[,2] <- y[((7*k)+1):(8*k),2]
alpha1[,3] <- y[((8*k)+1):(9*k),2]
alpha1[,4] <- y[((9*k)+1):(10*k),2]
alpha1[,5] <- y[((10*k)+1):(11*k),2]
alpha1[,6] <- y[((11*k)+1):(12*k),2]
alpha1[,7] <- y[((12*k)+1):(13*k),2]
alpha1[,8] <- y[((13*k)+1):(14*k),2]
alpha2[,1] <- y[((14*k)+1):(15*k),2]
alpha2[,2] <- y[((15*k)+1):(16*k),2]
alpha2[,3] <- y[((16*k)+1):(17*k),2]
alpha2[,4] <- y[((17*k)+1):(18*k),2]
alpha3[,1] <- y[((18*k)+1):(19*k),2]
alpha3[,2] <- y[((19*k)+1):(20*k),2]
alpha3[,3] <- y[((20*k)+1):(21*k),2]
alpha3[,4] <- y[((21*k)+1):(22*k),2]

alpha <- array(0,dim=c(25000,8,3))

for (j in 1:k){
  alpha[j,1,1] <- alpha1[j,1] 
  alpha[j,2,1] <- alpha1[j,2]
  alpha[j,3,1] <- alpha1[j,3] 
  alpha[j,4,1] <- alpha1[j,4]
  alpha[j,5,1] <- alpha1[j,5] 
  alpha[j,6,1] <- alpha1[j,6]
  alpha[j,7,1] <- alpha1[j,7] 
  alpha[j,8,1] <- alpha1[j,8]
  
  alpha[j,1,2] <- alpha2[j,1] 
  alpha[j,2,2] <- alpha2[j,2]
  alpha[j,3,2] <- alpha2[j,3] 
  alpha[j,4,2] <- alpha2[j,4]
  
  alpha[j,1,3] <- alpha3[j,1] 
  alpha[j,2,3] <- alpha3[j,2]
  alpha[j,3,3] <- alpha3[j,3] 
  alpha[j,4,3] <- alpha3[j,4]
}

beta.mean[,1] <- y[((22*k)+1):(23*k),2]
beta.mean[,2] <- y[((23*k)+1):(24*k),2]
beta.mean[,3] <- y[((24*k)+1):(25*k),2]
beta.sig[,1] <- y[((25*k)+1):(26*k),2]
beta.sig[,2] <- y[((26*k)+1):(27*k),2]
beta.sig[,3] <- y[((27*k)+1):(28*k),2]
beta1[,1] <- y[((28*k)+1):(29*k),2]
beta1[,2] <- y[((29*k)+1):(30*k),2]
beta1[,3] <- y[((30*k)+1):(31*k),2]
beta1[,4] <- y[((31*k)+1):(32*k),2]
beta1[,5] <- y[((32*k)+1):(33*k),2]
beta1[,6] <- y[((33*k)+1):(34*k),2]
beta1[,7] <- y[((34*k)+1):(35*k),2]
beta1[,8] <- y[((35*k)+1):(36*k),2]
beta2[,1] <- y[((36*k)+1):(37*k),2]
beta2[,2] <- y[((37*k)+1):(38*k),2]
beta2[,3] <- y[((38*k)+1):(39*k),2]
beta2[,4] <- y[((39*k)+1):(40*k),2]
beta3[,1] <- y[((40*k)+1):(41*k),2]
beta3[,2] <- y[((41*k)+1):(42*k),2]
beta3[,3] <- y[((42*k)+1):(43*k),2]
beta3[,4] <- y[((43*k)+1):(44*k),2]

beta <- array(0,dim=c(25000,8,4))

for (j in 1:k){
  beta[j,1,1] <- beta1[j,1] 
  beta[j,2,1] <- beta1[j,2]
  beta[j,3,1] <- beta1[j,3] 
  beta[j,4,1] <- beta1[j,4]
  beta[j,5,1] <- beta1[j,5] 
  beta[j,6,1] <- beta1[j,6]
  beta[j,7,1] <- beta1[j,7] 
  beta[j,8,1] <- beta1[j,8]
  
  beta[j,1,2] <- beta2[j,1] 
  beta[j,2,2] <- beta2[j,2]
  beta[j,3,2] <- beta2[j,3] 
  beta[j,4,2] <- beta2[j,4]
  
  beta[j,1,3] <- beta3[j,1] 
  beta[j,2,3] <- beta3[j,2]
  beta[j,3,3] <- beta3[j,3] 
  beta[j,4,3] <- beta3[j,4]
}

mu.alpha[] <- y[((50*k)+1):(51*k),2]
mu.beta[] <- y[((51*k)+1):(52*k),2]

sig.alpha[]  <- y[((52*k)+1):(53*k),2]
sig.beta[] <- y[((53*k)+1):(54*k),2]

sigma[,1] <- y[((54*k)+1):(55*k),2]
sigma[,2] <- y[((55*k)+1):(56*k),2]
sigma[,3] <- y[((56*k)+1):(57*k),2]
rm(y)

# read posterior samples

wrep1_1 <- matrix(0,nrow=k,ncol=Nperrat)
wrep1_2 <- matrix(0,nrow=k,ncol=Nperrat)
wrep1_3 <- matrix(0,nrow=k,ncol=Nperrat)
wrep1_4 <- matrix(0,nrow=k,ncol=Nperrat)
wrep1_5 <- matrix(0,nrow=k,ncol=Nperrat)
wrep1_6 <- matrix(0,nrow=k,ncol=Nperrat)
wrep1_7 <- matrix(0,nrow=k,ncol=Nperrat)
wrep1_8 <- matrix(0,nrow=k,ncol=Nperrat)

wrep2_1 <- matrix(0,nrow=k,ncol=Nperrat)
wrep2_2 <- matrix(0,nrow=k,ncol=Nperrat)
wrep2_3 <- matrix(0,nrow=k,ncol=Nperrat)
wrep2_4 <- matrix(0,nrow=k,ncol=Nperrat)

wrep3_1 <- matrix(0,nrow=k,ncol=Nperrat)
wrep3_2 <- matrix(0,nrow=k,ncol=Nperrat)
wrep3_3 <- matrix(0,nrow=k,ncol=Nperrat)
wrep3_4 <- matrix(0,nrow=k,ncol=Nperrat)

for (i in 1:k){
  mu1_1 <- alpha1[i,1] + beta1[i,1]*t[1:11]  
  wrep1_1[i,] <- rnorm(n=Nperrat,mean=mu1_1,sd=sigma[i,1])
  mu1_2 <- alpha1[i,2] + beta1[i,2]*t[12:22]  
  wrep1_2[i,] <- rnorm(n=Nperrat,mean=mu1_2,sd=sigma[i,1])
  mu1_3 <- alpha1[i,3] + beta1[i,3]*t[23:33]  
  wrep1_3[i,] <- rnorm(n=Nperrat,mean=mu1_3,sd=sigma[i,1])
  mu1_4 <- alpha1[i,4] + beta1[i,4]*t[34:44]  
  wrep1_4[i,] <- rnorm(n=Nperrat,mean=mu1_4,sd=sigma[i,1])
  mu1_5 <- alpha1[i,5] + beta1[i,5]*t[1:11]  
  wrep1_5[i,] <- rnorm(n=Nperrat,mean=mu1_5,sd=sigma[i,1])
  mu1_6 <- alpha1[i,6] + beta1[i,6]*t[12:22]  
  wrep1_6[i,] <- rnorm(n=Nperrat,mean=mu1_6,sd=sigma[i,1])
  mu1_7 <- alpha1[i,7] + beta1[i,7]*t[23:33]  
  wrep1_7[i,] <- rnorm(n=Nperrat,mean=mu1_7,sd=sigma[i,1])
  mu1_8 <- alpha1[i,8] + beta1[i,8]*t[34:44]  
  wrep1_8[i,] <- rnorm(n=Nperrat,mean=mu1_8,sd=sigma[i,1])
  
  mu2_1 <- alpha2[i,1] + beta2[i,1]*t[1:11]  
  wrep2_1[i,] <- rnorm(n=Nperrat,mean=mu2_1,sd=sigma[i,2])
  mu2_2 <- alpha2[i,2] + beta2[i,2]*t[12:22]  
  wrep2_2[i,] <- rnorm(n=Nperrat,mean=mu2_2,sd=sigma[i,2])
  mu2_3 <- alpha2[i,3] + beta2[i,3]*t[23:33]  
  wrep2_3[i,] <- rnorm(n=Nperrat,mean=mu2_3,sd=sigma[i,2])
  mu2_4 <- alpha2[i,4] + beta2[i,4]*t[34:44]  
  wrep2_4[i,] <- rnorm(n=Nperrat,mean=mu2_4,sd=sigma[i,2])

  mu3_1 <- alpha3[i,1] + beta3[i,1]*t[1:11]  
  wrep3_1[i,] <- rnorm(n=Nperrat,mean=mu3_1,sd=sigma[i,3])
  mu3_2 <- alpha3[i,2] + beta3[i,2]*t[12:22]  
  wrep3_2[i,] <- rnorm(n=Nperrat,mean=mu3_2,sd=sigma[i,3])
  mu3_3 <- alpha3[i,3] + beta3[i,3]*t[23:33]  
  wrep3_3[i,] <- rnorm(n=Nperrat,mean=mu3_3,sd=sigma[i,3])
  mu3_4 <- alpha3[i,4] + beta3[i,4]*t[34:44]  
  wrep3_4[i,] <- rnorm(n=Nperrat,mean=mu3_4,sd=sigma[i,3])  
  
}  
  wrep <- cbind(wrep1_1,wrep1_2,wrep1_3,wrep1_4,
                wrep1_5,wrep1_6,wrep1_7,wrep1_8,
                wrep2_1,wrep2_2,wrep2_3,wrep2_4,
                wrep3_1,wrep3_2,wrep3_3,wrep3_4)

# head(wrep)

color_scheme_set("brightblue")
ind <- sample(1:k, 25, replace=FALSE)
ppc_dens_overlay(wt, wrep[ind, ])
ind2 <- sample(1:k, 10, replace=FALSE)
ppc_boxplot(wt, wrep[ind2, ],notch=FALSE)
