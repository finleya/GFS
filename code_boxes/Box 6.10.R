rm(list=ls())
set.seed(17)
#
#  Be sure to set working directory!!
#
#setwd("  ")

k <- 25000  # k = posterior sample size
alpha_mu <- matrix(0,nrow=k,ncol=1)
alpha_sig <- matrix(0,nrow=k,ncol=1)
beta_mu <- matrix(0,nrow=k,ncol=1)
beta_sig <- matrix(0,nrow=k,ncol=1)
sigma <- matrix(0,nrow=k,ncol=1)
muhat <- matrix(0,nrow=k,ncol=1)
# read posterior samples

y <- read.table("rat diet 2.out",header=FALSE,row.names = NULL)

alpha_mu[] <- y[((4*k)+1):(5*k),2]
alpha_sig[] <- y[((5*k)+1):(6*k),2]
beta_mu[] <- y[((10*k)+1):(11*k),2]
beta_sig[] <- y[((11*k)+1):(12*k),2]
sigma <- y[((12*k)+1):(13*k),2]
rm(y)
yhat <- rep(0, nrow = k, ncol = 1)
x <- 70

for (i in 1:k){
alpha <- rnorm(n=1,mean=alpha_mu[i],sd=alpha_sig[i])
beta <- rnorm(n=1,mean=beta_mu[i],sd=beta_sig[i])
muhat[i] <- alpha + beta*x  
yhat[i] <- rnorm(n=1,mean=muhat[i],sd=sigma[i])
}

stats <- matrix(0,nrow=1,ncol=4)
stats[1,1] <- mean(yhat)
stats[1,2] <- sqrt(var(yhat))
stats[1,3] <- quantile(yhat,0.025)
stats[1,4] <- quantile(yhat,0.975)

est <- data.frame(stats)
colnames(est) <- c("mean","std dev","0.025 pct","0.975 pct")
est

d <- density(yhat,from=quantile(yhat,0.01),
             to=quantile(yhat,0.99),kernel="optcosine")

plot(d, ylab = "density",
     xlab= "predicted weight, age 70 for new, independent rat",main=" ",
     xlim=c(0,1200))