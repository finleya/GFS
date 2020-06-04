rm(list=ls())
#
#  Be sure to set working directory!!
#
setwd("//Users//edwingreen//Documents//Bayes text")

k = 25000  # k = posterior sample size
alpha = matrix(0,nrow=k,ncol=1)
beta = matrix(0,nrow=k,ncol=1)
sigma = matrix(0,nrow=k,ncol=1)

# read posterior samples

y = read.table("rat diet 2.out",header=FALSE,row.names = NULL)
alpha[] = y[1:k,2]
beta[] = y[((6*k)+1):(7*k),2]
sigma = y[((12*k)+1):(13*k),2]
rm(y)
yhat = rep(0, nrow = k, ncol = 1)

x = 70

#for (i in 1:k){
  muhat = alpha[,1] + beta[,1]*x  
  yhat = rnorm(n=k,mean=muhat[],sd=sigma[])
#}

stats = matrix(0,nrow=1,ncol=4)
stats[1,1]=mean(yhat)
stats[1,2]=sqrt(var(yhat))
stats[1,3]= quantile(yhat,0.025)
stats[1,4]= quantile(yhat,0.975)

est = data.frame(stats)
colnames(est)=c("mean","std dev","0.025 pct","0.975 pct")
est
 
d = density(yhat)
plot(d, ylab = "density",
     xlab= "predicted weight at age 70 for rat 1",main=" ")
