rm(list=ls())
par(mfrow=c(1,1))
set.seed(17)
d = 15
h = 75
d2h = d*d*h
k =15000  # posterior sample size
# read posterior samples
# set to directory where posterior samples are stored

setwd("  ") 
# read in posterior samples (change file name as 
# appropriate)
y = read.table("trees.out",header=FALSE,row.names = NULL)
a = y[1:k,2]
b = y[(k+1):(2*k),2]
#rsq = y[((2*k)+1):(3*k),2]
sigma = y[((3*k)+1):(4*k),2]
rm(y)
ynew=rep(0,k)
for (i in 1:k){
  mu = a[i] + b[i]*(d2h)    
  ynew[i] = rnorm(1,mean=mu,sd=sigma[i])
}
d <- density(ynew)
plot(d, ylab = "density",
     xlab= expression(paste( "Volume"," (","ft"^"3",")")),
     main=" ")
stats=matrix(0, nrow=1, ncol=5)
stats[1,1] = mean(ynew)
stats[1,2] = sqrt(var(ynew))
stats[1,3] = quantile(ynew,0.025)
stats[1,4] = quantile(ynew,0.5)
stats[1,5] = quantile(ynew,0.975)
est = data.frame(stats)
colnames(est)=c("mean","std dev","0.025 pct","median","0.975 pct")
est