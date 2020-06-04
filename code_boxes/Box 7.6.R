rm(list=ls())
par(mfrow=c(1,1))

x = c(0.245, 0.247,	0.285,	0.299,	0.327,  0.347,	0.356,	0.36,	  0.363,	0.364,	
      0.398, 0.4,	  0.409,	0.421,	0.432,	0.473,	0.509,	0.529,	0.561,	0.569,	
      0.594, 0.638,	0.656,	0.816,	0.853,	0.938,	1.036,  1.045)
y = c(0,0,1,1,1,1,0,1,0,1,0,1,0,1,0,1,1,1,0,0,1,1,1,1,1,1,1,1)
xgrid = seq(from=min(x), to=max(x),length=101)
corr_x = xgrid-mean(xgrid)
alpha = 1.177
beta = 6.294 
p_logit = exp(alpha + beta*corr_x)/(1+exp(alpha + beta*corr_x))
plot(xgrid,p_logit,xlab="grain size (mm)",ylab="P(presence)",pch=20,type = "l",
     ylim=c(0,1))
points(x,y,pch=19)