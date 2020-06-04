rm(list=ls())
alpha_1 = 2
beta_1 = 2
alpha_2 = 5
beta_2 = 5
alpha = c(alpha_1,alpha_2)
beta = c(beta_1,beta_2)

mean = alpha/beta
mean
var = alpha/(beta*beta)
var
x_lo = mean - 5*sqrt(var)
x_lo = max(c(0.001,min(x_lo)))
x_up = mean + 5*sqrt(var)
x_up = max(x_up)
x = seq(from = x_lo, to = x_up, by = 0.01)

y1 = dgamma(x, shape = alpha[1], rate = beta[1], log = FALSE)
y_up1 = max(y1)
y_lo1 = min(y1)

y2 = dgamma(x, shape=alpha[2], rate=beta[2], log=FALSE)

y_up2 = max(y2)
y_lo2 = min(y2)

y_up = max(c(y_up1,y_up2))
y_lo = min(c(y_lo1,y_lo2))

plot(x,y1,type="l",xlab="x",ylab="f(x)",
     xlim=c(x_lo,x_up),ylim=c(y_lo,y_up),
    main="gamma distrbutions")   # distribution 1 is black

lines(x,y2,col=3)     # distribution 2 is green
