x <- seq(from=-5,to=5,by=0.1)
var <- 1.0E6
std<-sqrt(var)
fx <- dnorm(x, mean=0, sd=std, log = FALSE)
par(mar=c(5.1,5.1,4.1,2.1),mfrow=c(2,2))
plot(x,fx,ylab = "f(x)",xlab= "x",type="l",ylim=c(0,0.01),
     main=expression(paste("N","(0, ", "1.0 x ", 10^6,") prior")))
alpha0 <- 0.001 
beta0 <- 0.001
y <- seq(from=0.001, to=25.0, by=0.01)
fy <- dgamma(y, shape = alpha0, rate = beta0, log = FALSE)
y_up <- max(y)
y_lo <- min(y)
plot(y,fy,type="l",xlab="x",ylab="f(x)",
     xlim =c(0,25),ylim=c(y_lo, y_up),
     main=expression(paste("Ga(0.001, 0.001) prior")))
