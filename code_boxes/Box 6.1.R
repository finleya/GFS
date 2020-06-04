x <- seq(from=1, to=10, by=0.02)
n <- length(x)
alpha <- 1
beta <- 2
yhat <- alpha + beta*x
res1 <- rnorm(n, mean=0, sd=1)
y1 <- yhat + res1
res2 <- rnorm(n, mean=0, sd=1)*0.5*x
y2 <- yhat + res2
maxy = max(c(max(y1),max(y2)))
par(mfrow=c(1,2))
plot(x,y1,ylab="Y", xlab="X", pch=".",
     main="homogeneous variance",ylim=c(0,maxy))
plot(x,y2,ylab="Y", xlab="X", pch=".",
     main="heterogeneous variance",ylim=c(0,maxy))

