rm(list=ls())
alpha <- 143.001  # change to desired value  
beta <- 10.001    # change to desired value 
mean <- alpha/beta
sd <- sqrt(alpha/(beta^2))
lo_lim_x <- mean - 5 * sd
up_lim_x <- mean + 5 * sd 
inc <- (up_lim_x - lo_lim_x)/10000
x <- seq(from = lo_lim_x, to = up_lim_x, by = inc)
ly <- (alpha-1)*log(x) + alpha*log(beta) - x*beta - lgamma(alpha)
y <- exp(ly)
plot(x, y, type ="l", xlab ="x", ylab="f(x)", xlim =
       c(lo_lim_x,up_lim_x))