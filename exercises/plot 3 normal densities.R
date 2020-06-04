rm(list=ls())
# Specify x values for which f(x) will be computed
x=seq(from=-10, to= 50, by = 0.01)
# specify mean and standard deviation for 1st dist.

# Note: the 1st dist. should be the one with 
# the smallest variance

mean_1 = 10
mean_2 = 10
mean_3 = 15

variance_1 = 1
variance_2 = 5
variance_3 = 10

mu = c(mean_1, mean_2, mean_3)
sigma = sqrt(c(variance_1, variance_2, variance_3))

p1 = dnorm(x, mean=mu[1], sd=sigma[1])
p2 = dnorm(x, mean=mu[2], sd=sigma[2])
p3 = dnorm(x, mean=mu[3], sd=sigma[3])

plot(x,p1,type="l",
   xlab=expression(italic(x)),
   ylab=expression(italic(f(x))),
   main = "normal distributions") # dist 1 is black
lines(x,p2,col=3)                # dist 2 is green
lines(x,p3,col=10)               # dist 3 is red

