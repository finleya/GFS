rm(list=ls())
library(viridis)
library(MBA)
library(fields)

matern <- function(d, phi, nu){
    r <- (d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi, nu=nu)
    r[d == 0] <- 1
    r
}

eff.rang <- function(r, phi, nu){

    f <- function(d, phi, nu){
        r - (d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi, nu=nu)
    }
    
    uniroot(f, c(0.01, 1000), phi=phi, nu=nu)$root
 
}

#####################################
##Code for Figure 8.1
#####################################
d <- seq(0, 1, length.out=1000)

phi <- seq(3/1, 3/0.1, length.out=6)

nu <- c(0.25,0.5,1.5)

for(j in nu){
    
    
    phi.nu <- expand.grid(phi, j)
    
    col <- magma(nrow(phi.nu)+1); col <- rev(col[1:nrow(phi.nu)])
    r.eff <- rep(0, nrow(phi.nu))
    
    #pdf(file=paste0("figures/matern-cor-nu",gsub("[.]","_",j),".pdf"))
    #par(xpd=TRUE) 
    plot(0, 0, xlim=c(0,1), ylim=c(0,1), typ="n", xlab="Distance", ylab="", cex.lab=1.5, cex=1.5, cex.axis=1.5)
    text(-0.15, 0.5, expression(rho), cex=1.75)
    
    for(i in 1:nrow(phi.nu)){
        r <- matern(d, phi.nu[i,1], phi.nu[i,2])
        r.eff[i] <- eff.rang(0.05, phi.nu[i,1], phi.nu[i,2])
        lines(d, r, col=col[i], lw=3)
    }
    legend(0.6, 0.9, legend=round(phi,1), col=col, lty=1, lw=3, bty="n", cex=1.5)
    
    legend(0.83, 0.9, legend=round(r.eff,2), bty="n", cex=1.5)
    
    text(0.8, 0.9, expression(phi1), cex=1.75)
    text(0.97, 0.95, "Eff.", cex=1.5)
    text(0.96, 0.9, "range", cex=1.5)
    #dev.off()

}

#####################################
##Code for Box 8.1 and Figure 8.2
#####################################

## rmvn <- function(n, mu=0, V = matrix(1)){
##   p <- length(mu)
##   if(any(is.na(match(dim(V),p))))
##     stop("Dimension problem!")
##   D <- chol(V)
##   t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
## }

##Simplified version of rmvn to keep the box uncluttered
rmvn  <- function(mu, Sigma){
    mu + t(chol(Sigma))%*%rnorm(length(mu))
}

set.seed(1)

n <- 2000
coords <- cbind(runif(n,0,1), runif(n,0,1))
sigma.sq <- 1
phi <- 3/0.1

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(rep(0,n), sigma.sq*R)

surf <- mba.surf(cbind(coords, w), no.X=100, no.Y=100, extend=TRUE)$xyz.est

#pdf(file=paste0("figures/sim-surf-effr0_1.pdf"))
image.plot(surf, col=viridis(100), axes=F, useRaster=TRUE, xaxs="r", yaxs="r",  xlab="Easting", ylab="Northing", cex.lab=1.5, cex=1.5, cex.axis=1.5)
axis(1)
axis(2)
#dev.off()

#
phi <- 3/1
D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(rep(0,n), sigma.sq*R)

surf <- mba.surf(cbind(coords, w), no.X=100, no.Y=100, extend=TRUE)$xyz.est

#pdf(file=paste0("figures/sim-surf-effr1_0.pdf"))
image.plot(surf, col=viridis(100), axes=F, useRaster=TRUE, xaxs="r", yaxs="r",  xlab="Easting", ylab="Northing", cex.lab=1.5, cex=1.5, cex.axis=1.5)
axis(1)
axis(2)
#dev.off()
