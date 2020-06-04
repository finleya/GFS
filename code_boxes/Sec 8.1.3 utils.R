library(Matrix)

crps <- function(y, y.hat, y.var){
    sd <- sqrt(y.var)
    y.std <- (y-y.hat)/sd
    mean(-sd*(1/sqrt(pi) - 2*dnorm(y.std) - y.std*(2*pnorm(y.std) - 1)))
}

rmspe <- function(y.hat, y){
    sqrt(mean((y-y.hat)^2))
}

svcDiag <- function(y, X, Z, beta, theta, w){

    ##where Z is the coluns of X that are space-varying
    q <- ncol(Z)
    Z.tild <- t(bdiag(as.list(as.data.frame(t(Z)))))
        
    n.samples <- nrow(theta)
    n <- length(y)
    
    out <- list()
    
    ##############
    ##DIC
    ##############  
    beta.mu <- apply(beta, 2, mean)
    w.mu <- matrix(apply(w, 1, mean), nrow(w), 1)
    
    tau.sq <- theta[,"tau.sq"]
    tau.sq.mu <- mean(tau.sq) 
    
    L <- sum(dnorm(y, X%*%beta.mu+(Z.tild%*%w.mu)[,1], sqrt(tau.sq.mu), log=TRUE))
    
    llSum <- 0
    for(i in 1:n.samples){
        llSum <- llSum+sum(dnorm(y, X%*%t(beta[i,,drop=FALSE])+(Z.tild%*%w[,i,drop=FALSE])[,1], sqrt(tau.sq[i]), log=TRUE))
    }
    
    P <- 2*(L-(1/n.samples*llSum))
    DIC <- data.frame(value=c(-2*(L-P), P, L))
    rownames(DIC) <- c("DIC", "pD", "L")
    out$DIC <- DIC
    
    ##############
    ##WAIC
    ##############
    LPPD <- 0
    P1 <- 0
    P2 <- 0
    
    for(i in 1:n){
        tau.sq <- theta[,"tau.sq"]
        L <- dnorm(y[i], as.vector(X[i,]%*%t(beta)+Z[i,]%*%w[((i-1)*q+1):((i-1)*q+q),]), sqrt(tau.sq), log=FALSE)       
        LPPD <- LPPD+log(mean(L))##log pointwise predictive density
        
        a <- log(mean(L))
        b <- mean(log(L))
        
        P1 <- P1 + 2*(a-b)
        P2 <- P2 + var(log(L))
    }
    
    WAIC.1 <- -2*(LPPD-P1)
    WAIC.2 <- -2*(LPPD-P2)
    
    out$WAIC <- data.frame(value=c(WAIC.1, WAIC.2, P1, P2, LPPD))
    row.names(out$WAIC) <- c("WAIC.1","WAIC.2","P.1","P.2","LPPD")
    
    out
}


nonspDiag <- function(y, X, beta, tau.sq){


    n.samples <- nrow(beta)
    n <- length(y)
    
    out <- list()
    
    ##############
    ##DIC
    ##############  
    beta.mu <- apply(beta, 2, mean)
    tau.sq.mu <- mean(tau.sq) 
    
    L <- sum(dnorm(y, X%*%beta.mu, sqrt(tau.sq.mu), log=TRUE))
    
    llSum <- 0
    for(i in 1:n.samples){
        llSum <- llSum+sum(dnorm(y, X%*%t(beta[i,,drop=FALSE]), sqrt(tau.sq[i]), log=TRUE))
    }
    
    P <- 2*(L-(1/n.samples*llSum))
    DIC <- data.frame(value=c(-2*(L-P), P, L))
    rownames(DIC) <- c("DIC", "pD", "L")
    out$DIC <- DIC
    
    ##############
    ##WAIC
    ##############
    LPPD <- 0
    P1 <- 0
    P2 <- 0
    
    for(i in 1:n){

        L <- dnorm(y[i], as.vector(X[i,]%*%t(beta)), sqrt(tau.sq), log=FALSE)       
        LPPD <- LPPD+log(mean(L))##log pointwise predictive density
        
        a <- log(mean(L))
        b <- mean(log(L))
        
        P1 <- P1 + 2*(a-b)
        P2 <- P2 + var(log(L))
    }
    
    WAIC.1 <- -2*(LPPD-P1)
    WAIC.2 <- -2*(LPPD-P2)
    
    out$WAIC <- data.frame(value=c(WAIC.1, WAIC.2, P1, P2, LPPD))
    row.names(out$WAIC) <- c("WAIC.1","WAIC.2","P.1","P.2","LPPD")
    
    out
}
