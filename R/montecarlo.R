plain.mc.quadratic <- function(N, mixture.param, mc.cores){
    get.T <- function(y, Sigma){
        return(sum(drop(y%*%Sigma)*y))
    }
    m0 <- mixture.param$m0
    m1 <- mixture.param$m1
    plus <- sqrt(m0/(m0 + m1)/m1)
    minus <- - sqrt(m1/(m0 + m1)/m0)
    x0 <- c(rep(minus, m0), rep(plus, m1))
    threshold <- 10^7
    if(N > threshold){
        N.sig <- rep(0, floor(N/threshold))
        if(length(N.sig) > 0){
            for(j in 1:length(N.sig)){
                N.sig[j] <- sum(mclapply(1:threshold, function(i) {x <- sample(x0, m0+m1); get.T(x, mixture.param$Sigma)}, mc.cores = mc.cores) >= mixture.param$T.obs)
            }
        }
        if(N - threshold*floor(N/threshold)>0){
            N.sig <- c(N.sig, sum(mclapply(1:(N - threshold*floor(N/threshold)), function(i) {x <- sample(x0, m0+m1); get.T(x, mixture.param$Sigma)}, mc.cores = mc.cores) >= mixture.param$T.obs))
        }
    }else{
        N.sig <- sum(mclapply(1:N, function(i) {x <- sample(x0, m0+m1); get.T(x, mixture.param$Sigma)}, mc.cores = mc.cores) >= mixture.param$T.obs)
    }
    return(sum(N.sig)/N)
}


plain.mc.linear <- function(N, lineardata, mc.cores){
    x0 <- lineardata$x0
    x1 <- lineardata$x1
    x <- c(x0, x1)
    m0 <- length(x0)
    m1 <- length(x1)
    plus <- sqrt(m0/(m0 + m1)/m1)
    minus <- - sqrt(m1/(m0 + m1)/m0)
    get.T <- function(x.perm){
        return(sum(x.perm[1:m0])*minus + sum(x.perm[(m0+1):(m0+m1)])*plus)
    }
    lineardata$T.obs <- get.T(c(x0, x1))

    threshold <- 10^7
    if(N > threshold){
        N.sig <- rep(0, floor(N/threshold))
        if(length(N.sig) > 0){
            for(j in 1:length(N.sig)){
                N.sig[j] <- sum(abs(unlist(mclapply(1:threshold, function(i) {x.perm <- sample(x, m0+m1); get.T(x.perm)}, mc.cores = mc.cores))) >= abs(lineardata$T.obs))
            }
        }
        if(N - threshold*floor(N/threshold)>0){
            N.sig <- c(N.sig, sum(abs(unlist(mclapply(1:threshold, function(i) {x.perm <- sample(x, m0+m1); get.T(x.perm)}, mc.cores = mc.cores))) >= abs(lineardata$T.obs)))
        }
    }else{
        N.sig <- sum(abs(unlist(mclapply(1:N, function(i) {x.perm <- sample(x, m0+m1); get.T(x.perm)}, mc.cores = mc.cores))) >= abs(lineardata$T.obs))
    }
    return(sum(N.sig)/N)
}
