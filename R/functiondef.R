normalize <- function(a){
    a <-  a - mean(a)
    a <- a/sqrt(sum(a^2))
    return(a)
}


get.exact.regular <- function(m0, m1, Sigma, T.obs){
    get.T <- function(y, Sigma){
        return(sum(drop(y%*%Sigma)*y))
    }
    A <- combn(1:(m0+m1), m1)
    plus <- sqrt(m0/(m0 + m1)/m1)
    minus <- - sqrt(m1/(m0 + m1)/m0)
    y.mat <- matrix(minus, choose(m0+m1, m1), m0 + m1)
    for(i in 1:dim(A)[2]){
        y.mat[i, A[,i]] <- plus
    }
    T <- apply(y.mat, 1, function(x) get.T(x, Sigma))
    pval <- length(which(abs(T) >= abs(T.obs)))/length(T)
    return(pval)
}



get.mapping <- function(r, n, b){
    comb <- rep(NA, r)
    m <- r
    i <- 1
    while(m > 0){
        if(b > choose(n-i, m - 1)){
            b <- b - choose(n-i, m - 1)
            comb[i] <- 1
        }else{
            comb[i] <- 0
            m <- m - 1
        }
        i <- i + 1
    }
    if(i <= n)
        comb[i:n] <- 1
    comb <- which(comb == 0)
    return(comb)
}



next.comb <- function(r, n, s){
    m <- r
    max.val <- n
    while(m > 0 && s[m] == max.val){
        m <- m - 1
        max.val <- max.val - 1
    }
    if(m <= 0){
        return(FALSE)
    }
    s[m] <- s[m] + 1
    if(m+1 <= r){
        ## s[(m+1):r] <- s[m:(r-1)] + 1
        for(j in (m+1):r){
            s[j] <- s[j-1] + 1
        }
    }
    return(s)
}





get.exact.par <- function(m0, m1, k, N.par, Sigma, T.obs){
    get.T <- function(y, Sigma){
        return(sum(drop(y%*%Sigma)*y))
    }
    N <- choose(m0 + m1, m1)
    if(k == 1){
        start <- 1
    }else{
        start <- floor((k - 1)/N.par*N) + 1
    }

    if(k == N.par){
        end <- N
    }else{
        end <- floor(k/N.par*N)
    }

    plus <- sqrt(m0/(m0 + m1)/m1)
    minus <- - sqrt(m1/(m0 + m1)/m0)
    comb <- get.mapping(m0, m0 + m1, start)
    T <- rep(NA, end - start + 1)

    for(i in 1:length(T)){
        y <- rep(plus, m0 + m1)
        y[comb] <- minus
        T[i] <- get.T(y, Sigma)
        comb <- next.comb(m0, m0 + m1, comb)
    }
    return(list(N = end - start + 1, sig.count = length(which(abs(T) >= abs(T.obs)))))
}




get.a <- function(S, C){
    N <- nrow(S)
    d <- ncol(C) + 1
    Scumprod <- matrix(1, N, d)
    for(j in 2:(d-1)){
        ## Scumprod[,j] is the prod_{i=1}1^{j-1}
        Scumprod[,j] <- Scumprod[,j-1]*S[,j-1]
    }
    Scumprod[,d] <- Scumprod[,d-1]*S[,d -1]
    a <- matrix(NA, N, d)
    for(j in 1:(d-1)){
        a[,j] <- Scumprod[,j]*C[,j]
    }
    a[,d] <- Scumprod[,d]
    return(a)
}

get.mixture.param <- function(delta = 1, m0 = 10, m1 = 10, L = 4, G = 5, s = 0, exact = FALSE, N.par = 1, mc.cores = 20, dist = "norm"){
    id <- 0
    seed <- 2
    ## dist <- 1
    is.rep <- 3
    N.mc <- 10^5
    mc.rep <- 0
    N.is <- 10^5
    param.test <- matrix(as.numeric(c(id, seed, delta, m0, m1, G, N.mc, mc.rep, N.mc, is.rep)), nrow = 1)
    param.test <- data.frame(param.test)
    names(param.test) <- c("id", "seed", "delta", "m0", "m1", "G", "N.mc", "mc.rep", "N.is", "is.rep")
    mixture.param <- list(s = 0,m0=m0, m1=m1, p = m0 + m1, d = m0 + m1,  df = rep(10, L), J = L + 1)
    if(dist == "norm"){
        x0 <- matrix(nrow = m0, rnorm(m0*G, 0, 1), byrow = TRUE)
        x1 <- matrix(nrow = m1, rnorm(m1*G, delta, 1), byrow = TRUE)
    }else if(dist == "exp"){
        x0 <- matrix(nrow = m0, rexp(m0*G, 1), byrow = TRUE)
        x1 <- matrix(nrow = m1, rexp(m1*G, 1 + delta), byrow = TRUE)

    }
    d <- m0 + m1
    x <- rbind(x0, x1)
    x <- apply(x, 2, normalize)
    Sigma <- x%*%t(x)
    tmp <- svd(Sigma)
    lambda <- tmp$d
    eta <- tmp$u
    result <- list()
    plus <- sqrt(m0/(m0 + m1)/m1)
    minus <- - sqrt(m1/(m0 + m1)/m0)
    get.T <- function(y){
        return(sum(drop(y%*%Sigma)*y))
    }
    y.obs <- c(rep(minus, m0), rep(plus, m1))
    T.obs <- get.T(y.obs)
    lambda[lambda < 10^(-10)] <- 0
    mixture.param$lambda <- lambda
    mixture.param$eta <- eta
    mixture.param$Sigma <- Sigma
    mixture.param$T.obs <- T.obs
    mixture.param$denominator <- 1
    ## mixture.param$s <- s
    if(exact){
        if(N.par > 1){
            ## print("estimated hours")
            ## system.time(mclapply(1:20, function(i) get.exact.par(m0, m1, i, 10^4, Sigma, T.obs), mc.cores = 20))*500/60/60
            t0 <- proc.time()[[3]]
            a <- mclapply(1:N.par, function(i) get.exact.par(m0, m1, i, N.par, Sigma, T.obs), mc.cores = mc.cores)
            aa <- colSums(matrix(unlist(a), ncol = 2, byrow = TRUE))
            exactpval <- aa[2] / aa[1]
            time <- (proc.time()[[3]] - t0)/60
        }else{
            t0 <- proc.time()[[3]]
            exactpval <- get.exact.regular(m0, m1, Sigma, T.obs)
            time <- (proc.time()[[3]] - t0)/60
        }
        mixture.param$exact <- exactpval
        mixture.param$time <- time
    }
    return(mixture.param)
}

rtrunc.chisq <- function(n, a, b, df){
    u <- runif(n, 0, 1)
    return(qchisq(u*(pchisq(b, df) - pchisq(a, df)) + pchisq(a, df), df))
}

## Requires R >= 3.2.3
rtrunc.beta <- function(n, a, b, shape1, shape2){
    if(b<a){
        print(paste0("a = ", a, " b = ", b))
        warning("In truncated beta lower limit > upper limit")
    }
    u <- runif(n, 0, 1)
    if(b == 1){
        y <- qbeta(log(1 - u) + pbeta(a, shape1, shape2, lower.tail = FALSE, log.p = TRUE),shape1, shape2, lower.tail = FALSE, log.p = TRUE)
    }else if(a == 0){
        y <- qbeta(log(1 - u) + pbeta(b, shape1, shape2, log.p = TRUE),shape1, shape2, log.p = TRUE)
    }else{
        stop("need to write formula for sampling general truncated beta")
    }
    return(y)
}

project.y.to.discrete <- function(z, m0, m1){
    plus <- sqrt(m0/(m0 + m1)/m1)
    minus <- - sqrt(m1/(m0 + m1)/m0)
    plus.index <- matrix(apply(z, 1, function(x) order(x, decreasing = TRUE)[1:m1]), nrow = m1)
    y.mat <- matrix(minus, dim(z)[1], m0 + m1)
    ## qz suspects this is slow because of the for loop
    ## try matrix index...
    for(i in 1:dim(z)[1]){
        y.mat[i, plus.index[,i]] <- plus
    }
    return(y.mat)
}

test.rtrunc.beta <- function(){
    a <- 0.9
    b <- 1
    y <- rtrunc.beta(10^5, a, b, 1/2, 10/2)
    hist(y, freq = FALSE, breaks = 100)
    x <- seq(0,1,length.out=10000)
    lines(x, dbeta(x, 1/2, 10/2)/(pbeta(b, 1/2, 10/2) - pbeta(a, 1/2, 10/2)))
}





matrix.power <- function(A,k){
    a <- svd(A)
    return(a$u %*% diag((a$d)^k) %*% t(a$v))
}


rp.large <- function(N, mixture.param){
    d <- mixture.param$d
    z <- matrix(NA, N, d - 2)
    for(j in 1:(d-2)){
        z[,j] <- rbeta(N, 1/2, (d-j)/2)
    }
    signs.random <- matrix(sample(c(-1, 1), N*(d-2), replace = TRUE), N, d-2)
    z <- sqrt(z)*signs.random
    C <- matrix(NA, N, d - 1)
    S <- matrix(NA, N, d - 1)
    C[,1:(d-2)] <- z
    S[,1:(d-2)] <- sqrt(1 - z^2)
    phi <- runif(N, 0, 1)
    S[,d-1] <- sin(2*pi*phi)
    C[,d-1] <- cos(2*pi*phi)
    a <- get.a(S, C)
    return(a)
}


dp.large <- function(a, mixture.param){
    N <- nrow(a)
    d <- mixture.param$d
    r <- t(apply(matrix(a[,d:1]^2, nrow = nrow(a)), 1, cumsum))[,d:2]
    C <- a[,1:(d-1)]/sqrt(r)
    C[which(is.nan(C))] <- 0
    C <- matrix(C, nrow = nrow(a))
    density <- apply(matrix(sapply(1:(d-2), function(j) pi/beta(1/2,(d-j)/2)*(1-C[,j]^2)^((d-j-1)/2)), nrow = nrow(a)), 1, prod)
    return(density)
}


rp.small <- function(N, mixture.param){
    d <- mixture.param$d
    z <- matrix(NA, N, d)
    r <- matrix(NA, N, d)
    a <- matrix(NA, N, d)
    r[,d] <- 1
    for(j in d:3){
        z[,j] <- rbeta(N, 1/2, (j-1)/2)
        signs.random <- sample(c(-1, 1), N, replace = TRUE)
        z[,j] <- sqrt(z[,j])*signs.random
        a[,j] <- r[,j]*z[,j]
        r[,j-1] <- sqrt(r[,j]^2 - a[,j]^2)
    }
    phi <- runif(N, 0, 1)
    a[,2] <- r[,2]*cos(2*pi*phi)
    a[,1] <- r[,2]*sin(2*pi*phi)
    return(a)
}

dp.small <- function(a, mixture.param){
    N <- nrow(a)
    d <- mixture.param$d
    r <- t(apply(matrix(a^2,nrow = nrow(a)), 1, cumsum))[,2:d]
    C <- matrix(NA, N, d)
    C[,2:d] <- a[,2:d]/sqrt(r)
    C[which(is.nan(C))] <- 0
    C <- matrix(C, nrow = nrow(a))
    density <- apply(matrix(sapply(3:d, function(k) pi/beta(1/2,(k-1)/2)*(1-C[,k]^2)^((k-1-1)/2)), nrow = nrow(a)), 1, prod)
    return(density)
}


rp.small.reduced <- function(N, mixture.param){
    d <- mixture.param$d
    r <- length(which(mixture.param$lambda != 0))
    if(r <= 1 | r >= d){
        stop("this method is for 1 < r < d")
    }
    S1 <- rchisq(N, r)
    S2 <- rchisq(N, d - r)
    b <- matrix(NA, N, r)
    C <- matrix(NA, N, r)
    R <- matrix(NA, N, r)
    if(r >= 3){
        for(k in r:3){
            b[,k] <- rbeta(N, 1/2, (k-1)/2)
            signs.random <- sample(c(-1, 1), N, replace = TRUE)
            C[,k] <- sqrt(b[,k])*signs.random
        }
    }
    phi <- runif(N, 0, 1)
    C[,2] <- cos(2*pi*phi)
    return(list(S1 = S1, S2 = S2, C = C))
}




dp.small.reduced <- function(sample.list, mixture.param){
    d <- mixture.param$d
    r <- length(which(mixture.param$lambda != 0))
    if(r <= 1 | r >= d){
        stop("this method is for 1 < r < d")
    }
    S1 <- sample.list$S1
    S2 <- sample.list$S2
    C <- sample.list$C
    N <- length(S1)
    density <- dchisq(S1, r)*dchisq(S2, d-r)
    if(r >= 3){
        density <- density*apply(matrix(sapply(3:r, function(k) 1/beta(1/2,(k-1)/2)*(1-C[,k]^2)^((k-1)/2 -1)), nrow = N), 1, prod)
    }
    density <- density*1/(2*pi*sqrt(1-C[,2]^2))
    return(density)
}




rq.large <- function(N, j, mixture.param){
    ## s <- mixture.param$s[j]
    d <- mixture.param$d
    lambda <- mixture.param$lambda/(mixture.param$T.obs/mixture.param$denominator[j])
    w <- matrix(0, N, d)
    if(lambda[1] - lambda[2] > 10^(-10)){
        w[,1] <- (1 - lambda[2])/(lambda[1] - lambda[2])
    }
    w[which(w[,1]<0),1] <- 0
    z <- matrix(NA, N, d)
    a <- matrix(NA, N, d)
    r <- matrix(NA, N, d)
    r[,1] <- 1
    for(k in 1:(d-2)){
        z[,k] <- sapply(1:N, function(i) rtrunc.beta(1, a = w[i, k], b = 1, shape1 = 1/2, shape2 = (d-k)/2))
        signs.random <- sample(c(-1, 1), N, replace = TRUE)
        z[,k] <- sqrt(z[,k])*signs.random
        a[,k] <- r[,k]*z[,k]
        r[,k+1] <- sqrt(r[,k]^2 - a[,k]^2)
        numerator <- (1 - colSums(matrix(t(a[,1:k]^2)*lambda[1:k], nrow = k)) - (1 - rowSums(matrix(a[,1:k]^2, ncol = k)))*lambda[k + 2])
        denominator <- (lambda[k+1] - lambda[k+2])*r[,k+1]^2## abs(1-rowSums(matrix(a[,1:k]^2, ncol = k)))
        index <- which(abs(denominator) > 10^(-10))
        w[index, k+1] <- numerator[index]/denominator[index]
        w[which(w[,k+1]<0),k+1] <- 0
    }
    phi <- sapply(1:N, function(i) runif(1, -acos(sqrt(w[i,d-1]))/(2*pi), acos(sqrt(w[i,d-1]))/(2*pi)))
    index <- sample(1:N, round(N/2))
    phi[index] <- 1/2 + phi[index]
    a[,d-1] <- r[,d-1]*cos(2*pi*phi)
    a[,d] <- r[,d-1]*sin(2*pi*phi)
    return(a)
}


dq.large <- function(a, j, mixture.param){
    density <- rep(0, nrow(a))
    lambda <- mixture.param$lambda/(mixture.param$T.obs/mixture.param$denominator[j])
    index <- which(apply(t(a^2)*lambda, 2, sum) >= 1)
    if(length(index) > 1){
        anew <- matrix(a[index, ], nrow = length(index))
        ## s <- mixture.param$s[j]
        N <- nrow(anew[, ])
        d <- mixture.param$d
        lambda <- mixture.param$lambda/(mixture.param$T.obs/mixture.param$denominator[j])
        r <- t(apply(matrix(anew[,d:1]^2,nrow = nrow(anew)), 1, cumsum))[,d:2]
        C <- anew[,1:(d-1)]/sqrt(r)
        C[which(is.nan(C))] <- 0
        C <- matrix(C, nrow = nrow(anew))
        w <- matrix(NA, N, d - 1)
        w[,1] <- (1 - lambda[2])/(lambda[1] - lambda[2])
        w[which(is.nan(w[,1])|w[,1]<0),1] <- 0
        w <- matrix(w, nrow = nrow(anew))
        for(k in 1:(d-2)){
            numerator <- (1 - colSums(matrix(t(anew[,1:k]^2)*lambda[1:k], nrow = k)) - (1 - rowSums(matrix(anew[,1:k]^2, ncol = k)))*lambda[k + 2])
            numerator[which(numerator < 10^(-10))] <- -10^(-10) #avoid the case where numerator is almost zero
            w[, k+1] <- numerator/(lambda[k+1] - lambda[k+2])/abs(1-rowSums(matrix(anew[,1:k]^2, ncol = k)))
            w[which(is.nan(w[,k+1])|w[,k+1]<0),k+1] <- 0
        }
        density[index] <- apply(matrix(sapply(1:(d-2), function(k) pi/(beta(1/2,(d-k)/2)*(1-pbeta(w[,k], 1/2, (d-k)/2)))*(1-C[,k]^2)^((d-k-1)/2)), nrow = nrow(anew)), 1, prod) * pi/(2*acos(sqrt(w[,d-1])))
    }
    return(density)
}

rq.small <- function(N, j, mixture.param){
    if(N == 0){
        return(NULL)
    }
    ## s <- mixture.param$s
    d <- mixture.param$d
    lambda <- mixture.param$lambda/(mixture.param$T.obs/mixture.param$denominator[j])
    w <- matrix(1, N, d)
    if(lambda[1] - lambda[d] > 10^(-10)){
        w[,d] <- (1 - lambda[1])/(lambda[d] - lambda[1])
    }
    w[which(w[,d]>1),d] <- 1
    z <- matrix(NA, N, d)
    a <- matrix(NA, N, d)
    r <- matrix(NA, N, d)
    r[,d] <- 1
    for(k in d:3){
        if(any(w[,k] < 0)){
            stop("w[,k]< 0")
        }
        z[,k] <- sapply(1:N, function(i) rtrunc.beta(1, a = 0, b = w[i, k], shape1 = 1/2, shape2 = (k-1)/2))
        if(any(z[,k] < 0)){
            stop("z[,k] < 0")
        }
        signs.random <- sample(c(-1, 1), N, replace = TRUE)
        z[,k] <- sqrt(z[,k])*signs.random
        a[,k] <- r[,k]*z[,k]
        r[,k-1] <- sqrt(r[,k]^2 - a[,k]^2)
        numerator <- (1 - colSums(matrix(t(a[,k:d]^2)*lambda[k:d], nrow = d-k + 1)) - (1 - rowSums(matrix(a[,k:d]^2, ncol = d-k+1)))*lambda[1])
        if(length(which(numerator > 10^(-10))) > 0){
            print(paste0("k = ", k))
            index <- which(numerator > 0)
            print("a[index, k:d]")
            print(a[index, k:d])
            print("numerator[index]")
            print(numerator[index])
            stop("numerator > 0")
        }
        ## numerator[numerator > 10^(-10)] <- 0
        ## numerator[which(numerator < 10^(-10))] <- -10^(-10) #avoid the case where numerator is almost zero
        denominator <- (lambda[k-1] - lambda[1])*r[,k-1]^2## abs(1-rowSums(matrix(a[,k:d]^2, ncol = d-k+1)))
        index <- which(abs(denominator) > 10^(-10))
        w[index, k-1] <- numerator[index]/denominator[index]
        if(any(w[, k-1] < 0)){
            print(paste0("k = ", k))
            index <- which(w[,k-1] < 0)
            print("denominator[index]")
            print(denominator[index])
            print("numerator[index]")
            print(numerator[index])
        }
        w[which(w[,k-1]>1),k - 1] <- 1
    }
    phi <- sapply(1:N, function(i) runif(1, acos(sqrt(w[i,2]))/(2*pi), 1/2 - acos(sqrt(w[i,2]))/(2*pi)))
    index <- sample(1:N, round(N/2))
    phi[index] <- 1/2 + phi[index]
    a[,2] <- r[,2]*cos(2*pi*phi)
    a[,1] <- r[,2]*sin(2*pi*phi)
    return(a)
}

dq.small <- function(a, j, mixture.param){
    if(is.null(a)){
        return(NULL)
    }
    density <- rep(0, nrow(a))
    lambda <- mixture.param$lambda/(mixture.param$T.obs/mixture.param$denominator[j])
    nonzeroindex <- which(apply(t(a^2)*lambda, 2, sum) >= 1)
    if(length(nonzeroindex) > 1){
        anew <- matrix(a[nonzeroindex, ], nrow = length(nonzeroindex))
        ## s <- mixture.param$s
        N <- nrow(anew[, ])
        d <- mixture.param$d
        r <- t(apply(matrix(anew^2,nrow = nrow(anew)), 1, cumsum))[,2:d]
        C <- matrix(NA, N, d)
        C[,2:d] <- anew[,2:d]/sqrt(r)
        C[which(is.nan(C))] <- 0
        C <- matrix(C, nrow = nrow(anew))
        w <- matrix(1, N, d)
        if(lambda[1] - lambda[d] > 10^(-10)){
            w[,d] <- (1 - lambda[1])/(lambda[d] - lambda[1])
        }
        w <- matrix(w, nrow = nrow(anew))
        for(k in d:3){
            numerator <- (1 - colSums(matrix(t(anew[,k:d]^2)*lambda[k:d], nrow = d-k + 1)) - (1 - rowSums(matrix(anew[,k:d]^2, ncol = d-k+1)))*lambda[1])
            ## numerator[which(numerator < 10^(-10))] <- -10^(-10) #avoid the case where numerator is almost zero
            denominator <- (lambda[k-1] - lambda[1])*abs(1-rowSums(matrix(anew[,k:d]^2, ncol = d-k+1)))
            index <- which(abs(denominator) > 10^(-10))
            w[index, k-1] <- numerator[index]/denominator[index]
            w[which(w[,k-1]>1),k - 1] <- 1
        }
        density[nonzeroindex] <- apply(matrix(sapply(3:d, function(k) pi/(beta(1/2,(k-1)/2)*(pbeta(w[,k], 1/2, (k-1)/2)))*(1-C[,k]^2)^((k-1-1)/2)), nrow = nrow(anew)), 1, prod) * pi/(pi - 2*acos(sqrt(w[,2])))
    }
    return(density)
}
rq.small.reduced <- function(N, j, mixture.param){
    d <- mixture.param$d
    lambda <- mixture.param$lambda/(mixture.param$T.obs/mixture.param$denominator[j])
    ## s <- mixture.param$s
    r <- length(which(mixture.param$lambda != 0))
    if(r <= 1 | r >= d){
        stop("this method is for 1 < r < d")
    }
    S1 <- rchisq(N, r)
    S2 <- rtrunc.chisq(N, 0, (lambda[1] - 1)*S1, d - r)
    lp <- matrix(rep(lambda[1:r], N), ncol = r, byrow = TRUE)
    lp <- lp/(S1 + S2)
    b <- matrix(NA, N, r)
    C <- matrix(NA, N, r)
    R <- matrix(NA, N, r)
    z <- matrix(NA, N, r)
    w <- matrix(1, N, r)
    numerator <- (1 - S1*lp[,1])
    denominator <- (lp[,r] - lp[,1])*S1
    index <- which(abs(denominator) > 10^(-15))
    w[index,r] <- numerator[index]/denominator[index]
    w[which(w[,r] >= 1), r] <- 1
    R[,r] <- sqrt(S1)
    if(r >= 3){
        for(k in r:3){
            b[,k] <- sapply(1:N, function(i) rtrunc.beta(1, a = 0, b = w[i, k], shape1 = 1/2, shape2 = (k-1)/2))
            signs.random <- sample(c(-1, 1), N, replace = TRUE)
            C[,k] <- sqrt(b[,k])*signs.random
            z[,k] <- R[,k]*C[,k]
            R[,k-1] <- sqrt(R[,k]^2 - z[,k]^2)
            numerator <- (1 - rowSums(z[,k:r, drop = FALSE]^2*lp[,k:r, drop = FALSE]) - (S1 - rowSums(z[,k:r, drop = FALSE]^2))*lp[,1, drop = FALSE])
            denominator <- (lp[,k-1] - lp[,1])*R[,k-1]^2
            index <- which(abs(denominator) > 10^(-15))
            w[index, k-1] <- numerator[index]/denominator[index]
            w[which(w[,k-1] >= 1), k-1] <- 1
        }
    }
    phi <- sapply(1:N, function(i) runif(1, acos(sqrt(w[i,2]))/(2*pi), 1/2 - acos(sqrt(w[i,2]))/(2*pi)))
    index <- sample(1:N, round(N/2))
    phi[index] <- 1/2 + phi[index]
    C[,2] <- cos(2*pi*phi)
    return(list(S1 = S1, S2 = S2, C= C))
}

dq.small.reduced <- function(sample.list, j, mixture.param){
    d <- mixture.param$d
    lambda <- mixture.param$lambda/(mixture.param$T.obs/mixture.param$denominator[j])
    r <- length(which(mixture.param$lambda != 0))
    ## s <- mixture.param$s
    if(r <= 1 | r >= d){
        stop("this method is for 1 < r < d")
    }

    S1 <- sample.list$S1
    S2 <- sample.list$S2
    C <- sample.list$C
    N <- length(S1)
    lp <- matrix(rep(lambda[1:r], N), ncol = r, byrow = TRUE)
    lp <- lp/(S1 + S2)
    N <- length(S1)
    z <- matrix(NA, N, r)
    R <- matrix(NA, N, r)
    w <- matrix(1, N, r)
    ## get w[,r]
    numerator <- (1 - S1*lp[,1])
    denominator <- (lp[,r] - lp[,1])*S1
    index <- which(abs(denominator) > 10^(-15))
    w[index,r] <- numerator[index]/denominator[index]
    w[which(w[,r] >= 1), r] <- 1

    density <- dchisq(S1, r)*dchisq(S2, d-r)/pchisq((lambda[1] - 1)*S1, d - r)
    R[,r] <- sqrt(S1)
    if(r >= 3){
        for(k in r:3){
            z[,k] <- R[,k]*C[,k]
            R[,k-1] <- sqrt(R[,k]^2 - z[,k]^2)
            numerator <- (1 - rowSums(z[,k:r, drop = FALSE]^2*lp[,k:r, drop = FALSE]) - (S1 - rowSums(z[,k:r, drop = FALSE]^2))*lp[,1, drop = FALSE])
            denominator <- (lp[,k-1] - lp[,1])*R[,k-1]^2
            index <- which(abs(denominator) > 10^(-15))
            w[index, k-1] <- numerator[index]/denominator[index]
            w[which(w[,k-1] >= 1), k-1] <- 1
        }
        density <- density*apply(matrix(sapply(3:r, function(k) 1/(beta(1/2,(k-1)/2)*pbeta(w[,k],1/2, (k-1)/2))*(1-C[,k]^2)^((k-1)/2 -1)), nrow = N), 1, prod)
    }
    density <- density*1/(2*pi*sqrt(1-C[,2]^2))*pi/(pi - 2*acos(sqrt(w[,2])))
    return(density)
}

get.indicator.regular <- function(sample, j, mixture.param){
    get.T <- function(y){
        return(sum(drop(y%*%diag(mixture.param$lambda))*y))
    }
    return(as.numeric(apply(sample, 1, get.T) >= mixture.param$T.obs))

}

get.indicator.reduced <- function(sample.list, j, mixture.param){
    d <- mixture.param$d
    lambda <- mixture.param$lambda/mixture.param$T.obs
    r <- length(which(mixture.param$lambda != 0))
    get.T <- function(y){
        return(sum(drop(y%*%diag(mixture.param$lambda[1:r]))*y))
    }
    if(r <= 1 | r >= d){
        stop("this method is for 1 < r < d")
    }
    S1 <- sample.list$S1
    S2 <- sample.list$S2
    N <- length(S1)
    C <- sample.list$C
    z <- matrix(NA, N, r)
    R <- matrix(NA, N, r)
    w <- matrix(NA, N, r)
    R <- matrix(NA, N, r)
    R[,r] <- sqrt(S1)
    if(r >= 3){
        for(k in r:3){
            z[,k] <- R[,k]*C[,k]
            R[,k-1] <- sqrt(R[,k]^2 - z[,k]^2)
        }
    }
    z[,2] <- R[,2]*C[,2]
    z[,1] <- R[,2]*sqrt(1 - C[,2]^2)
    z <- z/sqrt(S1 + S2)
    return(as.numeric(apply(z, 1, get.T) >= mixture.param$T.obs))

}
