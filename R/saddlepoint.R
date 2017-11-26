get.exact.saddle <- function(x0, x1){
    m <- length(x0)
    n <- length(x1)
    N <- m + n
    p <- m/N
    q <- n/N
    w <- N*p*q/(N - 1)
    xc <- c(x0, x1)
    a <- (xc - mean(xc))/sqrt(sum((xc - mean(xc))^2))
    Tstar <- 1/sqrt(w)*sum(a[1:m])
    A <- combn(a, m)
    T <- 1/sqrt(w)*apply(A, 2, sum)
    pval <- length(which(abs(T) >= abs(Tstar)))/length(T)
    return(pval)
}


get.saddle <- function(x0, x1){
    m0 <- length(x0)
    m1 <- length(x1)
    if(min(x0) > max(x1) | max(x0) < min(x1)){
        if(m0 == m1)
            return(list(saddle1 = 2/choose(m0 + m1, m0)))
        else
            return(list(saddle1 = 1/choose(m0 + m1, m0)))
    }
    K <- function(x, p, q){
        a = p*exp(q*x)
        b = q*exp(-p*x)
        if(a == Inf){
            return(log(p) + q*x)
        }else if(b == Inf){
            return(log(q) - p*x)
        }else{
            return(log(a + b))
        }
    }

    K1 <- function(x,p ,q){
        if(x > 0){
            return(1/(p + q*exp(-(p+q)*x))*(p*q - q*p*exp(-(p+q)*x)))
        }else{
            return(1/(p*exp((p+q)*x) + q)*(p*q*exp((p+q)*x) - q*p))
        }
    }

    K2 <- function(x, p, q){
        if(x > 0){
            return((p^2*q*exp(-(p + q) * x) + p*q^2)/(p+ q*exp(-(p+q)*x)) - (p*q - p*q*exp(-(p+q)*x))^2/(p + q*exp(-(p+q)*x))^2)
        }else{
            return((p^2*q + p*q^2*exp((p + q) * x))/(p*exp((p+q)*x)+ q) - (p*q*exp((p+q)*x) - p*q)^2/(p*exp((p+q)*x) + q)^2)
        }
    }

    K3 <- function(x, p, q){
        ## (p*q^3*exp(q*x) - p^3*q*exp(-p*x))/(p*exp(q*x) + q*exp(-p*x)) - 3*(p^2*q*exp(-p*x) + p*q^2*exp(q*x))*(p*q*exp(q*x) - p*q*exp(-p*x))/(p*exp(q*x) + q*exp(-p*x))^2 + 2*(p*q*exp(q*x) - p*q*exp(-p*x))^3/(p*exp(q*x) + q*exp(-p*x))^3
        if(x > 0){
            (p*q^3 - p^3*q*exp(-(p+q)*x))/(p + q*exp(-(p+q)*x)) - 3*(p^2*q*exp(-(p+q)*x) + p*q^2)*(p*q - p*q*exp(-(p+q)*x))/(p + q*exp(-(p+q)*x))^2 + 2*(p*q - p*q*exp(-(p+q)*x))^3/(p + q*exp(-(p+q)*x))^3
        }else{
            (p*q^3*exp((p+q)*x) - p^3*q)/(p*exp((p+q)*x) + q) - 3*(p^2*q + p*q^2*exp((p+q)*x))*(p*q*exp((p+q)*x) - p*q)/(p*exp((p+q)*x) + q)^2 + 2*(p*q*exp((p+q)*x) - p*q)^3/(p*exp((p+q)*x) + q)^3
        }
    }


    mf <- function(a, u, w, p, q){
        N <- length(a)
        alpha <- solve.alpha(a, u, w, p, q)
        return(1/sqrt(w)*sum(a* sapply(1:N, function(i) K1(u*a[i]/sqrt(w) + alpha, p, q))))
    }


    solve.alpha <- function(a, u, w, p, q){
        k <- 100
        N <- length(a)
        tmp <- try(result <- uniroot(function(alpha) sum(sapply(1:N, function(i) K1(u*a[i]/sqrt(w) + alpha, p, q))), tol = .Machine$double.eps^0.5, interval = c(-k, k))$root, silent = TRUE)
        count <- 0
        while(class(tmp) == "try-error"&count < 10){
            count <- count + 1
            k <- k*5
            tmp <- try(result <- uniroot(function(alpha) sum(sapply(1:N, function(i) K1(u*a[i]/sqrt(w) + alpha, p, q))), tol = .Machine$double.eps^0.5, interval = c(-k, k))$root, silent = TRUE)
        }
        if(class(tmp) == "try-error"){
            return(NA)
        }
        return(result)
    }


    solve.u <- function(t, a, w, p, q){
        k <- 100
        tmp <- try(result <- uniroot(function(u) mf(a, u, w, p ,q) - t, tol = .Machine$double.eps^0.5, interval = c(-k, k))$root, silent = TRUE)
        count <- 0
        while(class(tmp) == "try-error"&count<10){
            count <- count + 1
            k <- k*2
            tmp <- try(result <- uniroot(function(u) mf(a, u, w, p, q) - t, tol = .Machine$double.eps^0.5, interval = c(-k, k))$root, silent = TRUE)
        }
        if(class(tmp) == "try-error"){
            return(NA)
        }
        return(result)
    }


    Q <- function(a, u, w, alpha, p ,q){
        x <- u*a/sqrt(w) + alpha
        N <- length(a)
        Kseq <- sapply(1:N, function(i) K(x[i], p, q))
        K2seq <- sapply(1:N, function(i) K2(x[i], p, q))
        return((1/(N*p*q)*sum(K2seq))^(-1/2)*exp(sum(Kseq)))
    }

    logQ <- function(a, u, w, alpha, p ,q){
        x <- u*a/sqrt(w) + alpha
        N <- length(a)
        Kseq <- sapply(1:N, function(i) K(x[i], p, q))
        K2seq <- sapply(1:N, function(i) K2(x[i], p, q))
        return(-1/2*log(1/(N*p*q)*sum(K2seq)) + sum(Kseq))
    }

    lognormtail <- function(x){
        log(1 - exp(-1.4 * x)) - log(x) - x^2/2 -1.04557
    }

    W1 <- function(v){
        dnorm(v)/(1-pnorm(v)) - v
    }


    W3 <- function(v){
        (v^2 - 1)*dnorm(v)/(1-pnorm(v)) - v^3
    }

    m <- length(x0)
    n <- length(x1)
    N <- m + n
    p <- m/N
    q <- n/N
    w <- N*p*q/(N - 1)
    xc <- c(x0, x1)
    a <- (xc - mean(xc))/sqrt(sum((xc - mean(xc))^2))

    Tstar <- 1/sqrt(w)*sum(a[1:m])
    ## print("yay")
    ## save(a, w, p, q, Tstar, file = "tmp.rda")
    if(Tstar < 0){
        return(get.saddle(x1, x0))
    }
    u <- solve.u(Tstar, a, w, p, q)
    alpha <- solve.alpha(a, u, w, p, q)
    if(is.na(u) | is.na(alpha)){
        return(list(saddle1 = NA, saddle2 = NA))
    }
    x <- u*a/sqrt(w) + alpha
    K1seq <- sapply(1:N, function(i) K1(x[i], p, q))
    K2seq <- sapply(1:N, function(i) K2(x[i], p, q))
    K3seq <- sapply(1:N, function(i) K3(x[i], p, q))
    m.n <- 1/sqrt(w)*sum(a*K1seq)
    sigmau2 <- 1/w*(sum(a^2*K2seq) - sum(a*K2seq)^2/sum(K2seq))
    A1 <- exp(logQ(a, u, w, alpha, p, q) -u*m.n + 1/2*u^2*sigmau2)*(1-pnorm(u*sqrt(sigmau2))) #saddle point estimation 1
    if(is.nan(A1)){
        A1 <- logQ(a, u, w, alpha, p, q) + (-u*m.n + 1/2*u^2*sigmau2) + lognormtail(u*sqrt(sigmau2))
        A1 <- exp(A1)
    }
    H <- sum(a*K2seq)/sum(K2seq)
    k1 <- (sum(a*K3seq)/sum(K2seq) - sum(K3seq)*sum(a*K2seq)/sum(K2seq)^2)/sqrt(sigmau2)
    k3 <- (sum(a^3*K3seq) - 3*sum(a^2*K3seq*H) + 3*sum(a*K3seq*H^2 - H^3))/sqrt(sigmau2)^3
    B1 <- A1*(1- 1/2*k1*W1(u*sqrt(sigmau2)) + 1/6*k3*W3(u*sqrt(sigmau2)))
    return(list(saddle1 = A1*2, saddle2 = B1*2))
}
