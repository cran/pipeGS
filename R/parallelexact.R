get.exact <- function(x0, x1, statistic = c("linear", "quadratic"), parallel = FALSE, estimate.time = FALSE, mc.cores = 100){
    get.exact.par.help <- function(m0, m1, k, N.par, T.obs, statfunc = "get.T"){
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
            T[i] <- eval(parse(text = statfunc))(y)
            comb <- next.comb(m0, m0 + m1, comb)
        }
        if(T.obs < 0){
            return(list(N = end - start + 1, sig.count.1 = length(which(T <= T.obs)), sig.count.2 = length(which(abs(T) >= abs(T.obs)))))
        }else{
            return(list(N = end - start + 1, sig.count.1 = length(which(T >= T.obs)), sig.count.2 = length(which(abs(T) >= abs(T.obs)))))
        }
    }

    if(!all(dim(x0)[2] == dim(x1)[2])){
        return("can't do test with inputs with different columns")
    }else if(is.null(dim(x0))){
        m0 <- length(x0)
        m1 <- length(x1)
        x <- c(x0, x1)
        x <- normalize(x)
        if(statistic == "linear"){
            get.T <- function(y){
                return(cor(x, y))
            }
        }else if(statistic == "quadratic"){
            get.T <- function(y){
                return(sum(x*y)^2)
            }
        }
    }else{
        m0 <- dim(x0)[1]
        m1 <- dim(x1)[1]
        p <- dim(x0)[2]
        x <- rbind(x0, x1)
        x <- apply(x, 2, normalize)
        if(statistic == "linear"){
            x.vec <- rowSums(x)
            get.T <- function(y){
                return(cor(x.vec, y))
            }
        }else if(statistic == "quadratic"){
            Sigma <- x%*%t(x)
            get.T <- function(y){
                return(sum(drop(y%*%Sigma)*y))
            }
        }
    }
    plus <- sqrt(m0/(m0 + m1)/m1)
    minus <- - sqrt(m1/(m0 + m1)/m0)
    y.obs <- c(rep(minus, m0), rep(plus, m1))
    T.obs <- get.T(y.obs)
    if(estimate.time){
        print("estimated running time in hours")
        print(system.time(mclapply(1:10, function(i) get.exact.par.help(m0, m1, i, 10^4, T.obs), mc.cores = mc.cores))*1000/60/60)
    }
    if(parallel){
        a <- mclapply(1:mc.cores, function(i) get.exact.par.help(m0, m1, i, mc.cores, T.obs), mc.cores = mc.cores)
    }else{
        a <- lapply(1:mc.cores, function(i) get.exact.par.help(m0, m1, i, mc.cores, T.obs))
    }
    aa <- colSums(matrix(unlist(a), ncol = 3, byrow = TRUE))
    p.value.onesided <- aa[2] / aa[1]
    p.value.twosided <- aa[3]/ aa[1]
    return(list(p.star = p.value.twosided))
    ## return(list(pval.onesided = p.value.onesided, pval.twosided = p.value.twosided, N = choose(m0 + m1, m0), N.sig.1sided = aa[2], N.sig.2sided = aa[3]))
}




get.mc <- function(x0, x1, statistic = c("linear", "quadratic"), N.mc = 10^6, parallel = FALSE, estimate.time = FALSE, mc.cores = 100){
    if(!all(dim(x0)[2] == dim(x1)[2])){
        return("can't do test with inputs with different columns")
    }else if(is.null(dim(x0))){
        m0 <- length(x0)
        m1 <- length(x1)
        x <- c(x0, x1)
        x <- normalize(x)
        if(statistic == "linear"){
            get.T <- function(y){
                return(cor(x, y))
            }
        }else if(statistic == "quadratic"){
            get.T <- function(y){
                return(sum(x*y)^2)
            }
        }
    }else{
        m0 <- dim(x0)[1]
        m1 <- dim(x1)[1]
        p <- dim(x0)[2]
        x <- rbind(x0, x1)
        x <- apply(x, 2, normalize)
        if(statistic == "linear"){
            x.vec <- rowSums(x)
            get.T <- function(y){
                return(cor(x.vec, y))
            }
        }else if(statistic == "quadratic"){
            Sigma <- x%*%t(x)
            get.T <- function(y){
                return(sum(drop(y%*%Sigma)*y))
            }
        }
    }
    plus <- sqrt(m0/(m0 + m1)/m1)
    minus <- - sqrt(m1/(m0 + m1)/m0)
    y.obs <- c(rep(minus, m0), rep(plus, m1))
    T.obs <- get.T(y.obs)
    get.T.mc <- function(seed){
        a <- sample(1:(m0 + m1), m0)
        y <- rep(minus, m0 + m1)
        y[a] <- plus
        T <- get.T(y)
        return(T)
    }
    s <- sample.int(1e7, N.mc)
    T <- unlist(mclapply(1:length(s), function(i) get.T.mc(s[i]), mc.cores = 4))
    pval <- length(which(abs(T) >= abs(T.obs)))/length(T)
    return(pval)
}
