uniform.mcmc.linear <- function(phi, x0, x1, B, N){
    m0 <- length(x0)
    m1 <- length(x1)
    plus <- sqrt(m0/(m0 + m1)/m1)
    minus <- - sqrt(m1/(m0 + m1)/m0)
    get.T <- function(x0, x1){
        return(sum(x1)*plus + sum(x0)*minus)
    }
    delta <- plus - minus
    x0old <- x0
    x1old <- x1
    Told <- get.T(x0old, x1old)

    update.T <- function(Told, x0old, x1old, k1, k2){
        Tnew <- Told + delta*x0old[k1] - delta*x1old[k2]
        return(Tnew)
    }
    for(k in 0:B){
        k1 <- sample(1:m0, 1)
        k2 <- sample(1:m1, 1)
        Tnew <- update.T(Told, x0old, x1old, k1, k2)
        if(abs(Tnew) >= abs(phi)){
            tmp <- x0old[k1]
            x0old[k1] <- x1old[k2]
            x1old[k2] <- tmp
            Told <- Tnew
        }
    }
    A <- matrix(0, N, m0 + m1)
    A[1,] <- c(x0old, x1old)
    if(1 <= N - 1){
        for(i in 1:(N-1)){
            k1 <- sample(1:m0, 1)
            k2 <- sample(1:m1, 1)
            Tnew <- update.T(Told, x0old, x1old, k1, k2)
            if(abs(Tnew) >= abs(phi)){
                tmp <- x0old[k1]
                x0old[k1] <- x1old[k2]
                x1old[k2] <- tmp
                Told <- Tnew
                A[i+1, ] <- c(x0old, x1old)
            }else{
                A[i+1, ] <- A[i,]
            }
        }
    }
    return(A)
}


nested.estimation.linear <- function(n, q = 0.2, B = 10, lineardata, do.sd = TRUE){
    x0 <- lineardata$x0
    x1 <- lineardata$x1
    m0 <- length(x0)
    m1 <- length(x1)
    plus <- sqrt(m0/(m0 + m1)/m1)
    minus <- - sqrt(m1/(m0 + m1)/m0)
    get.T <- function(x0, x1){
        return(sum(x1)*plus + sum(x0)*minus)
    }
    lineardata$T.obs <- get.T(x0, x1)
    xbind <- c(x0, x1)
    xbind <- normalize(xbind)
    x0 <- xbind[1:m0]
    x1 <- xbind[(m0+1):(m0+m1)]
    A.outer <- matrix(NA, n, m0 + m1)
    A.outer <- t(apply(A.outer, 1, function(x) {xbind[sample(1:(m0+m1), m0+m1)]}))
    Ts <- sapply(1:dim(A.outer)[1], function(i) get.T(A.outer[i,1:m0], A.outer[i, (m0+1):(m0+m1)]))
    Tq <- quantile(abs(Ts), probs = 1 - q)
    if(abs(Tq) >= abs(lineardata$T.obs) - .Machine$double.eps){
        P.star <- sum(abs(Ts) >= abs(lineardata$T.obs) - .Machine$double.eps)/length(Ts)
        N <- n
        N.burnin = 0
        l = 1
        ns <- n
        Ks <- 0
        result <- list(phat = P.star, N = sum(ns), N.burnin = sum(Ks)*B, L = l, ns = ns, Ks = Ks)
        if(do.sd){
            sigma.hat <- sqrt(P.star*(1-P.star)/length(Ts))
            result[["sigma.hat"]] <- sigma.hat
        }

        if(!is.null(lineardata$exact)){
            result$exact <- lineardata$exact
        }
        return(result)
    }
    indices <- which(abs(Ts) >= abs(Tq))
    ns <- nrow(A.outer)
    ss <- 0
    Ks <- length(indices)
    Ps <- Ks[1]/ns[1]
    deltas <- sqrt((1 - Ps[1])/(Ps[1]*ns[1]))
    l <- 1
    while(abs(Tq) < abs(lineardata$T.obs) - .Machine$double.eps){
        l <- l + 1
        ss <- c(ss, ceiling(n/Ks[l-1]))
        A.outer.list <- lapply(1:Ks[l-1], function(j) {uniform.mcmc.linear(Tq, A.outer[indices[j], 1:m0,drop = FALSE], A.outer[indices[j], (m0+1):(m0+m1),drop = FALSE],B, ss[l])})
        A.outer <- do.call(rbind, A.outer.list)
        Ts <- sapply(1:dim(A.outer)[1], function(i) get.T(A.outer[i,1:m0], A.outer[i, (m0+1):(m0+m1)]))
        Ts.mat <- matrix(Ts, nrow = ss[l], ncol = Ks[l-1])
        Tq <- quantile(abs(Ts), probs = 1 - q)
        if(Tq > abs(lineardata$T.obs)){
            Tq <- abs(lineardata$T.obs)
        }
        I.mat <- matrix(as.numeric(abs(Ts.mat) >= abs(Tq)), nrow = ss[l], ncol = Ks[l-1])
        indices <- which(abs(Ts) >= abs(Tq))
        Ks <- c(Ks, length(indices))
        ns <- c(ns, nrow(A.outer))
        Ps <- c(Ps, Ks[l]/ns[l])
        if(do.sd){
            R0 <- Ps[l]*(1-Ps[l])
            Rs <- rep(0, ss[l] - 1)
            for(j in 1:(ss[l] - 1)){
                Rs[j] <- 1/(ns[l] - j*Ks[l-1])*sum(apply(I.mat, 2, function(x) sum(x[1:(ss[l] - j)]*x[(1+j):ss[l]]))) - Ps[l]^2
            }
            rhos <- Rs/R0
            gamma <- 2*sum((1 - (1:(ss[l]-1))/ss[l])*rhos)
            if(Ps[l] == 1){
                delta <- 0
            }else{
                delta <- sqrt((1-Ps[l])/(Ps[l]*ns[l])*(1+gamma))
            }
            deltas <- c(deltas, delta)
        }
    }
    P.star <- prod(Ps)
    result <- list(phat = P.star, N = sum(ns), N.burnin = sum(Ks)*B, L = l, ns = ns, Ks = Ks)
    if(do.sd){
        delta <- sqrt(sum(deltas^2))
        sigma.hat <- delta*P.star
        result[["sigma.hat"]] <- sigma.hat
    }
    if(!is.null(lineardata$exact)){
        result$exact <- lineardata$exact
    }
    return(result)
}



uniform.mcmc <- function(phi, x0, B, N, mixture.param){
    get.T <- function(y, Sigma){
        return(sum(drop(y%*%Sigma)*y))
    }
    m0 <- mixture.param$m0
    m1 <- mixture.param$m1
    plus <- sqrt(m0/(m0 + m1)/m1)
    minus <- - sqrt(m1/(m0 + m1)/m0)
    delta <- plus - minus
    xold <- x0
    Told <- get.T(xold, mixture.param$Sigma)
    minus.set <- which(xold == minus)
    plus.set <- which(xold == plus)
    if(length(minus.set) == 0){
        print("paste(plus, minus)")
        print(paste(plus, minus))
        print("xold")
        print(xold)
        print("minus.set")
        print(minus.set)
        print("plus.set")
        print(plus.set)
    }

    update.T <- function(Told, xold, Sigma, k1, k2){
        Tnew <- Told + 2*sum(xold*(Sigma[,k1]-Sigma[,k2]))*delta + (Sigma[k2, k2] - Sigma[k1, k2] -
                                                                        Sigma[k2, k1] + Sigma[k1, k1])*delta^2
        if(is.na(Tnew)){
            ## print(paste0("Told", Told))
            ## print(paste0("delta", delta))
            ## print("c(k1, k2)")
            ## print(c(k1, k2))
            ## print("sum(xold*(Sigma[,k1]-Sigma[,k2])")
            ## print(sum(xold*(Sigma[,k1]-Sigma[,k2])))
            ## print("(Sigma[k2, k2] - Sigma[k1, k2] - Sigma[k2, k1] + Sigma[k1, k1])")
            ## print((Sigma[k2, k2] - Sigma[k1, k2] - Sigma[k2, k1] + Sigma[k1, k1]))
        }
        return(Tnew)
    }
    for(k in 0:B){
        k1 <- sample(1:m0, 1)
        k2 <- sample(1:m1, 1)
        if(is.na(minus.set[k1])){
            print("minus.set")
            print(minus.set)
            print("length(minus.set)")
            print(length(minus.set))
            print("k1")
            print(k1)
        }
        Tnew <- update.T(Told, xold, mixture.param$Sigma, minus.set[k1], plus.set[k2])
        if(Tnew >= phi){
            xold[minus.set[k1]] <- plus
            xold[plus.set[k2]] <- minus
            tmp <- plus.set[k2]
            plus.set[k2] <- minus.set[k1]
            minus.set[k1] <- tmp
            Told <- Tnew
        }
    }
    A <- matrix(0, N, length(xold))
    A[1,] <- xold
    if(1 <= N - 1){
        for(i in 1:(N-1)){
            k1 <- sample(1:m0, 1)
            k2 <- sample(1:m1, 1)
            Tnew <- update.T(Told, A[i,], mixture.param$Sigma, minus.set[k1], plus.set[k2])
            if(Tnew >= phi){
                A[i+1, ] <- A[i, ]
                A[i+1, minus.set[k1]] <- plus
                A[i+1, plus.set[k2]] <- minus
                tmp <- plus.set[k2]
                plus.set[k2] <- minus.set[k1]
                minus.set[k1] <- tmp
                Told <- Tnew
            }else{
                A[i+1, ] <- A[i,]
            }
        }
    }
    return(A)
}



nested.estimation <- function(n, q = 0.2, B = 10, mixture.param, do.sd = TRUE){
    m0 <- mixture.param$m0
    m1 <- mixture.param$m1
    plus <- sqrt(m0/(m0 + m1)/m1)
    minus <- - sqrt(m1/(m0 + m1)/m0)
    A.outer <- matrix(minus, n, m0 + m1)
    A.outer <- t(apply(A.outer, 1, function(x) {x[sample(1:(m0+m1), m1)] <- plus; x}))
    get.T <- function(y, Sigma){
        return(sum(drop(y%*%Sigma)*y))
    }
    Ts <- apply(A.outer, 1, get.T, mixture.param$Sigma)
    Tq <- quantile(Ts, probs = 1 - q)
    if(Tq >= mixture.param$T.obs - .Machine$double.eps){
        P.star <- sum(Ts >= mixture.param$T.obs - .Machine$double.eps)/length(Ts)
        N <- n
        N.burnin = 0
        l = 1
        ns <- n
        Ks <- 0
        result <- list(phat = P.star, N = sum(ns), N.burnin = sum(Ks)*B, L = l, ns = ns, Ks = Ks)
        if(do.sd){
            sigma.hat <- sqrt(P.star*(1-P.star)/length(Ts))
            result[["sigma.hat"]] <- sigma.hat
        }
        if(!is.null(mixture.param$exact)){
            result$exact <- mixture.param$exact
        }
        return(result)
    }
    indices <- which(Ts >= Tq)
    ns <- nrow(A.outer)
    ss <- 0
    Ks <- length(indices)
    Ps <- Ks[1]/ns[1]
    deltas <- sqrt((1 - Ps[1])/(Ps[1]*ns[1]))
    l <- 1
    while(Tq < mixture.param$T.obs - .Machine$double.eps){
        l <- l + 1
        ss <- c(ss, ceiling(n/Ks[l-1]))
        A.outer.list <- lapply(1:Ks[l-1], function(j) {uniform.mcmc(Tq, A.outer[indices[j], ,drop = FALSE], B, ss[l], mixture.param)})
        A.outer <- do.call(rbind, A.outer.list)
        Ts <- apply(A.outer, 1, get.T, mixture.param$Sigma)
        Ts.mat <- matrix(Ts, nrow = ss[l], ncol = Ks[l-1])
        Tq <- quantile(Ts, probs = 1 - q)
        if(Tq > mixture.param$T.obs){
            Tq <- mixture.param$T.obs
        }
        I.mat <- matrix(as.numeric(Ts.mat >= Tq), nrow = ss[l], ncol = Ks[l-1])
        indices <- which(Ts >= Tq)
        Ks <- c(Ks, length(indices))
        ns <- c(ns, nrow(A.outer))
        Ps <- c(Ps, Ks[l]/ns[l])
        if(do.sd){
            R0 <- Ps[l]*(1-Ps[l])
            Rs <- rep(0, ss[l] - 1)
            for(j in 1:(ss[l] - 1)){
                Rs[j] <- 1/(ns[l] - j*Ks[l-1])*sum(apply(I.mat, 2, function(x) sum(x[1:(ss[l] - j)]*x[(1+j):ss[l]]))) - Ps[l]^2
            }
            rhos <- Rs/R0
            gamma <- 2*sum((1 - (1:(ss[l]-1))/ss[l])*rhos)
            ## sigma <- sqrt(Ps[l]*(1-Ps[l])/ns[l]*(1+gamma))
            if(Ps[l] == 1){
                delta <- 0
            }else{
                delta <- sqrt((1-Ps[l])/(Ps[l]*ns[l])*(1+gamma))
            }
            deltas <- c(deltas, delta)
        }

    }
    P.star <- prod(Ps)

    result <- list(phat = P.star, N = sum(ns), N.burnin = sum(Ks)*B, L = l, ns = ns, Ks = Ks)
    if(do.sd){
        delta <- sqrt(sum(deltas^2))
        sigma.hat <- delta*P.star
        result[["sigma.hat"]] <- sigma.hat
    }
    if(!is.null(mixture.param$exact)){
        result$exact <- mixture.param$exact
    }
    return(result)
}
