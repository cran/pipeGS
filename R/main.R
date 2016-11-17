#' Get the pvalue for the linear statistic
#' @export
#' @param x binary vector of treatment assignment
#' @param y gene expression measurement matrix
#' @param method method to estimate pvalue for the linear statistic
#' @param N.mc number of mc samples is method == "mc"
#' @param mc.cores number of cores to use in mcapply
#' @param N.level number of samples in each level if method == "nested"
#' @param q progression quantile if method == "nested"
#' @param B number of burn-ins to use if method == "nested"
#' @param do.sd logical to indicate whether calculating sd or not
#' @return a list containing phat from different methods
#' @examples
#' x <- c(rep(0, 4), rep(1, 4))
#' y <- c(rnorm(4, 0, 1), rnorm(4, 2 ,1))
#' get.phat.linear(x, y, method = "saddlepoint")
#' get.phat.linear(x, y, method = "phat1")
#' get.phat.linear(x, y, method = "phat2")
#' get.phat.linear(x, y, method = "phat3")
#' get.phat.linear(x, y, method = "nested")
#' get.phat.linear(x, y, method = "mc")
#' get.phat.linear(x, y, method = "exact")
#' get.phat.linear(x, y, method = "saddlepoint", do.sd = TRUE)
#' get.phat.linear(x, y, method = "phat1", do.sd = TRUE)
#' get.phat.linear(x, y, method = "phat2", do.sd = TRUE)
#' get.phat.linear(x, y, method = "phat3", do.sd = TRUE)
#' get.phat.linear(x, y, method = "nested", do.sd = TRUE)
#' get.phat.linear(x, y, method = "mc", do.sd = TRUE)
#' get.phat.linear(x, y, method = "exact", do.sd = TRUE)
get.phat.linear <- function(x, y, method = c("saddlepoint", "phat1","phat2", "phat3", "nested", "mc", "exact"), N.mc = 10^3,mc.cores = 1, N.level = 1000, q = 0.2, B = 5*length(x), do.sd = FALSE){
    if(!all(unique(x) %in% c(0, 1))){
        stop("x needs to be binary 0/1 vector")
    }
    m0.ind <- which(x == 0)
    m0 <- length(m0.ind)
    m1.ind <- which(x == 1)
    m1 <- length(m1.ind)
    if(is.null(dim(y))){
        y <- matrix(y, nrow = length(y))
    }
    if(dim(y)[1] != m0 + m1){
        stop("input y has different number of rows than length of x")
    }
    y <- apply(y, 2, normalize)
    y <- apply(y, 1, sum)
    y.norm <- normalize(y)
    x.norm <- normalize(x)
    x0 <- y.norm[m0.ind]
    x1 <- y.norm[m1.ind]
    rho.hat <- sum(x.norm*y.norm)
    method <- match.arg(as.character(method),
                        c("saddlepoint", "phat1","phat2", "phat3", "nested", "mc", "exact"))

    if(method == "saddlepoint"){
        result <- get.saddle(x0, x1)
    }else if(method == "phat1"){
        result <- get.stolarsky(m0, m1, y.norm, rho.hat, two.sided = TRUE, rel.tol = .Machine$double.eps^0.5, do.phat1 = TRUE, do.phat2 = FALSE, do.phat3 = FALSE, do.sd = do.sd)
    }else if(method == "phat2"){
        result <- get.stolarsky(m0, m1, y.norm, rho.hat, two.sided = TRUE, rel.tol = .Machine$double.eps^0.5, do.phat1 = FALSE, do.phat2 = TRUE, do.phat3 = FALSE, do.sd = do.sd)
    }else if(method == "phat3"){
        result <- get.stolarsky(m0, m1, y.norm, rho.hat, two.sided = TRUE, rel.tol = .Machine$double.eps^0.5, do.phat1 = FALSE, do.phat2 = FALSE, do.phat3 = TRUE, do.sd = do.sd)
    }else if(method == "nested"){
        result <- nested.estimation.linear(N.level, q, B, list(x0 = x0, x1 = x1), do.sd)
    }else if(method == "mc"){
        result <- plain.mc.linear(N.mc, list(x0 = x0, x1 = x1), mc.cores)
    }else if(method == "exact"){
        result <- get.exact(x0, x1, statistic = c("linear"), parallel = TRUE, estimate.time = FALSE, mc.cores = mc.cores)
    }
    return(result)
}

#' Get the pvalue for the quadratic statistic
#' @export
#' @param x binary vector of treatment assignment
#' @param y gene expression measurement matrix
#' @param method method to estimate pvalue for the linear statistic
#' @param N.mc number of mc samples is method == "mc"
#' @param N.is number of mc samples is method == "is"
#' @param N.level number of samples in each level if method == "nested"
#' @param q progression quantile if method == "nested"
#' @param B number of burn-ins to use if method == "nested"
#' @param mc.cores number of cores to use in mcapply
#' @param do.sd logical to indicate whether calculating sd or not
#' @return a list containing phat from different methods
#' @examples
#' x <- c(rep(0, 4), rep(1, 4))
#' y <- matrix(c(rnorm(4*5, 0, 1), rnorm(4*5, 2 ,1)), nrow = 8)
#' get.phat.quadratic(x, y, method = "is")
#' get.phat.quadratic(x, y, method = "nested")
#' get.phat.quadratic(x, y, method = "mc")
#' get.phat.quadratic(x, y, method = "exact")
#' get.phat.quadratic(x, y, method = "is", do.sd = TRUE)
#' get.phat.quadratic(x, y, method = "nested", do.sd = TRUE)
#' get.phat.quadratic(x, y, method = "mc", do.sd = TRUE)
#' get.phat.quadratic(x, y, method = "exact", do.sd = TRUE)
get.phat.quadratic <- function(x, y, method = c("is", "nested", "mc", "exact"), N.mc = 10^3, N.is = 10^3, N.level = 1000, q = 0.2, B = 5*length(x), mc.cores = 1, do.sd = FALSE){
    if(!all(unique(x) %in% c(0, 1))){
        stop("x needs to be binary 0/1 vector")
    }
    m0.ind <- which(x == 0)
    m0 <- length(m0.ind)
    m1.ind <- which(x == 1)
    m1 <- length(m1.ind)
    mixture.param <- list(m0 = m0, m1 = m1, p = m0 + m1, d = m0 + m1)
    if(is.null(dim(y))){
        y <- matrix(y, nrow = length(y))
    }
    if(dim(y)[1] != m0 + m1){
        stop("input y has different number of rows than length of x")
    }
    y <- apply(y, 2, normalize)
    G <- dim(y)[2]
    Sigma <- y%*%t(y)
    tmp <- svd(Sigma)
    lambda <- tmp$d
    eta <- tmp$u
    result <- list()
    plus <- sqrt(m0/(m0 + m1)/m1)
    minus <- - sqrt(m1/(m0 + m1)/m0)
    get.T <- function(x){
        return(sum(drop(x%*%Sigma)*x))
    }
    x.obs <- c(rep(minus, m0), rep(plus, m1))
    T.obs <- get.T(x.obs)
    lambda[lambda < 10^(-10)] <- 0
    mixture.param$lambda <- lambda
    mixture.param$eta <- eta
    mixture.param$Sigma <- Sigma
    mixture.param$T.obs <- T.obs
    mixture.param$denominator <- 1

    method <- match.arg(as.character(method),
                        c("is", "nested", "mc", "exact"))
    if(method == "is"){
        if(G < m0 + m1){
            rp = rp.small.reduced; dp = dp.small.reduced; rq = rq.small.reduced; dq = dq.small.reduced; get.indicator = get.indicator.reduced
        }else{
            rp = rp.small; dp = dp.small; rq = rq.small; dq = dq.small; get.indicator = get.indicator.regular
        }
        a <- rq(N.is, 1, mixture.param)
        phat.is <- mean(get.indicator(a, 1, mixture.param)*dp(a,mixture.param)/dq(a, 1, mixture.param), na.rm = TRUE)
        result <- list(phat.is = phat.is)
        if(do.sd){
            sd.is <- sd(get.indicator(a, 1, mixture.param)*dp(a,mixture.param)/dq(a, 1, mixture.param), na.rm = TRUE)/sqrt(N.is)
            result[["sd.is"]] <- sd.is
        }

    }else if(method == "nested"){
        result <- nested.estimation(N.level, q, B, mixture.param, do.sd = TRUE)
    }else if(method == "exact"){
        result <- get.exact(y[m0.ind, ], y[m1.ind, ], statistic = c("quadratic"), parallel = TRUE, estimate.time = FALSE, mc.cores = mc.cores)
    }else if(method == "mc"){
        p.mc <- plain.mc.quadratic(N.mc, mixture.param, mc.cores)
        result <- list(p.mc = p.mc)
    }
    return(result)
}
