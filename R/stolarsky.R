get.stolarsky <- function(m0, m1, x, rho.hat, x0 = NULL, x1 = NULL, two.sided = TRUE, rel.tol = .Machine$double.eps^0.5, do.phat1 = TRUE, do.phat2 = TRUE, do.phat3 = TRUE, do.sd = FALSE){
    ## ' compute log of the surface area of the sphere S^d
    ## ' @param d R^(d + 1) euclidean space
    ## ' @return log of the surface area of the sphere S^d
    ## ' that is embedded in the R^n euclidean space
    lnsphereS <- function(d){
        log(2) + ((d+1)/2)*log(pi) - lgamma((d+1)/2)
    }

    ## ' compute the surface area of the sphere S^d
    ## ' @param d R^(d+1) euclidean space
    ## ' @return compute the surface area of the sphere S^d
    ## ' that is embedded in the R^d euclidean space
    nsphereS <- function(d){
        2*exp(((d+1)/2)*log(pi) - lgamma((d+1)/2))
    }

    ## ' compute phat from rho.hat
    ## ' @param rho.hat rho.hat
    ## ' @param d the intrinsic dimension of the sphere S^d
    ## ' @return phat
    get.phat <- function(rho.hat, d){
        get.phat.help <- function(rho.hat, d){
            if(rho.hat >= 0 & rho.hat <= 1){
                1/2*pbeta(1-rho.hat^2, d/2, 1/2)
            }else if(rho.hat >= -1 & rho.hat <0){
                1 - 1/2*pbeta(1-rho.hat^2, d/2, 1/2)
            }else if(rho.hat > 1){
                return(0)
            }else if(rho.hat < -1){
                return(1)
            }
        }
        sapply(rho.hat, get.phat.help, d=d)
    }


    ## ' compute the dimension of S^d that all points lie on
    ## ' @param m0 number of controls
    ## ' @param m1 number of treatments
    ## ' @return the dimension of the sphere that the
    ## ' centered and normalized x, y lied in
    get.d <- function(m0, m1){
        return(m0 + m1 - 2)
    }

    ## ' return inner product for normalized vectors with m0 + and m1 - with r transpositions
    get.ip <- function(m0, m1, r){
        return(1 - r*(1/m0 + 1/m1))
    }

    A.norm <- function(theta, d){
        1/2*pbeta(sin(theta)^2, d/2, 1/2)
    }

    ## center and normalize a
    normalize <- function(a){
        a <-  a - mean(a)
        a <- a/sqrt(sum(a^2))
        return(a)
    }

    get.Ep.step1 <- function(rho.hat, d, two.sided){
        if(two.sided){
            return(get.phat(rho.hat, d)*2)
        }else{
            return(get.phat(abs(rho.hat), d))
        }

    }


    V2 <- function(m0, m1, ur, t, two.sided){
        V2.help <- function(m0, m1, ur, t){
            if(t<0)
                print("t < 0 in V2")
            d <- m0 + m1 -2
            if(abs(ur - 1) <= .Machine$double.eps){
                return(get.phat(abs(t), d))
            }else if(abs(ur + 1) <= .Machine$double.eps){
                return(0)
            }else{
                if(t >= 0){
                    return(nsphereS(d-1)/nsphereS(d)*integrate(function(s) (1-s^2)^(d/2 - 1)*get.phat((t - s*ur)/sqrt(1-s^2)/sqrt(1-ur^2), d-1), t, 1, rel.tol = rel.tol)$value)
                }else{
                    return(nsphereS(d-1)/nsphereS(d)*integrate(function(s) (1-s^2)^(d/2 - 1)*(1-get.phat((t - s*ur)/sqrt(1-s^2)/sqrt(1-ur^2), d-1)), -1, t, rel.tol = rel.tol)$value)
                }
            }
        }
        if(two.sided){
            return(2*(V2.help(m0, m1, ur, abs(t)) + V2.help(m0, m1, -ur, abs(t))))
        }else{
            if(t > 0){
                return(V2.help(m0, m1, ur, t))
            }else{
                return(V2.help(m0, m1, ur, abs(t)))
                ## return(V2.help(m0, m1, ur, t))
            }

        }

    }

    get.Ep2.step1 <- function(rho.hat, m0, m1, two.sided){
        S <- 0
        for(r in 0:min(m0, m1)){
            ur <- 1 - r*(1/m0 + 1/m1)
            S <- S + choose(m0, r)*choose(m1, r)* V2(m0, m1, ur, rho.hat, two.sided)

        }
        S <- S/choose(m0 + m1, m0)
        return(S)
    }

    ## get.MSE <- function(m0, m1, phat0, two.sided){
    ##     d <- m0 + m1 - 2
    ##     rho.hat <- get.rho.hat(phat0, d)
    ##     A <- sapply(1:length(rho.hat), function(i) get.Ep2.step1(rho.hat[i], m0, m1, two.sided) - phat0[i]^2)
    ##     return(A)
    ## }


    P1 <- function(m0, m1, ur, rho.tilde, rho.hat, two.sided){
        P1.help <- function(m0, m1, ur, rho.tilde, rho.hat){
            d <- m0 + m1 -2
            if(abs(ur - 1)<= .Machine$double.eps | abs(ur + 1)<= .Machine$double.eps | abs(rho.tilde + 1)<= .Machine$double.eps|abs(rho.tilde + 1)<= .Machine$double.eps){
                return(as.numeric(rho.tilde*ur >= rho.hat))
            }else{
                return(get.phat((rho.hat - rho.tilde*ur)/sqrt(1-rho.tilde^2)/sqrt(1-ur^2), d-1))
            }
        }
        if(two.sided){
            ## print("P1.help(m0, m1, ur, rho.tilde, abs(rho.hat))")
            ## print(P1.help(m0, m1, ur, rho.tilde, abs(rho.hat)))
            ## print(P1.help(m0, m1, -ur, rho.tilde, abs(rho.hat)))
            return(P1.help(m0, m1, ur, rho.tilde, abs(rho.hat)) + P1.help(m0, m1, -ur, rho.tilde, abs(rho.hat)))
        }else{
            if(rho.hat > 0){
                return(P1.help(m0, m1, ur, rho.tilde, rho.hat))
            }else if(rho.hat <= 0){
                return(P1.help(m0, m1, -ur, rho.tilde, -rho.hat))
            }

        }
    }


    get.Ep.step2 <- function(m0, m1, rho.tilde, rho.hat, two.sided){
        N <- choose(m0 + m1, m0)
        Ep <- 0
        for(r in 0:min(m0, m1)){
            ur <- 1 - r*(1/m0 + 1/m1)
            Ep <-  Ep + choose(m0, r)*choose(m1, r)*P1(m0, m1, ur, rho.tilde, rho.hat, two.sided)
        }
        Ep <- Ep/N
        return(Ep)
    }

    P2 <- function(m0, m1, u1, u2, uk, rho.tilde, rho.hat, two.sided){
        P2.help <- function(m0, m1, u1, u2, uk, rho.tilde, rho.hat){
            d <- m0 + m1 -2
            if(u1 == u2 & u2 == uk & uk == 1){
                return(as.numeric(rho.tilde >= rho.hat))
            }else if(u1 == 1& u2 !=1  & uk !=1){
                return(as.numeric(rho.tilde >= rho.hat)*P1(m0, m1, u2, rho.tilde, rho.hat, two.sided = FALSE))
            }else if(u2 == 1& u1 != 1 & uk != 1){
                return(as.numeric(rho.tilde >= rho.hat)*P1(m0, m1, u1, rho.tilde, rho.hat, two.sided = FALSE))
            }else if(u1 == u2 & u2 != 1 & uk == 1){
                return(P1(m0, m1, u2, rho.tilde, rho.hat, two.sided = FALSE))
            }else{
                if(u1 == -1){
                    return(as.numeric(-rho.tilde >= rho.hat)*P1(m0, m1, u2, rho.tilde, rho.hat, two.sided = FALSE))
                }
                if(u2 == -1){
                    return(as.numeric(-rho.tilde >= rho.hat)*P1(m0, m1, u1, rho.tilde, rho.hat, two.sided = FALSE))
                }

                ukstar <- (uk - u1*u2)/sqrt(1 - u1^2)/sqrt(1 - u2^2)
                if(abs(rho.tilde + 1) <= .Machine$double.eps || abs(rho.tilde - 1) <= .Machine$double.eps){
                    return(as.numeric(rho.tilde*u1 >= rho.hat)*as.numeric(rho.tilde*u2 >= rho.hat))
                }
                rho1 <- (rho.hat - rho.tilde * u1)/sqrt(1 - rho.tilde^2)/sqrt(1 - u1^2)
                rho2 <- (rho.hat - rho.tilde * u2)/sqrt(1 - rho.tilde^2)/sqrt(1 - u2^2)
                if(abs(ukstar - 1) < 10^(-10)){
                    lower <- max(rho1, rho2, -1)
                    if(lower < 1){
                        return(integrate(function(t) (1-t^2)^(d/2-3/2), lower, 1, rel.tol= rel.tol)$value*nsphereS(d-2)/nsphereS(d-1))
                    }else{
                        return(0)
                    }
                }else if(abs(ukstar + 1) < 10^(-10)){
                    lower <- max(-1, rho1)
                    upper <- min(1, -rho2)
                    if(lower < upper){
                        return(integrate(function(t) (1-t^2)^(d/2-3/2), lower, upper, rel.tol= rel.tol)$value*nsphereS(d-2)/nsphereS(d-1))
                    }else{
                        return(0)
                    }
                }else{
                    lower <- max(-1, rho1)
                    if(lower < 1){
                        return(integrate(function(t) (1-t^2)^(d/2-3/2)*get.phat((rho2 - t*ukstar)/sqrt(1-t^2)/sqrt(1-ukstar^2), d-2), lower, 1, rel.tol= rel.tol)$value*nsphereS(d-2)/nsphereS(d-1))
                    }else{
                        return(0)
                    }
                }
            }
        }
        if(two.sided){
            ## print("values")
            ## print(c(P2.help(m0, m1, u1, u2, uk, rho.tilde, abs(rho.hat)), P2.help(m0, m1, -u1, u2, -uk, rho.tilde, abs(rho.hat)), P2.help(m0, m1, u1, -u2, -uk, rho.tilde, abs(rho.hat)), P2.help(m0, m1, -u1, -u2, uk, rho.tilde, abs(rho.hat))))
            return(P2.help(m0, m1, u1, u2, uk, rho.tilde, abs(rho.hat)) +
                       P2.help(m0, m1, -u1, u2, -uk, rho.tilde, abs(rho.hat)) +
                           P2.help(m0, m1, u1, -u2, -uk, rho.tilde, abs(rho.hat)) +
                               P2.help(m0, m1, -u1, -u2, uk, rho.tilde, abs(rho.hat)))
        }else{
            if(rho.hat > 0){
                return(P2.help(m0, m1, u1, u2, uk, rho.tilde, rho.hat))
            }else{
                return(P2.help(m0, m1, -u1, -u2, uk, rho.tilde, abs(rho.hat)))
            }

        }
    }

    ## test.P2 <- function(){
    ##     m0 <- 10
    ##     m1 <- 10
    ##     r1 <- 4
    ##     r2 <- 4
    ##     k <- 1
    ##     rho.tilde <- 0.8
    ##     rho.hat <- 0.8
    ##     print(P2(m0, m1, r1, r2, k, rho.tilde, rho.hat))
    ## }



    K4 <- function(m0, m1, rho.tilde, rho.hat, two.sided){
        s <- 0
        for(r1 in 1:min(m0, m1)){
            u1 <- 1 - r1*(1/m0 + 1/m1)
            for(r2 in 1:min(m0, m1)){
                u2 <- 1 - r2*(1/m0 + 1/m1)
                klow <- max(1, (r1 + r2 - 2*min(r1, r2)))
                kup <- min(r1+r2,m0,m1,m0+m1-r1-r2)
                if(klow <= kup){
                    for(k in klow:kup){
                        uk <- 1 - k*(1/m0 + 1/m1)
                        coef <- 0
                        delta1low <- max(max(0, r1+r2 - m0), r1 + r2 - k - min(r1, r2))
                        delta1up <- min(min(r1, r2), r1+r2 - k -max(0, r1 + r2 - m1))
                        if(delta1low <= delta1up){
                            for(delta1 in delta1low:delta1up){
                                delta2 <- r1 + r2 - k - delta1
                                coef <- coef + crdelta(m0, m1, r1, r2, delta1, delta2)
                            }

                            s <- s + coef*P2(m0, m1, u1, u2, uk, rho.tilde, rho.hat, two.sided)
                        }
                    }
                }
            }
        }
        return(s)
    }


    get.Ep2.step2 <- function(m0, m1, rho.tilde, rho.hat, two.sided){
        N <- choose(m0 + m1, m0)
        K1 <- P2(m0, m1, 1,1,1,rho.tilde, rho.hat, two.sided)
        K2 <- 0
        for(r2 in 1:min(m0, m1)){
            u2 <- 1 - r2*(1/m0 + 1/m1)
            K2 <- K2 + 2*choose(m0, r2)*choose(m1, r2)*P2(m0, m1, 1, u2, u2, rho.tilde, rho.hat, two.sided)
        }
        K3 <- 0
        for(r2 in 1:min(m0, m1)){
            u2 <- 1 - r2*(1/m0 + 1/m1)
            K3 <- K3 + choose(m0, r2)*choose(m1, r2)*P2(m0, m1, u2, u2, 1,rho.tilde, rho.hat, two.sided)
        }
        tmp <- c(K1, K2, K3, K4(m0, m1, rho.tilde, rho.hat, two.sided))
        ## print("tmp/N^2")
        ## print(tmp/N^2)
        Ep2 <- sum(tmp)/N^2
        return(Ep2)
    }




    crdelta <- function(m0, m1, r1, r2, delta1, delta2){
        if(delta1 > m0 | delta2 > m1 | r1 > m0 | r1 > m1 |r1  +r2 >m0 + delta1|r1 + r2 > m1  +delta2){
            print("m0, m1, r1, r2, delta1, delta2")
            print(c(m0, m1, r1, r2, delta1, delta2))
            stop("binomial coefficient is wrong")
        }
        choose(m0, delta1)*choose(m1, delta2)*choose(m0-delta1, r1 - delta1)*choose(m1 - delta2, r1 - delta2)*choose(m0 - r1, r2 - delta1)*choose(m1 - r1, r2 - delta2)
    }


    if(missing(m0) & missing(m1) & missing(rho.hat) &!missing(x0) &!missing(x1)){
        m0 <- length(x0)
        m1 <- length(x1)
        x <- c(x0, x1)
        y <- c(rep(0, m0), rep(1, m1))
        rho.hat <- cor(x, y)
    }else if(missing(m0) | missing(m1) | missing(rho.hat)){
        stop("Either input x0, x1, or input m0, m1, rho.hat")
    }
    d <- m0 + m1 -2


    if(do.phat1){
        phat1 <- get.Ep.step1(abs(rho.hat), d, two.sided)
        if(do.sd){
            sd1 <- sqrt(get.Ep2.step1(rho.hat, m0, m1, two.sided) - phat1^2)
        }
    }
    if(do.phat2){
        phat2 <- get.Ep.step2(m0, m1, rho.hat, rho.hat, two.sided)
        if(do.sd){
            sd2 <- sqrt(get.Ep2.step2(m0, m1, rho.hat, rho.hat, two.sided) - phat2^2)
        }
    }
    if(do.phat3){
        y.min.1 <- rep(0, m0 + m1)
        y.min.2 <- rep(0, m0 + m1)
        y.min.1[order(x, decreasing = TRUE)[1:m1]] <- 1
        y.min.2[order(-x, decreasing = TRUE)[1:m1]] <- 1
        rho.tilde <- ifelse(abs(cor(x, y.min.1)) >=abs(cor(x, y.min.2)), cor(x, y.min.1), cor(x, y.min.2))
        ## y.min <- rep(0, m0 + m1)
        ## y.min[order(x, decreasing = TRUE)[1:m1]] <- 1
        ## y.min <- normalize(y.min)
        ## rho.tilde <- sum(x*y.min)

        phat3 <- get.Ep.step2(m0, m1, rho.tilde, rho.hat, two.sided)
        if(do.sd){
            sd3 <- sqrt(get.Ep2.step2(m0, m1, rho.tilde, rho.hat, two.sided) - phat3^2)
        }

    }
    result <- list(two.sided = two.sided, rho.hat = rho.hat)
    if(do.phat1){
        result[["phat1"]] <- phat1
        if(do.sd){
            result[["sd1"]] <- sd1
        }
    }
    if(do.phat2){
        result[["phat2"]] <- phat2
        if(do.sd){
            result[["sd2"]] <- sd2
        }
    }
    if(do.phat3){
        result[["phat3"]] <- phat3
        if(do.sd){
            result[["sd3"]] <- sd3
        }
    }
    result[["rho.hat"]] <- rho.hat
    return(result)
}
