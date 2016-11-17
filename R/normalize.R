get.rho.hat <- function(phat, d){
    get.rho.hat.help <- function(phat, d){
        if(phat > 1 | phat < 0){
            stop("phat value out of range")
        }
        sqrt(1 - qbeta(phat, d/2, 1/2))
    }
    sapply(phat, get.rho.hat.help, d=d)
}
