# Main Author for this function:  Paul Blanche from previous code from the bpcp package.

# Description :
# This function gives NPMLE of the survival function
# under the constraint that the survival function is 'pstar'
# at 'tstar'.
#
# This is a slightly modified version of the funtion KMconstrainedfrom
# the bpcp package (bpcp_1.3.4, acessed in july 2016),
# We included some stuff of the kmgw.calc function of the same packge

KMconstrained <- function (tstar, pstar,  data) {
    # {{{ data management steps (just fo rme)
    data <- data[order(data$time),]
    time <- data$time
    status <- data$status
    # }}}
    # {{{ From kmgw.calc function packge bpcp
    N <- length(time)
    tab <- table(time, status)
    dstat <- dimnames(tab)$status
    if (length(dstat) == 2) {
        di <- tab[, 2]      
        ci <- tab[, 1]
    }
    else if (length(dstat) == 1) {
        if (dstat == "1") {
            di <- tab[, 1]
            ci <- rep(0, length(tab))
        }
        else if (dstat == "0") {
            ci <- tab[, 1]
            di <- rep(0, length(tab))
        }
    }
    else stop("status should be 0 or 1")
    y <- as.numeric(dimnames(tab)[[1]])
    k <- length(y)
    ni <- c(N, N - cumsum(ci)[-k] - cumsum(di)[-k])
    names(ni) <- names(di)
    ni <- as.numeric(ni)
    di <- as.numeric(di)
    ci <- as.numeric(ci)
    KM <- cumprod((ni - di)/ni)
    gw <- KM^2 * cumsum(di/(ni * (ni - di)))
    gw[KM == 0] <- 0
    I <- di > 0
    x <- list(time = y[I], # times (uncensored observations only)
              di = di[I],  # numbers of observed event at the times
              ni = ni[I],  # numbers of subjects at risk at the times
              KM = KM[I],  # KM estimates at the times
              gw = gw[I]   # Greenwood estimate of the variance of KM estimates at the times
              )
    # }}}
    # {{{ From KMconstrained of bpcp package
    II <- x$time <= tstar
    if (any(II)) {
        nj <- x$ni[II]
        dj <- x$di[II]
        rootfunc <- function(lambda) {
            nl <- length(lambda)
            out <- rep(NA, nl)
            for (i in 1:nl) {
                out[i] <- prod((nj + lambda[i] - dj)/(nj + lambda[i])) - 
                    pstar
            }
            out
        }
        lambda <- uniroot(rootfunc, c(-min(nj - dj), 10^3 * nj[1]))$root
        ## browser()
        pbar <- (nj + lambda - dj)/(nj + lambda)
        Sbar <- cumprod(pbar)
        out <- list(time = x$time[II], Sc = Sbar, x=x, tstar=tstar, pstar=pstar)
    }
    else {
        out <- list(time = tstar, Sc = pstar, x=x, tstar=tstar, pstar=pstar)
    }
    class(out) <- "constrainedKM"
    out
    # }}}    
}
