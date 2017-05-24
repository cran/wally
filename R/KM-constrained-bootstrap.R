# Main Author for this function:  Paul Blanche

# Description :
#
# This function generates 'M' new data sets from NPMLE under the constraint
# that the survival function is 'pstar' at 'tstar'.



GenSurvDataKMCons <- function (tstar, pstar, M = 1000, data, myseed=140786) {
    old <- .Random.seed
    # {{{ set seed
    set.seed(myseed)
    # }}}
    # number of observations to generate
    n <- nrow(data)
    # {{{ To generate time to event
    # fit KM with and without constraint
    x <- KMconstrained(tstar, pstar,  data)
    survTimes <- c(x$time, Inf)
    dSurv <- -diff(c(1, x$Sc))
    dSurv <- abs(c(dSurv, 1 - sum(dSurv))) # absolute to avoid problems without long doubles
    # }}}
    # {{{ to generate censoring
    # fit KM for censoring
    if(sum(data$status==0)>0){
        fitcens <- prodlim(Hist(time,status)~1,data=data,reverse=TRUE)
        # censored observation times
        censTimes <- sort(unique(data$time[data$status==0]))
        # compute probability with which the observed times are drawn
        survCens <- predict(fitcens,times=censTimes,type="surv")
        censTimes <- c(censTimes, max(data$time) + 1)
        dCens <- -diff(c(1, survCens))
        dCens <- abs(c(dCens, 1 - sum(dCens))) # absolute to avoid problems without long doubles
    }
    # }}}
    # {{{ function to generate new data    
    ftoloop <- function(i){
        ## browser()
        if (sum(data$status==0)>0) {
            Ci <- sample(censTimes, n, prob = dCens, replace = TRUE)
        }
        else{
            Ci <- rep(max(data$time) + 1, n)
        }
        Ti <- sample(survTimes, n, prob = dSurv, replace = TRUE)
        Time <- pmin(Ti,Ci)
        Status <- as.numeric(Ti<=Ci)
        d <- data.frame(time=Time,
                        status=Status,
                        eventtime=Ti,
                        censtime=Ci
                        )
        d <- d[order(d$time),]
        d
    }
    ## browser()
    ## ftoloop(1)
    # }}}    
    # {{{ output
    out <- lapply(1:M,ftoloop)
    # }}}
    .Random.seed <<- old
    out
}
