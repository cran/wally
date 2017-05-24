# Main Author for this function:  Paul Blanche

# Description :
#
# This function generates 'M' new data sets from NPMLE under the constraint
# that the cumulative incidence function of event 1 is 'pstar' at 'tstar' 

GenCRDataAJCons <- function (tstar, pstar, M = 1000, data, myseed=140786) 
{
    old <- .Random.seed
    # {{{ set seed
    set.seed(myseed)
    # }}}
    # number of observations to generate
    n <- nrow(data)
    # {{{ To generate time to event and cause of event
    # CR fit with and without constraint
    x <- AJconstrained(tstar, pstar,  data)
    survTimes <- c(x$time, Inf)
    Surv <- 1 - x$CIF1c - x$CIF2c
    Survminus <- c(1,Surv[-length(Surv)])
    dSurv <- -diff(c(1, Surv))
    dSurv <- abs(c(dSurv, 1 - sum(dSurv))) # absolute to avoid problems without long doubles
    diffcumInc1 <- diff(c(0,x$CIF1c))
    diffcumInc2 <- diff(c(0,x$CIF2c))   
    lambda1 <- diffcumInc1/Survminus
    lambda2 <- diffcumInc2/Survminus
    lambda <- lambda1 + lambda2
    ## cbind(survTimes,c(Surv,0),dSurv,c(0,x$CIF1c),c(0,x$CIF2c)) #diffcumInc1,diffcumInc2)
    # }}}
    ## browser()
    # {{{ to generate censoring
    # fit KM for censoring
    if(sum(data$status==0)>=1){ # only if there is censoring
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
        if(sum(data$status==0)>=1){ # if there are censored observations
            Ci <- sample(censTimes, n, prob = dCens, replace = TRUE)
        }
        else{
            Ci <- rep(max(data$time) + 1, n)
        }
        # {{{ generate time and cause of event
        # first generate the random number of times that each observed time is selected
        HowMany <- rmultinom(1, size=n, prob=dSurv)        
        # then create the vector of observed times by repeated each observed time the suitable number of times
        Ti <- rep(survTimes,times=HowMany)
        # generate cause
        probtogenevent <- rep(c((lambda2/(lambda1 + lambda2)),0.5),times=HowMany) # 0.5 because it does not matter after time point of interest       
        events <- rbinom(length(Ti),size=1,prob=probtogenevent ) + 1       
        # }}}
        Time <- pmin(Ti,Ci)
        Status <- as.numeric(Ti<=Ci)*events
        d <- data.frame(time=Time,
                        status=Status,
                        event=events,
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
    out <- list(l=lapply(1:M,ftoloop),x=x)
    # }}}
    .Random.seed <<- old
    out
}
