# Main Author for this function:  Paul Blanche

# Description :
#
# This function generates NPMLE under the constraint
# that the cumulative incidence function of event 1 is 'pstar'
# at 'tstar'.
# This current version assumes there are no ties !

AJconstrained <- function (tstar, CIF1star,  data, mytol=10^(-5), plotUNBUG = FALSE) 
{
    # {{{ data management steps (just fo rme)
    data <- data[order(data$time),]
    time <- data$time
    status <- data$status
    # {{{ help for unbugging
    ## plot(prodlim(Hist(time,status)~1,data=data))
    ## plot(prodlim(Hist(time,status)~1,data=data),cause=2,col="red",add=TRUE)
    ## abline(v=tstar,lty=2)
    # }}}
    # }}}
    # {{{ STOP : if ties
    if(length(unique(time))!=length(time)){
        stop("This function does not handle ties to compute the constrained version of the Aalen-Johansen estimates.")
    }
    # }}}
    # {{{ handle the case where no main event is observed
    if(sum(data$time<=tstar & data$status==1)==0){
        # we should had an observed event in that case to allow a jump of the CIF at that time !
        # first we check if an observed time corresponds to tstar, if yes, we had an event.
        # If no, we add a time to the data
        if(!(tstar %in% data$time)){
            time <- c(time,tstar)
            status <- c(status,1)
            status <- status[order(time)]
            time <- time[order(time)]
            ## cbind(time,status)
        }else{
            status[which(tstar==time)] <- 1
        }    
        # Be careful, when compute for lambda=0, then we want the original data !
    }
    ## browser()
    # }}}
    # {{{ useful quantities
    N <- length(time) 
    ## tab <- table(time, status)
    tab <- table(time,factor(status,levels=c(0,1,2)))
    d1i <- tab[, 2]      
    d2i <- tab[, 3]      
    ci <- tab[, 1]
    ti <- as.numeric(dimnames(tab)[[1]])
    k <- length(ti)
    ni <- c(N, N - cumsum(ci)[-k] - cumsum(d1i)[-k] - cumsum(d2i)[-k])
    # keep only for ti before tstar
    II <- ti <= tstar    
    nj <- ni[II]
    d1j <- d1i[II]
    d2j <- d2i[II]
    cj <- ci[II]
    kj <- sum(II)
    ## cbind(ti[II],nj,d1j,d2j,cj)
    ## browser()
    # }}}
    # {{{ function to compute estimates given lamda
    CIFlambda <- function(lambda) {
        if(kj==0){
            out <- list(forRoot=0,
                        Surv=1,
                        a=0,
                        a1=0,
                        a2=0,
                        CIF1=0,
                        CIF2=0,
                        alla1a2Positive=TRUE,
                        alla1Positive=TRUE,
                        alla2Positive=TRUE,
                        allOK=TRUE)            
            return(out)
        }
        # initialize vectors
        a1 <- rep(0,kj) # cause-1-specific hazard
        a2 <- rep(0,kj) # cause-2-specific hazard
        a <-  rep(0,kj) # all-causes hazard
        Surv <- rep(1,kj) # survival : P(T > t_i)
        CIF1 <- rep(0,kj) # CIF 1:  P(T<=t_i, eta=1)
        CIF2 <- rep(0,kj) # CIF 1:  P(T<=t_i, eta=2)
        # loop over ti
        for (i in 1:kj) {
            if(i==1){
                if(d1j[i]>0){
                    # event at ti is event 1
                    a1[i] <- 1/(nj[i] + lambda*(1-0-CIF1star) )
                }else{
                    if(d2j[i]>0){
                        # event at ti is event 2
                        a2[i] <- 1/(nj[i] + lambda*(0-CIF1star))
                    }
                } # else, i.e. if event at ti is censoring, then we keep a2[i]=a1[i]=0
                a[i] <- a1[i] + a2[i]    
                Surv[i] <- 1*(1- a[i])
                CIF1[i] <- 0 + 1*a1[i]
                CIF2[i] <- 0 + 1*a2[i]                
            }else{
                if(d1j[i]>0){
                    # event at ti is event 1
                    a1[i] <-  1/(nj[i] + lambda*(1 - CIF2[i-1] - CIF1star) )
                }else{
                    if(d2j[i]>0){
                        # event at ti is event 2
                        a2[i] <- 1/(nj[i] + lambda*(CIF1[i-1] - CIF1star))
                    }
                }
                # else, i.e. if event at ti is censoring, then we keep a2[i]=a1[i]=0
                a[i] <- a1[i] + a2[i]    
                Surv[i] <- Surv[i-1]*(1-a[i])
                CIF1[i] <- CIF1[i-1] + Surv[i-1]*a1[i]
                CIF2[i] <- CIF2[i-1] + Surv[i-1]*a2[i]                               
            }
        }
        # output
        out <- list(forRoot=CIF1[kj],
                    Surv=Surv,
                    a=a,
                    a1=a1,
                    a2=a2,
                    CIF1=CIF1,
                    CIF2=CIF2,
                    alla1a2Positive=all(c(a1>=0,a2>=0)),
                    alla1Positive=all(a1>=0),
                    alla2Positive=all(a2>=0),
                    allOK=all(c(a1>=0,
                        a2>=0,
                        CIF1<=1,
                        CIF1>=0,
                        CIF2<=1,
                        CIF2>=0,
                        CIF1 + CIF2<=1
                        ))
                    )
        out
    }
    # }}}
    # {{{ find lambda
    ## browser()
    # {{{ we find interval for lambda with positive values for all aki in which to use uniroot
    ## print(paste("lambda.min : start"))
    thelambdamin <- -N*1000
    try(thelambdamin <- uniroot(f=function(l){as.numeric(CIFlambda(l)$allOK) - 0.5},
                                lower=-N*1000, # How to choose that ????
                                upper=0,
                                tol=mytol)$root,silent=TRUE)
    ## browser()
    if(!CIFlambda(thelambdamin)$allOK){
        nwhilellop <- 0
        while(!CIFlambda(thelambdamin)$allOK & nwhilellop<30){
            ## print(paste("While loop (1), step ",nwhilellop,": thelambdamin=",thelambdamin,"start"))
            nwhilellop <- nwhilellop + 1
            thelambdamin <- uniroot(f=function(l){as.numeric(CIFlambda(l)$allOK) - 0.5},
                                    lower=thelambdamin,
                                    upper=0,
                                    ## extendInt = c("downX"),
                                    tol=mytol)$root
            ## print(paste("While loop (1), step ",nwhilellop,": thelambdamin=",thelambdamin,"stop"))            
        }
    }
    ## print(paste("lambda.min : done"))
    ## print(paste("thelambdamin=",thelambdamin))
    ## browser()
    if(sum(data$status==2)>=1){ # only if there are cometing events
        ## browser()
        thelambdamax <- N*1000
        try(thelambdamax <- uniroot(f=function(l){as.numeric(CIFlambda(l)$allOK) - 0.5},
                                    lower=0,
                                    upper=N*1000, # How to choose that ????
                                    ## extendInt = c("upX"),
                                    tol=mytol)$root,silent=TRUE)
        ## print(paste("lambda.max : first done"))    
        if(!CIFlambda(thelambdamax)$allOK){
            nwhilellop2 <- 0
            while(!CIFlambda(thelambdamax)$allOK & nwhilellop2<30){
                ## print(paste("While loop (2), step ",nwhilellop2,": thelambdamax=",thelambdamax))
                nwhilellop2 <- nwhilellop2 + 1
                thelambdamax <- uniroot(f=function(l){as.numeric(CIFlambda(l)$allOK) - 0.5},
                                        lower=0,
                                        upper=thelambdamax,
                                        tol=mytol)$root
                ## print(paste("while loop :",nwhilellop2))
            }
        }
    }else{
        thelambdamax <- N*1000
    }
    ## print(paste("thelambdamax=",thelambdamax))
    ## print(paste("lambda.max : done"))
    # }}}
    # {{{ to understand bugs    
    ## browser()
    if(plotUNBUG){
        lengthxxx <- 1000
        ## xxx <- seq(from=thelambdamin-1,to=thelambdamax+1,length.out=lengthxxx)
        xxx <- seq(from=thelambdamin,to=thelambdamax,length.out=lengthxxx)
        yyy <- sapply(xxx,function(x){o <- CIFlambda(x);o$forRoot - CIF1star})
        plot(xxx,yyy,ylim=c(-1,1),col=ifelse(yyy<=0,"blue","red"))
        abline(h=0)
        abline(v=thelambdamin)
        abline(v=thelambdamax)
        ## abline(v=(-N/(1-CIF1star)),col="red")
        ## abline(v=(N/CIF1star),col="red")
        #    
        ## yyy[1]
        ## yyy[lengthxxx]
        ## CIFlambda(xxx[1])$alla1a2Positive
        ## CIFlambda(xxx[lengthxxx])$alla1a2Positive
        ## CIFlambda(xxx[1] + 5)$alla1a2Positive
        ## CIFlambda(xxx[1] + 100/N)$forRoot - CIF1star
        ## CIFlambda(xxx[lengthxxx])$alla1a2Positive
        #
        ## CIFlambda(thelambdamin)$alla1a2Positive
        ## CIFlambda(thelambdamax)$alla1a2Positive
        ## CIFlambda(thelambdamax)$forRoot
        ## CIFlambda(thelambdamax)$forRoot - CIF1star
        ## CIFlambda(thelambdamin)$forRoot - CIF1star
        ## points(thelambdamin,CIFlambda(thelambdamin)$forRoot - CIF1star,col="red")
        ## points(thelambdamax,CIFlambda(thelambdamax)$forRoot - CIF1star,col="red")
        ## CIFlambda(thelambdamax)$allOK
        ## CIFlambda(thelambdamax)$CIF1
        ## CIFlambda(thelambdamax)$CIF2
        ## CIFlambda(thelambdamax)$a
    }
    # }}}
    ## print(paste("Int lambda.min , lambda.max : done"))
    # {{{ find lambda within interval
    ## browser()
    thelambda <- uniroot(f=function(x){CIFlambda(x)$forRoot - CIF1star},
                         lower=(thelambdamin + 1/N),
                         upper=(thelambdamax -1/N))$root

        # {{{ to understand bugs
        if(plotUNBUG){
            abline(v=thelambda,col="red")
        }
        # }}}

    # }}}
    # }}}
    # {{{ compute everything with the right value for lambda
    outCIFlambda <- CIFlambda(thelambda)
    timeout <- time[II] # save the time points for we observe some jumps in CIF1 or CIF2 (as modify it just after)
    if(sum(data$time<=tstar & data$status==1)==0){
        # in that case we need to redefine what is used by CIFlambda
        data <- data[order(data$time),]
        time <- data$time
        status <- data$status      
        # {{{ useful quantities
        N <- length(time) 
        ## tab <- table(time, status)
        tab <- table(time,factor(status,levels=c(0,1,2)))
        d1i <- tab[, 2]      
        d2i <- tab[, 3]      
        ci <- tab[, 1]
        ti <- as.numeric(dimnames(tab)[[1]])
        k <- length(ti)
        ni <- c(N, N - cumsum(ci)[-k] - cumsum(d1i)[-k] - cumsum(d2i)[-k])
        # keep only for ti before tstar
        II <- ti <= tstar    
        nj <- ni[II]
        d1j <- d1i[II]
        d2j <- d2i[II]
        cj <- ci[II]
        kj <- sum(II)
        ## cbind(ti[II],nj,d1j,d2j,cj)
        # }}}       
    }
    outCIFlambda0 <- CIFlambda(0)
    # }}}
    # {{{ output    
    out <- list(n=N,
                time = timeout,
                CIF1c = outCIFlambda$CIF1,
                CIF2c = outCIFlambda$CIF2,
                tstar=tstar,
                CIF1star=CIF1star,
                outCIFlambda=outCIFlambda,
                outCIFlambda0=outCIFlambda0,
                lambda=thelambda,
                lambdaSeachInt=c(thelambdamin,thelambdamax)
                )
    class(out) <- "constrainedAJ"
    # }}}
    ## print(paste("CIF1star = ",CIF1star,": DONE"))
    out
}
