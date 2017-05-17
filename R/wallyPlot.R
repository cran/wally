# {{{ Wally plot function
# {{{ header

##' ##' Wally plots to assess calibration of a risk or survival prediction
##'
##' @title Wally plots to assess calibration of a risk or survival prediction
##' @param object Probabilistic survival predictions or probabilistic event risk predictions
##' evaluated at \code{time} for the subjects in \code{data}. Either
##' given in form of a numeric vector of probabilistic predictions or
##' as an object which  has a \code{predictRisk} method
##' @param time Time interest for evaluating calibration of the
##' predictions.
##' @param formula A survival or event history formula. The left hand
##' side is used to compute the expected event status. If
##' \code{formula} is \code{missing}, try to extract a formula from
##' the first element in object.
##' @param data A data frame in which to validate the prediction
##' models and to fit the censoring model. If \code{data} is missing,
##' try to extract a data set from the first element in object.
##' @param cause For competing risks settings the cause of interest.
##' @param q The number of quantiles. Defaults to 10.
##' @param ylim Limits of y-axis. If missing the function tries to
##' find appropriate limits based on the simulated and real data.
##' @param hanging If \code{TRUE}, hang bars corresponding to observed
##' frequencies at the value of the corresponding prediction.
##' @param seed A seed value to make results reproducible.
##' @param mar Plot margins passed to par.
##' @param colbox Color of the box which identifies the real data
##' calibration plot.
##' @param verbose If \code{TRUE} warn about missing formula and data.
##' @param col Colour of the bars. 
##' @param xlab Label for x-axis
##' @param labels Label below the bars. Either \code{"quantiles"} or \code{"quantiles.label"}
##' @param ... Further arguments passed to the subroutine \code{wallyCalPlot} and if \code{hanging}
##' is \code{TRUE} also to subroutine \code{lines}.
##' @return List of simulated and real data.
##' @references
##' Blanche P, Gerds T A, Ekstrom C T (2017). The Wally plot approach to assess the
##' calibration of clinical prediction models, submitted.
##'
##'
##' 
##' @examples
##'
##' # Survival setting
##' library(prodlim)
##' library(survival)
##' set.seed(180)
##' d = SimSurv(180)
##' f = coxph(Surv(time,status)~X1+X2,data=d,x=TRUE)
##' \dontrun{
##' wallyPlot(f,
##'           time=4,
##'           q=10,
##'           data=d,
##'           formula=Surv(time,status)~1)
##'  wallyPlot(f,
##'           time=4,
##'           q=10,
##'           hanging=TRUE,
##'           data=d,
##'           formula=Surv(time,status)~1)
##' }
##' 
##' # Competing risks setting
##' library(prodlim)
##' library(survival)
##' library(riskRegression)
##' set.seed(180)
##' d2 = SimCompRisk(180)
##' f2 = CSC(Hist(time,event)~X1+X2,data=d2)
##' \dontrun{
##' wallyPlot(f2,
##'           time=5,
##'           q=3,
##'           hanging=TRUE,
##'           data=d2,
##'           formula=Hist(time,event)~1)
##'           
##' }
##' 
##' # Reproduce Wally plots presented in Blanche et al. (2017)
##' \dontrun{
##' data(threecity)
##' wallyPlot(threecity$pi,
##' time=5,
##' hanging=TRUE,
##' formula=Hist(time,status)~1,
##' data=threecity,
##' ylim=c(-.1,.25),
##' seed= 511,
##' hline.lwd=3,
##' mar=c(1.01, 4.1, 1.15, 2))
##' }
##'
##' \dontrun{
##' data(divat)
##' wallyPlot(divat$pi,
##' time=5,
##' hanging=TRUE,
##' formula=Hist(time,status)~1,
##' data=divat,
##' ylim=c(-.1,.60),
##' seed= 123459,
##' hline.lwd=3,
##' mar=c(1.01, 4.1, 1.15, 2))
##' }
##' 
##' 
##' @export 
##' @author Paul F. Blanche <paulfblanche@@gmail.com> and Thomas A. Gerds <tag@@biostat.ku.dk>
# }}}
wallyPlot <- function(object,
                      time,
                      formula,
                      data,
                      cause=1,
                      q=10,
                      ylim,
                      hanging=FALSE,
                      seed=NULL,
                      mar=c(4.1,4.1, 2, 2),
                      colbox="red",
                      verbose=TRUE,
                      col=c("grey90","grey30"),
                      xlab="Risk groups",
                      labels="quantiles.labels",
                      ...){
    # {{{ prepare graphics output
    oldpar <- par(no.readonly = TRUE)
    par(mar = mar)
    par(mfrow=c(3,3),oma=c(3,0,0,0))
    # }}}
    # {{{ data & formula
    if (missing(data)){
        trydata <- try(data <- eval(object$call$data),silent=TRUE)
        if (("try-error" %in% class(trydata))|| match("data.frame",class(data),nomatch=0)==0)
            stop("Argument data is missing.")
        else
            if (verbose)
                warning("Argument data is missing. I use the data from the call to the first model instead.")
    }
    if (data.table::is.data.table(data)) data <- as.data.frame(data)
    if (missing(formula)){
        tryformula <- try(as.character(object$call$formula),silent=TRUE)
        if (("try-error" %in% class(tryformula))||length(grep("~",as.character(object$call$formula)))==0){
            stop(paste("Argument formula is missing and first model has no usable formula."))
        } else{
            ftry <- try(formula <- eval(object$call$formula),silent=TRUE)
            if ((class(ftry)=="try-error") || match("formula",class(formula),nomatch=0)==0)
                stop("Argument formula is missing and first model has no usable formula.")
            else if (verbose)
                warning("Formula missing. Using formula from first model")
            ## remove covariates
            formula <- update.formula(formula,".~1")
        }
    }
    m <- model.frame(formula,data,na.action=na.fail)
    response <- model.response(m)
    if (match("Surv",class(response),nomatch=FALSE))
        model.type <- "survival"
    else
        model.type <- attr(response,"model")
    if (is.null(model.type) & length(unique(response))==2)
        model.type <- "binary"
    if (!(model.type=="binary")){
        neworder <- order(response[,"time"],-response[,"status"])
        response <- response[neworder,,drop=FALSE]
        Y <- response[,"time"]
        status <- response[,"status"]
        if (model.type=="competing.risks"){
            event <- as.numeric(response[,"event"])
            event[status==0] <- 0
            status <- event
            rm(event)
        }
        data <- data[neworder,]
        if (missing(time)) time <- median(Y)
        else if (length(time)>1)
            stop("Please specify only one time point.")
    }
    # }}}
    # {{{ extract predicted risks
    if (class(object)[1] %in% c("numeric","double")) object <- matrix(object,ncol=1)
    predict.args <- list(object,newdata=data)
    if (model.type=="survival") predict.args <- c(predict.args,times=time)
    if (model.type=="competing.risks") predict.args <- c(predict.args,times=time,cause=cause)
    pred <- as.vector(do.call(riskRegression::predictRisk,predict.args))
    ## if given as a vector we need to re-order
    if (class(object)[[1]]%in% c("matrix","numeric")) pred <- pred[neworder]
    if (any(is.na(pred))) stop("Missing values in prediction. Maybe time interest needs to be set earlier?")
    # }}}
    # {{{ define risk groups according to quantiles 
    quant <- quantile(pred,seq(0,1,1/q))
    riskGroups <- cut(pred,breaks=quant,labels=1:(length(quant)-1),include.lowest=TRUE)
    predictedRisk <- tapply(pred,riskGroups,mean)
    # }}}
    # {{{ create labels to appear below the bars
    if ((is.logical(labels[1]) && labels[1]==TRUE) || labels[1] %in% c("quantiles.labels","quantiles")){
        qq <- quant
        if (labels[1]=="quantiles.labels"){
            pp <- seq(0,1,1/q)
            barlabels <- paste0("(",
                                sprintf("%1.0f",100*pp[-length(pp)]),",",
                                sprintf("%1.0f",100*pp[-1]),
                                ")\n",
                                sprintf("%1.1f",100*qq[-length(qq)])," - ",
                                sprintf("%1.1f",100*qq[-1]))
        } else 
            barlabels <- paste0(sprintf("%1.1f",100*qq[-length(qq)])," - ",
                                sprintf("%1.1f",100*qq[-1]))
    }else{
        barlabels <- ""
    }
    # }}}
    # {{{ Test if time is ok
    if (any((tooshortFollowup <- tapply(Y,riskGroups,max))<time)){
        stop(paste0("The following risk groups have too short followup. Nth quantile of predicted risks: ",
                    paste(names(tooshortFollowup),collapse=", "),
                    "\nThe minimum of the maximal followup in the risk groups is: ",
                    signif(min(tooshortFollowup),2)))
    }
    # }}}
    # {{{ set seed for making the Wally plot repoducible
    set.seed(seed)
    Allseed <- round(runif(length(levels(riskGroups))+1)*10000)
    # }}}
    # {{{ loop to create the simulated data sets
    ## - we loop over the groups (of subjects of homogeneous predicted risks)
    ## - we call the function to create 8 Bootstrap new datasets under the calibration assumption
    ## - we merge data from all subgroups to create 8 new datasets which mimic the actual data under the calibration assumption.
    ListofListsPerGroups <- vector("list", length(levels(riskGroups)))
    for(g in 1:length(levels(riskGroups))){
        Predg <- mean(pred[riskGroups==(levels(riskGroups)[g])])
        datag <- data[riskGroups==(levels(riskGroups)[g]),all.vars(formula)]
        names(datag) <- c("time","status")      
        if(model.type!="competing.risks"){
            # {{{ without competing risks
            ListofListsPerGroups[[g]] <- GenSurvDataKMCons(tstar=time,
                                                           pstar=1-Predg, # we want risk not survival
                                                           M = 8,
                                                           data=datag,
                                                           myseed=Allseed[g])
            # }}}
        }else{
            # {{{ with competing risks
            # we should seperate the two cases : wether or not there are competing events
            # Here is the case where we observe the two competing events
            if(c(1) %in% unique(datag$status) & c(2) %in% unique(datag$status)){
                # The current version of the GenCRDataAJCons function does not handle ties,
                # so we add a little bit of jitter if there are ties
                # {{{ to handle ties with competing risks
                if(!all(!duplicated(datag$time))){
                    thelocations <- which(duplicated(datag$time))
                    set.seed(Allseed[g])
                    datag$time[thelocations] <- jitter(datag$time[thelocations],
                                                       amount=NULL,
                                                       factor=1/(length(datag$time)*2))
                }
                # }}}
                ListofListsPerGroups[[g]] <- GenCRDataAJCons(tstar=time,
                                                             pstar=Predg,
                                                             M = 8,
                                                             data=datag,
                                                             myseed=Allseed[g])$l
            }else{
                # Here is the case where we observe only one of the two competing events                
                # if only event 2 is observerved, we need a trick (transform to event 1 to generate the data
                # and then transform back to event 1)
                if(c(2) %in% unique(datag$status)){
                    # case in which we do not observe event 1 (the main event)
                    datag$status[datag$status==2] <- 1                                                                        
                    ListofListsPerGroups[[g]] <- GenSurvDataKMCons(tstar=time,
                                                                   pstar=1-Predg, # be careful here we want prob of event ! (not survival..)
                                                                   M = 8,
                                                                   data=datag,
                                                                   myseed=Allseed[g])
                    for(j in 1:8){
                        datag018 <- ListofListsPerGroups[[g]][[j]]
                        datag018$status[datag018$status==1] <- 2 # transform back 1 to 2
                        ListofListsPerGroups[[g]][[j]] <- datag018
                    }
                }else{
                    # case in which we do not observe event 2 (the competing event),
                    # we can do as without any competing risk
                    ListofListsPerGroups[[g]] <- GenSurvDataKMCons(tstar=time,
                                                                   pstar=1-Predg, # be careful here we want prob of event ! (not survival..)
                                                                   M = 8,
                                                                   data=datag,
                                                                   myseed=Allseed[g])
                }
            }
        }
        # }}}
        # {{{ Merge data to create the new data
        # add predictions to the table and remove useless column
        for(j in 1:8){
            datag018 <- ListofListsPerGroups[[g]][[j]]
            ## ListofListsPerGroups[[g]][[j]] <- cbind.data.frame(datag018[,c("time","status")],riskGroups=riskGroups[riskGroups==(levels(riskGroups)[g])])
            ListofListsPerGroups[[g]][[j]] <- cbind.data.frame(datag018[,c("time","status")],
                                                               riskGroups=riskGroups[riskGroups==(levels(riskGroups)[g])],
                                                               risk=pred[riskGroups==(levels(riskGroups)[g])])
        }
        # }}}    
    }
    DataList <- lapply(1:8,function(i){do.call("rbind", sapply(lapply(ListofListsPerGroups,"[",i),"[",1))})
    ## add real data at the last position
    realData <- cbind(data[,all.vars(update(formula,".~1")),drop=FALSE],riskGroups=riskGroups)
    DataList <- c(DataList,list(realData))
    # }}}  
    # {{{ sample the order of the data in the list
    pos <- 1:9
    set.seed(Allseed[length(Allseed)])
    pos <- sample(pos)        
    figpos <- order(pos)[9] # where is the actual plot (with the true data)
    # }}}
    # {{{ control of plot arguments
    superuser.defaults <- list(hide=TRUE,choice=6,zoom=FALSE)
    hline.defaults <- list(col=2, lwd=2.5, type="s")
    ## if (missing(ylim)) if (hanging) ylim <- c(-1,1) else ylim <- c(0,1)
    if (missing(ylim)) {
        user.ylim <- FALSE
        if (hanging) ylim <- c(-1,1) else ylim <- c(0,1)
    }else{
        user.ylim <- TRUE
    }
    axis2.DefaultArgs <- list(side=2,las=2,at=seq(0,ylim[2],ylim[2]/4),mgp=c(4,1,0))
    legend.DefaultArgs <- list(legend=c("Predicted risks","Observed frequencies"),col=col,cex=par()$cex,bty="n",x="topleft")
    
    barlabels.DefaultArgs <- list(cex=.7*par()$cex,
                                  y=c(-abs(diff(ylim))/15,-abs(diff(ylim))/25),
                                  text=barlabels)
    frequencies.DefaultArgs <- list(cex=.7*par()$cex,percent=FALSE,offset=0)
    lines.DefaultArgs <- list(type="l")
    abline.DefaultArgs <- list(lwd=1,col="red")
    barplot.DefaultArgs <- list(ylim = ylim,
                                col=col,
                                axes=FALSE,
                                ylab="",
                                xlab=xlab,
                                beside=TRUE,
                                legend.text=NULL,
                                cex.axis=par()$cex.axis,
                                cex.lab=par()$cex.lab)
    smartA <- prodlim::SmartControl(call= list(...),
                                    keys=c("superuser","barplot","legend","axis2","abline","barlabels","frequencies","hline"),
                                    ignore=NULL,
                                    ignore.case=TRUE,
                                    defaults=list("superuser"=superuser.defaults,
                                        "barplot"=barplot.DefaultArgs,
                                        "abline"=abline.DefaultArgs,
                                        "legend"=legend.DefaultArgs,
                                        "barlabels"=barlabels.DefaultArgs,
                                        "frequencies"=frequencies.DefaultArgs,
                                        "axis2"=axis2.DefaultArgs,
                                        "hline"=hline.defaults),
                                    forced=list("abline"=list(h=0)),
                                    verbose=TRUE)
    # }}}
    # {{{ loop to create the 9 plots
    TabList <- vector("list",9)
    printleg <- c(rep(FALSE,4),TRUE,rep(FALSE,4))
    for(i in pos){
        if (i==9) riskformula <- formula else riskformula <- formula("Hist(time,status)~1")
        # compute observed risk in groups
        xgroups <- (quant[-(length(quant))]+quant[-1])/2
        groupFormula <- update(riskformula,paste(".~riskGroups"))
        Obs <- as.vector(riskRegression::predictRisk(prodlim::prodlim(groupFormula,data=DataList[[i]]),
                                                     cause=cause,
                                                     newdata=data.frame(riskGroups=levels(riskGroups)),
                                                     times=time))
        if (any(is.na(Obs))) {
            Obs[is.na(Obs)] <- max(Obs,na.rm=TRUE)
            warning("Missing values in expected frequencies. Maybe too many quantiles relative to the number of observations? Or time interest set too late?")
        }
        TabList[[i]] <- list(Pred=predictedRisk,
                             Obs=Obs,
                             time=time,
                             cause=cause,
                             summary=summary,
                             control=smartA,
                             legend=legend,
                             diag=diag,
                             legend=FALSE,
                             labels=labels,
                             model.type=model.type,
                             hanging=hanging,
                             ## showFrequencies=1L,
                             showFrequencies=FALSE,
                             col=col,
                             ylim=ylim,
                             axes=1L)
    }
    if (hanging){
        ## FIXME
        minY <- min(sapply(TabList,function(x){min(x$Pred-x$Obs)}))
        maxY <- max(sapply(TabList,function(x){max(x$Pred)}))
        # rangeY <- range(sapply(TabList,function(x){c(min(x$Pred-x$Obs),max(x$Pred))}))
        minY <- min(0,seq(-1,1,0.05)[prodlim::sindex(eval.times=minY,jump.times=seq(-1,1,0.05))])
    } else{
        maxY <- max(sapply(TabList,function(x){max(x$Obs)}))
        maxY <- c(seq(0,1,0.05),1)[1+prodlim::sindex(eval.times=maxY,jump.times=seq(0,1,0.05))]
        minY <- 0
    }
    if (user.ylim){
        if (minY<0) minY <- min(ylim,minY)
        maxY <- ylim[2]
    }
    for (j in 1:9){
        i = pos[j]
        px <- TabList[[i]]
        px$control$barplot$ylim <- c(minY,maxY)
        px$control$barlabels$y <- c(-abs(diff(c(minY,maxY)))/15,-abs(diff(c(minY,maxY)))/25)
        ## need to round to avoid strange results 
        px$control$axis2$at <- sort(round(seq(minY,maxY,(maxY-minY)/4),3))
        ## place a legend below panel 8
        if (j==8) {
            px$legend <- TRUE
            px$control$legend$legend <- c("Predicted risk","Observed frequency")
            px$control$legend$xpd <- NA
            px$control$legend$x <- "bottom"
            px$control$legend$ncol=2
            px$control$legend$inset=c(0,-0.4)
            px$control$legend$cex <- 1.3
        }else{
            px$legend <- FALSE
        }
        ## add the plot to the grid
        px$labels <- FALSE
        calibrationBarplot(px)
        upleft <- par("usr")[c(1,4)]
        points(x=upleft[1],y=upleft[2],xpd=NA,pch=19,cex=5,col="orange")
        text(j,x=upleft[1],y=upleft[2],xpd=NA,col="black",cex=1.2)
    }
    ## mtext(side=1,line=1,"Can you find wally?",cex=1.5*par()$cex)
    # }}}
    # {{{ interact with user
    if (is.null(smartA$superuser$choice)){
        par(oldpar)
        invisible(TabList[order(pos)])
    }else{
        if (smartA$superuser$hide!=FALSE){
            # {{{ ask to show the actual plot
            xx <-  select.list(1:9,
                               multiple=FALSE,
                               title="\nWhere is Wally? Can you find the plot which is based on the real data?\nSelect an orange number: ")
            # }}}
        }else xx <- smartA$superuser$choice
        par(mfg = c(xx%/%3.1 + 1, xx - (xx%/%3.1) * 3))
        box(col = "green", lwd= 3)
        # {{{ show the actual plot by adding a red box

        if (xx==figpos) {
            mtext("Correct: real data!",col="green",cex=1.2,line=-3,xpd=NA)
        }else{
            mtext("Not correct!",col="red",cex=1.2,line=-3,xpd=NA)
            par(mfg = c(figpos%/%3.1 + 1, figpos - (figpos%/%3.1) * 3))
            mtext("Real data",col="red",cex=1.2,line=-3,side=3,xpd=NA)
            box(col = colbox, lwd = 5)
        }

        # }}}
        # {{{ ask to press a key to zoom in on the actual plot
        ## readline("Hit <Enter> to better show the original plot. ")
        if (smartA$superuser$zoom==FALSE && smartA$superuser$hide!=FALSE){
            zoom <- select.list(c("yes","no"),title="Zoom in on real data calibration plot? ")
        }else zoom <- "no"
        # }}}
        # {{{ show the actual plot

        if((smartA$superuser$zoom==TRUE) || (zoom=="yes")){
            par(mfrow=c(1,1),oma=c(2,2,2,2),mar=c(4.1,4.1, 4.1, 4.1))
            px <- TabList[[9]]
            px$control$barplot$ylim <- c(minY,maxY)
            px$control$barlabels$y <- c(-abs(diff(c(minY,maxY)))/15,-abs(diff(c(minY,maxY)))/25)
            ## need to round to avoid strange results like in
            px$control$axis2$at <- sort(round(seq(minY,maxY,(maxY-minY)/4),3))
            px$control$barplot$legend.text=c("Predicted risk","Observed frequency")
            px$control$legend$legend <- c("Predicted risk","Observed frequency")
            px$control$barplot$xlab <- "Risk groups"
            px$control$legend$xpd <- NA
            px$control$legend$x <- "top"
            px$control$legend$ncol <- 2
            px$control$legend$inset <- c(0,-0.2)
            px$showFrequencies <- TRUE
            px$legend=TRUE
            px$labels <- labels
            calibrationBarplot(px)
        }

        # }}}
        par(oldpar)
    }
    # }}}
    invisible(list(Tables=TabList[order(pos)],Data=DataList))
}
# }}}
