calibrationBarplot <- function(x,...){
    control <- x$control
    hanging <- x$hanging
    Pred <- x$Pred
    Obs <- x$Obs
    # {{{ legend
    if(is.logical(x$legend[1]) && x$legend[1]==FALSE){
        control$barplot$legend.text <- NULL
    }else{
        control$barplot$legend.text <- control$legend$legend
        control$barplot$args.legend <- control$legend
    }
    # }}}
    if (is.null(control$barplot$space)) control$barplot$space <- rep(c(1,0),length(Pred))
    # {{{ height of the bars
    PredObs <- c(rbind(Pred,Obs))
    control$barplot$height <- PredObs
    if (hanging){
        control$barplot$offset <- c(rbind(0,Pred-Obs))
        minval <- min(Pred-Obs)
        if (minval<0)
            negY.offset <- 0.05+seq(0,1,0.05)[prodlim::sindex(jump.times=seq(0,1,0.05),eval.times=abs(minval))]
        else
            negY.offset <- 0
        control$barplot$ylim[1] <- min(control$barplot$ylim[1],-negY.offset)
        control$barlabels$y <- control$barlabels$y-negY.offset
    }
    # }}}
    # {{{ call barplot
    coord <- do.call("barplot",control$barplot)
    # }}}
    # {{{ label bars (below)
    if (is.character(x$labels)){
        mids <- rowMeans(matrix(coord,ncol=2,byrow=TRUE))
        text(x=mids,
             y=rep(control$barlabels$y,length.out=length(mids)),
             labels=control$barlabels$text,
             xpd=NA,
             cex=control$barlabels$cex)
    }
    # }}}
    # {{{ horizontal line
    if (hanging){
        do.call("abline",control$abline)
    }
    # }}}
    # {{{ labels on top of the bars
    if (x$showFrequencies){
        if(hanging){
            text(x=coord,
                 cex=control$frequencies$cex,
                 pos=3,
                 y=(as.vector(rbind(Pred,Pred)) +rep(control$frequencies$offset,times=length(as.vector(coord))/2)),
                 paste(round(100*c(rbind(Pred,Obs)),0),ifelse(control$frequencies$percent,"%",""),sep=""),xpd=NA)
        }else{
            text(x=coord,
                 pos=3,
                 y=c(rbind(Pred,Obs))+control$frequencies$offset,
                 cex=control$frequencies$cex,
                 paste(round(100*c(rbind(Pred,Obs)),0),
                       ifelse(control$frequencies$percent,"%",""),
                       sep=""),xpd=NA)
        }
    }
    # }}}
    # {{{ Paulo add for ISCB
    if (hanging){
        offset <- control$barplot$offset
        hline.args <- c(list(x=c(0,1,ceiling(coord[,1]),ceiling(coord[,1])[length(ceiling(coord[,1]))],ceiling(coord[,1])[length(ceiling(coord[,1]))]+1),
                             y= c(offset[1],offset,offset[length(offset)],rep(offset[1],2))),
                        control$hline)
        do.call("lines",hline.args)
    }
    # }}}
    # {{{ axis
    if (x$axes){
        control$axis2$labels <- paste(100*control$axis2$at,"%")
        do.call("axis",control$axis2)
    }
    invisible(coord)
    # }}}
}
#----------------------------------------------------------------------
### plotCalibrationPlot.R ends here
