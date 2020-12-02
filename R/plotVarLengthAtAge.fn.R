#' Plot variability of length at age
#'
#' Plots the SD and CV of age at observed and predicted length
#'
#' @param dir A directory where you want to save the file if
#' \code{dopng = TRUE}. Inside this directory, an additional directory called
#' \code{"plots"} will be created if it doesn't already exist and the resulting
#' file will be saved there.
#' @param dat A data frame extracted from the NWFSC database, see
#' \code{\link{PullBio.fn}()}.
#' @param main File-path name that will be added to the beginning of
#' \code{"_VarLengthAtAge.png"} if \code{dopng = TRUE}.
#' If \code{dopng = TRUE} and \code{main = NULL}, which is the default,
#' the underscore is removed and nothing is added.
#' @param ageBin Currently fixed at 1, and does nothing in the code.
#' @param bySex TRUE/FALSE specifying if the plots should be sex-specific or not.
#' @param parStart Vector of starting parameters for Linf, K, L0, CV0, and CV1
#' in the estimation of von Bertalanffy growth parameters. Parameters are used
#' as initial vales if \code{estVB = TRUE} or
#' as estimates if \code{estBV = FALSE}. Provide the values in the exact order
#' as listed and in normal space because they will be logged transformed by
#' \code{PlotVarLengthAtAge} prior to using in \code{\link{GetVB.fn}()}.
#' @param estVB A logical value indicating whether or not to estimate the
#' von Bertalanffy growth parameters. The default is to estimate them, but if
#' \code{FALSE}, then the values in \code{parStart} will be used as input values
#' for calculating predicted lengths.
#' @param bins     The bins to put ages into. If NULL then simply uses the ages as recorded.
#' @param legX, legY The x- and y-coordinates to be used for the legend position.
#' The former can be specified using a keyword or coordinates, see the Details
#' section of \code{\link[grDevices]{xy.coords}} for more details.
#' The latter, \code{legY}, is \code{NULL} by default, which expects \code{legX}
#' to be a keyword such as \code{"bottomleft"}. This default value will
#' place the legend in the bottom left-hand corner of the figure.
#' @param dopng TRUE/FALSE specifying whether or not to save a png file, where
#' the default is the latter leading to nothing being saved.
#' @param ... Additional arguments passed to the four calls to
#' \code{\link{plot}()}.
#'
#' @author Allan Hicks and Chantel Wetzel
#' @export
#' @examples
#'\dontrun{
#' # SurveyName is only arg that has to be specified
#' bio_dat = PullBio.fn(SurveyName = "NWFSC.Combo")
#' PlotVarLengthAtAge.fn(dat = bio_dat, dopng = TRUE, dir = getwd())
#' dev.off()
#'}

PlotVarLengthAtAge.fn <- function(dir = NULL, dat, main = NULL, ageBin = 1,
                                  bySex = TRUE, parStart = c(52, 0.09, 1, 0.1, 0.1), estVB = TRUE,
                                  bins = NULL, legX = "bottomleft", legY = NULL,
                                  dopng = FALSE, ...)
{
    #calculate and plot the sd and cv for length at age
    #if you enter estVB=F, then it uses the parStart as the VB parameters

    dat <- dat[!is.na(dat$Length_cm), ]

    dat <- dat[!is.na(dat$Age),]
    dat <- dat[dat$Sex%in%c("F", "M"), ]
    

    datL <- dat[!is.na(dat$Age),]
    if(is.null(bins)) {datL$Age_2 <- datL$Age}
    if(!is.null(bins)) {datL$Age_2 <- findInterval(datL$Age,bins)}

    if(!bySex) {
        datL <- list(allSex=datL)
        nn <- 1
    }

    if(bySex) {
        datLf <- datL[datL$Sex=="F",]
        datLm <- datL[datL$Sex=="M",]
        datL <- list(female=datLf,male=datLm)            
        nn <- 2
    }
    
    if (dopng) {
        if (is.null(dir)){stop("Directory needs to be set.")} 
        if (!file.exists(dir)) { stop("The dir argument leads to a location", ",\ni.e., ", dir, ", that doesn't exist.") }
        plotdir <- file.path(dir, "plots")
        plotdir.isdir <- file.info(plotdir)$isdir
        if (is.na(plotdir.isdir) | !plotdir.isdir) {
            dir.create(plotdir)
        }
        fn.png <- file.path(dir, "plots",
          paste0(main, ifelse(is.null(main), "", "_"), "VarLengthAtAge.png"))
        if ( is.null(main)) { png(file.path(dir, "plots", "VarLengthAtAge.png"), height=7, width=7, units="in",res=300) }
        if (!is.null(main)) { png(file.path(dir, "plots", paste(main, "_VarLengthAtAge.png", sep = "")), height=7, width=7, units="in",res=300) }
        on.exit(dev.off(), add = TRUE)
    }

    par(mfcol=c(2,nn), mar =c(3,5,3,5))

    out <- vector(mode="list", length=nn)
    names(out) <- names(datL)
    outpars <- list()
    for(i in 1:length(datL)) {
        if(estVB) {
            xpar <- exp(stats::optim(par = log(c(parStart[2], parStart[1], parStart[3], parStart[4], parStart[5])),
                fn = GetVB.fn,
                Ages = datL[[i]]$Age, Lengths = datL[[i]]$Length_cm)$par)
            cat("Estimated VB parameters for",names(datL)[i],xpar,"\n")
        }
        if(!estVB) {
            xpar <- parStart
        }
        outpars[[length(outpars) + 1]] <- xpar
        predL <- GetVB.fn(Par = log(xpar), Ages = seq(max(datL[[i]]$Age)),
          Lengths = NULL, ReturnType = "Pred", sdFactor = 1)[, 2]
        names(predL) <- as.character(1:max(datL[[i]]$Age))

        x <- split(datL[[i]]$Length_cm,datL[[i]]$Age_2)
        xsd <- unlist(lapply(x,sd))
        xcv <- xsd/predL[names(xsd)]
        if(is.null(bins)) {ages <- as.numeric(names(xsd))}
        if(!is.null(bins)) {ages <- bins[as.numeric(names(xsd))]}
        out[[i]] <- data.frame(ages=ages,sd=xsd,cv=xcv)

        plot(ages,xsd,xlab="Age",ylab="SD of length-at-age",type="b",pch=16,lty=1,main=names(datL)[i],...)
        par(new = TRUE)
        plot(ages,xcv,xlab="",ylab="",yaxt="n",type="b",pch=3,lty=2,...)
        axis(4)
        mtext("CV",side=4,line=2.6)
        legend(x=legX,y=legY,c("SD","CV"),pch=c(16,3),lty=c(1,2))
        
        xvals <- GetVB.fn(Par = log(xpar), Ages = ages,
          Lengths = NULL, ReturnType = "Pred", sdFactor = 1)[, 2]
        plot(xvals, xsd, xlab = "Predicted length (cm) at age",
          ylab = "SD of length-at-age",
          type = "b", pch = 16, lty = 1, main = names(datL)[i], ...)
        par(new = TRUE)
        plot(xvals, xcv, xlab = "", ylab = "",
          yaxt = "n", type = "b", pch = 3, lty = 2, ...)
        axis(4)
        mtext("CV",side=4,line=2.6)
        legend(x=legX,y=legY,c("SD","CV"),pch=c(16,3),lty=c(1,2))
    }

    #### Save parameter values
    if (dopng) {
      fn.csv <- gsub("\\.png$", "_vbpars\\.csv", fn.png)
      outparsall <- do.call("rbind", outpars)
      colnames(outparsall) <- names(parStart)
      utils::write.csv(outparsall,
        file = fn.csv, row.names = FALSE)
    }
    return(out)
}