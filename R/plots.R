##=============================================================================
##
## Copyright (c) 2018 Paul McKeigue
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##=============================================================================

#' Plot the distribution of the weight of evidence in cases and in controls
#'
#' @param densities Densities object produced by \code{\link{Wdensities}}.
#' @param mask If not \code{NULL}, breaks the y axis to show more detail of
#'        lower end.
#' @param distlabels Character vector of length 2.
#'
#' @return
#' A ggplot object representing the distributions of crude and model-based
#' weights of evidence in cases and in controls.
#'
#' @examples
#' data("cleveland") # load example dataset
#' W <- with(cleveland, weightsofevidence(posterior.p, prior.p))
#' densities <- Wdensities(cleveland$y, W)
#' plotWdists(densities)
#'
#' # Example which requires fitting a mixture distribution
#' data("fitonly")
#' W <- with(fitonly, weightsofevidence(posterior.p, prior.p))
#' densities <- Wdensities(fitonly$y, W, in.spike=W < -2)
#'
#' # truncate spike
#' p <- plotWdists(densities)
#' p + ggplot2::scale_y_continuous(limits=c(0, 0.5), expand=c(0, 0))
#'
#' @importFrom reshape2 melt
#' @export
plotWdists <- function(densities, mask=NULL,
                       distlabels=c("Crude", "Adjusted")) {
    validate.densities(densities)
    dists.data <- data.frame(W=densities$x,
                             Controls=densities$f.ctrls.crude,
                             Cases=densities$f.cases.crude,
                             Controls.adj=densities$f.ctrls,
                             Cases.adj=densities$f.cases)
    dists.long <- melt(dists.data, id="W")
    names(dists.long)[2] <- "status"
    dists.long$adjusted <- ifelse(grepl(".adj", dists.long$status),
                                  distlabels[2], distlabels[1])
    dists.long$status <- gsub(".adj$", "", dists.long$status)
    if(!is.null(mask)) {
        max.value <- max(dists.long$value)
        ceiling.base <- 2
        base <- dists.long$value <= ceiling.base
        floor.peak <- floor(max.value - ceiling.base)
        peak <- dists.long$value > floor.peak
        dists.long$mask <- NA
        dists.long$mask[base] <- 1
        dists.long$mask[peak] <- 0
        print(table(dists.long$mask, exclude=NULL))

        max.value.base <- max(dists.long$value[base])
        min.value.peak <- min(dists.long$value[peak])
        ## rescale peak values to start from 0
        dists.long$value[peak] <- dists.long$value[peak] - floor.peak
        step <- 0.5
        ## breaks and labels must have the same length
        labels.base <- breaks.base <- c(seq(0, max.value.base, by=step))
        labels.peak <- c(seq(0, max.value, by=step))
        breaks.peak <- labels.peak - floor.peak
        labels <- c(labels.base, labels.peak)
        breaks <- c(breaks.base, breaks.peak)
        dists.long <- dists.long[!is.na(dists.long$mask), ]
    }

    expand <- c(0.005, 0.005)
    p <- ggplot(dists.long, aes_(x=quote(tobits(W)), y=~value,
                                linetype=~adjusted, colour=~status)) +
        geom_line(size=1.25) +
        scale_linetype_manual(values=c("dotted", "solid")) +
        scale_color_manual(values=c(Controls='#000000', Cases='#FF0000')) +
        scale_x_continuous(limits=c(min(dists.long$W), max(dists.long$W))) +
        theme_grey(base_size=20) +
        xlab("Weight of evidence case/control (bits)") +
        ylab("Probability density") +
        theme(legend.position=c(0.01, 0.99),
              legend.justification=c(0, 1), # top-left corner of legend box
              legend.title=element_blank()) +
        theme(aspect.ratio=1)

    if(!is.null(mask)) {
        p <- p + scale_y_continuous(breaks=breaks, labels=labels, expand=expand)
    } else {
        p <- p + scale_y_continuous(expand=expand)
    }
    return(p)
}

#' Plot the cumulative frequency distributions in cases and in controls
#'
#' @param densities Densities object produced by \code{\link{Wdensities}}.
#'
#' @return
#' A ggplot object representing the cumulative frequency distributions of the
#' smoothed densities of the weights of evidence in cases and in controls.
#'
#' @examples
#' data("cleveland") # load example dataset
#' W <- with(cleveland, weightsofevidence(posterior.p, prior.p))
#' densities <- Wdensities(cleveland$y, W)
#' plotcumfreqs(densities)
#'
#' @export
plotcumfreqs <- function(densities) {
    validate.densities(densities)
    cumfreqs.ctrls <- data.frame(W=densities$x, F=densities$cumfreq.ctrls,
                                 status="Controls")
    cumfreqs.cases <- data.frame(W=densities$x, F=densities$cumfreq.cases,
                                 status="Cases")
    cumfreqs <- rbind(cumfreqs.ctrls, cumfreqs.cases)
    thresh <- 1e-5
    keep <- which(cumfreqs.ctrls$F + cumfreqs.cases$F > thresh &
                  cumfreqs.ctrls$F + cumfreqs.cases$F < 2 - thresh)
    xlim <- 1.5 * range(cumfreqs.ctrls$W[keep], cumfreqs.cases$W[keep])

    breaks <- seq(0, 1, by=0.1)
    expand <- c(0.005, 0.005)
    p <- ggplot(cumfreqs, aes_(x=quote(tobits(W)), y=~F, colour=~status)) +
        geom_line(size=1.25) +
        scale_color_manual(values=c(Controls='#000000', Cases='#FF0000')) +
        scale_x_continuous(limits=xlim, expand=expand) +
        scale_y_continuous(limits=c(0, 1), breaks=breaks, expand=expand) +
        theme_grey(base_size=20) +
        xlab("Weight of evidence case/control (bits)") +
        ylab("Cumulative probability") +
        theme(legend.position=c(0.99, 0.01),
              legend.justification=c(1, 0), # bottom-right corner of legend box
              legend.title=element_blank()) +
        theme(aspect.ratio=1)
    return(p)
}

#' Plot crude and model-based ROC curves
#'
#' While the crude ROC curve can be non-concave and is generally not smooth,
#' the model-based ROC curve is always concave, as the corresponding densities
#' have been adjusted to be mathematically consistent.
#'
#' @param densities Densities object produced by \code{\link{Wdensities}}.
#'
#' @return
#' A ggplot object representing crude and model-based ROC curves.
#'
#' @examples
#' data("cleveland") # load example dataset
#' W <- with(cleveland, weightsofevidence(posterior.p, prior.p))
#' densities <- Wdensities(cleveland$y, W)
#' plotroc(densities)
#'
#' @importFrom pROC roc
#' @importFrom zoo rollmean
#' @export
plotroc <- function(densities) {
    validate.densities(densities)
    roc.model <- data.frame(x=1 - densities$cumfreq.ctrls,
                            y=1 - densities$cumfreq.cases)
    roc.model$calc <- "Model-based"
    cat("Model-based AUROC", auroc.model(densities), "\n")

    roc.crude <- roc(densities$y, densities$W, direction="<")
    roc.crude <- data.frame(x=1 - roc.crude$specificities,
                            y=roc.crude$sensitivities)
    roc.crude$calc <- "Crude"
    roc <- rbind(roc.model, roc.crude)

    breaks <- seq(0, 1, by=0.1)
    expand <- c(0.005, 0.005)
    p <- ggplot(roc, aes_(x=~x, y=~y, colour=~calc)) +
        geom_line(size=1.25) + coord_fixed() +
        scale_x_continuous(limits=c(0, 1), breaks=breaks, expand=expand) +
        scale_y_continuous(limits=c(0, 1), breaks=breaks, expand=expand) +
        theme_grey(base_size=20) +
        xlab("1 - Specificity") +
        ylab("Sensitivity") +
        theme(legend.position=c(0.99, 0.01),
              legend.justification=c(1, 0), # bottom-right corner of legend box
              legend.title=element_blank())
    return(p)
}

#' Plot log case/control density ratio against weight of evidence as a check that
#' the densities are mathematically consistent
#'
#' @param densities Densities object produced by \code{\link{Wdensities}}.
#'
#' @return
#' A ggplot object representing a plot of the natural log case/control density
#' ratio against the weight of evidence (should be a straight line of gradient
#' 1 passing through the origin).
#'
#' @examples
#' data("cleveland") # load example dataset
#' W <- with(cleveland, weightsofevidence(posterior.p, prior.p))
#' densities <- Wdensities(cleveland$y, W)
#' plotW(densities)
#' 
#' @export
plotW <- function(densities) {
    validate.densities(densities)
    densities.logratio <- log(densities$f.cases / densities$f.ctrls)
    wratios <- data.frame(Wdens=densities$x, Wratio=densities.logratio)
    axislimits <- 1.5 * range(densities$W)
    p <- ggplot(wratios, aes_(x=~Wdens, y=~Wratio)) +
        geom_line(size=1.25) + coord_fixed() +
        scale_x_continuous(limits=axislimits, expand=c(0, 0)) +
        scale_y_continuous(limits=axislimits, expand=c(0, 0)) +
        theme_grey(base_size=20) +
        xlab("Weight of evidence case/control (bits)") +
        ylab("Log ratio case density to control density")
    return(p)
}
