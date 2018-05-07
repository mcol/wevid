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
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
plotWdists <- function(Wdensities.unadj, Wdensities.adj, mask=NULL,
                       distlabels=c("Unadjusted", "Adjusted")) {

    dists.data <- data.frame(W=Wdensities.unadj$x,
                             Controls=Wdensities.unadj$f.ctrls,
                             Cases=Wdensities.unadj$f.cases,
                             Controls.adj=Wdensities.adj$f.ctrls,
                             Cases.adj=Wdensities.adj$f.cases)
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
    p <- ggplot(dists.long, aes(x=log2(exp(1))*dists.long$W, y=value,
                                linetype=adjusted, colour=status)) +
        geom_line(size=1.25) +
        scale_linetype_manual(values=c("dotted", "solid")) +
        scale_color_manual(values=c(Controls='#000000', Cases='#FF0000')) +
        scale_x_continuous(limit=c(min(dists.long$W), max(dists.long$W))) +
        scale_y_continuous() +
        theme_grey(base_size = 20) +
        xlab("Weight of evidence case/control (bits)") +
        ylab("Probability density") +
        theme(legend.position=c(0.01, 0.99),
              legend.justification=c(0, 1), # top-left corner of legend box
              legend.title = element_blank()) +
        theme(aspect.ratio = 1)

    if(!is.null(mask)) {
        p <- p + scale_y_continuous(breaks=breaks, labels=labels, expand=c(0, 0))
    } else {
        p <- p + scale_y_continuous(expand=c(0, 0))
    }
    return(p)
}
