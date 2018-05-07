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

#' Gaussian kernel for density estimation
#'
#' @param x scalar.
#' @param X vector.
#' @param h Bandwidth size.
#' @param n Length of vector \var{X}.
#' @keywords internal
fsmooth <- function(x, X, h, n) {
    ## density at x is weighted average of observed values
    ## with weights scaling with exp(-t^2)
    return(sum(dnorm((x - X) / h))/ (n * h))
}

#' Kullback-Leibler divergence of p from q
#'
#' @param p Probability distribution.
#' @param q Probability distribution.
#' @keywords internal
kl <- function(p, q) {
    kl <- p * log( p / q)
    kl[p==0] <- 0
    kl <- mean(kl)
    return(kl)
}

#' @export
wtrue.results <- function(studyname, y, posterior.p, prior.p) {
    auroc <- round(auc(y, posterior.p), 3)

    ## weight of evidence in favour of true status
    loglikrat <- (2 * y - 1) * (log(posterior.p) - log(1 - posterior.p) -
                                      log(prior.p / (1 - prior.p)))
    loglikrat.case <- loglikrat[y==1]
    loglikrat.ctrl <- loglikrat[y==0]
    ## test log-likelihood as difference from prior log-likelihood
    loglik <- y * log(posterior.p) + (1 - y) * log(1 - posterior.p) -
        (y * log(prior.p) + (1 - y) * log(1 - prior.p))

    results <- data.frame(model=studyname,
                          casectrlrat <- paste0(as.character(length(y[y==1])), " / ",
                                                as.character(length(y[y==0]))),
                          auroc=round(auroc, 3),
                          loglikrat.all=round(log2(exp(1)) * mean(loglikrat), 2), 
                          loglikrat.varmeanrat=round(var(loglikrat) / mean(loglikrat), 2),
                          loglikrat.case= round(log2(exp(1)) * mean(loglikrat.case), 2), 
                          loglikrat.ctrl=round(log2(exp(1)) * mean(loglikrat.ctrl), 2), 
                          test.loglik=round(log2(exp(1)) * sum(loglik), 2)
                          ) 
    names(results) <-
        c("Model", "Cases / controls",
          "C-statistic",
          "Average weight of evidence (bits)",
          "W variance / mean ratio (nats)",
          paste("Average W in cases (bits)"),
          paste("Average W in controls (bits)"),
          "Test log-likelihood (bits)")
    return(results)
}

#' Weights of evidence in nat log units
#' @export
weightsofevidence <- function(posterior.p, prior.p) {
 W <- (log(posterior.p) - log(1 - posterior.p) -
                     log(prior.p / (1 - prior.p)))
 return(W)
}

######### these functions take raw W values as arguments ######################

#' Calculate the unadjusted smoothed densities of W in cases and in controls
#'
#' @param y Binary outcome label (0 for controls, 1 for cases).
#' @param W Weight of evidence.
#' @param range.xseq Range of points where the curves should be sampled.
#' @param x.stepsize Distance between each point.
#' @param adjust.bw Bandwidth adjustment.
#'
#' @export
Wdensities.unadjusted <- function(y, W, range.xseq=c(-25, 25), x.stepsize=0.05,
                                  adjust.bw=1) {
    n.ctrls <- sum(y == 0)
    n.cases <- sum(y == 1)
    if (n.ctrls + n.cases != length(y))
        stop("y contains values different from 0 or 1")
    xseq <- seq(range.xseq[1], range.xseq[2], by=x.stepsize)
    fhat.cases.raw <- fhat.ctrls.raw <- numeric(length(xseq))
    W.ctrls <- W[y == 0]
    W.cases <- W[y == 1]
    bw.ctrls <- bw.SJ(W.ctrls) * adjust.bw
    bw.cases <- bw.SJ(W.cases) * adjust.bw
    for(i in 1:length(xseq)) {
        fhat.ctrls.raw[i] <- fsmooth(xseq[i], W.ctrls, h=bw.ctrls, n.ctrls)
        fhat.cases.raw[i] <- fsmooth(xseq[i], W.cases, h=bw.cases, n.cases)
    }
    return(data.frame(x=xseq, f.ctrls=fhat.ctrls.raw, f.cases=fhat.cases.raw))
}



##################################################################

#' @export
plotWdists <- function(Wdensities.unadj, Wdensities.adj,
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
        theme(legend.position=c(0.8, 0.7),
              legend.title = element_blank()) + 
            theme(aspect.ratio = 1)
    return(p)
}

density.spike.slab <- function(w, in.spike, range.xseq, x.stepsize) {
    xseq <- seq(range.xseq[1], range.xseq[2], by=x.stepsize)

    density.spike <- density(w[in.spike], bw="SJ", n=length(xseq),
                                    from=min(xseq), to=max(xseq))
    density.slab <- density(w[!in.spike], bw="SJ", n=length(xseq),
                                   from=min(xseq), to=max(xseq))
    wts.mix <- as.integer(table(in.spike))
    wts.mix <- wts.mix / sum(wts.mix)
    density.mix <- data.frame(x=xseq,
                              y=wts.mix[1] * density.slab$y +
                                  wts.mix[2] * density.spike$y)
    return(density.mix)
}

Wdensities.mix <- function(W, yobs, in.spike, range.xseq=c(-25, 25), x.stepsize=0.01) {
    xseq <- seq(range.xseq[1], range.xseq[2], by=x.stepsize)

    Wdensity.mix.ctrls <- density.spike.slab(W[yobs==0], in.spike[yobs==0],
                                             range.xseq, x.stepsize)
    Wdensity.mix.cases <- density.spike.slab(W[yobs==1], in.spike[yobs==1],
                                             range.xseq, x.stepsize)
    return(data.frame(x=xseq, f.ctrls=Wdensity.mix.ctrls$y,
                      f.cases=Wdensity.mix.cases$y))
}
