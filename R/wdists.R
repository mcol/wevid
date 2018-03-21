#!/usr/bin/Rscript

library(reshape2)
library(ggplot2)
library(pROC)

## gaussian kernel for density estimation
fsmooth <- function(x, X, h, n) {# scalar x, vector X, scalar h, scalar n
    ## n is length of vector X
    ## density at x is weighted average of observed values
    ## with weights scaling with exp(-t^2)
    return(sum(dnorm((x - X) / h))/ (n * h))
}

kl <- function(p, q) {#KL divergence of p from q
    kl <- p * log( p / q)
    kl[p==0] <- 0
    kl <- mean(kl)
    return(kl)
}

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

weightsofevidence <- function(posterior.p, prior.p) {# weights of evidence in nat log units
 W <- (log(posterior.p) - log(1 - posterior.p) -
                     log(prior.p / (1 - prior.p)))
 return(W)
}

######### these functions take raw W values as arguments ######################

Wdensities.unadjusted <- function(y, W, range.xseq=c(-25, 25), x.stepsize=0.05, adjust.bw=1) {
    n.ctrls <- length(W[y==0])
    n.cases <- length(W[y==1])
    xseq <- seq(range.xseq[1], range.xseq[2], by=x.stepsize)
    fhat.cases.raw <- fhat.ctrls.raw <- numeric(length(xseq))    
    for(i in 1:length(xseq)) {
        fhat.ctrls.raw[i] <- fsmooth(xseq[i], W[y==0], h=adjust.bw * bw.SJ(W[y==0]), n.ctrls)
        fhat.cases.raw[i] <- fsmooth(xseq[i], W[y==1], h=adjust.bw * bw.SJ(W[y==1]), n.cases)
    }
    return(data.frame(x=xseq, f.ctrls=fhat.ctrls.raw, f.cases=fhat.cases.raw))
}

error.fhat <- function(theta, w, range.xseq, x.stepsize, adjust.bw=1) {
     xseq <- seq(range.xseq[1], range.xseq[2], by=x.stepsize)
     weights <- ifelse(w < 0, exp(-theta * w^2), 1) # downweight negative values
     weights <- weights / sum(weights)
     fhat <- density(x=w, weights=weights, bw="SJ", adjust=adjust.bw,
                     from=range.xseq[1], to=range.xseq[2],
                     n=length(xseq))$y
     integral.fhat <- sum(exp(-xseq) * fhat * x.stepsize)
     return(abs(log(integral.fhat)))
}

weights.flipnormalized <- function(w, range.xseq, x.stepsize, adjust.bw=1) {
    ## reweight outliers so that smoothed exp(-x) * p(x) integrates to 1
    theta <- optimize(f=error.fhat, interval=c(-1, 1), w=w,
                      range.xseq=range.xseq, x.stepsize=x.stepsize, adjust.bw)
    weights <- ifelse(w < 0, exp(-theta$minimum * w^2), 1)
    weights <- weights / sum(weights)
    return(weights)
}

error.fhat2 <- function(par, w, range.xseq, x.stepsize, adjust.bw=1) {# returns objective function
    theta <- par[1]
    mu <- par[2]
    xseq <- seq(range.xseq[1], range.xseq[2], by=x.stepsize)
    weights <- exp(-theta * (w - mu)^2)
    weights <- weights / sum(weights)
    fhat <- density(w, bw=bw.SJ(w), adjust=adjust.bw, weights=weights,
                    n=length(xseq), from=min(xseq), to=max(xseq))$y
    #integral.fhat <- sum(fhat * x.stepsize)
    integral.fhat.flipped <- sum(exp(-xseq) * fhat * x.stepsize)
    mean.fhat <- sum(xseq * fhat * x.stepsize)
    return(abs(log(integral.fhat.flipped)) + 0.1 * abs(mean.fhat - mean(w)))
}

weights2.flipnormalized <- function(w, range.xseq, x.stepsize, adjust.bw=adjust.bw) {
    ## reweight outliers so that smoothed exp(-x) * p(x) integrates to 1
    optim.result <- optim(par=c(0, mean(w)), fn=error.fhat2, w=w,
                       range.xseq=range.xseq, x.stepsize=x.stepsize,
                       adjust.bw=adjust.bw)
    theta <- optim.result$par[1]
    mu <- optim.result$par[2]
    weights <- exp(-theta * (w - mu)^2)
    weights <- weights / sum(weights)
    return(weights)
}

Wdensities.reweighted <- function(y, W, range.xseq, x.stepsize, adjust.bw) {
    xseq <- seq(range.xseq[1], range.xseq[2], by=x.stepsize)

    weights.cases <- weights2.flipnormalized(w=W[y==1], range.xseq, x.stepsize,
                                             adjust.bw=adjust.bw)
    fhat.cases <- density(W[y==1], weights=weights.cases, bw="SJ", adjust=adjust.bw,
                          from=range.xseq[1],
                          to=range.xseq[2], n=length(xseq))$y

    ## for controls, reverse sign of w to calculate weights
    weights.ctrls <- weights2.flipnormalized(w=-W[y==0], range.xseq, x.stepsize, adjust.bw)
    fhat.ctrls <- density(W[y==0], weights=weights.ctrls, bw="SJ", adjust=adjust.bw,
                          from=range.xseq[1],
                          to=range.xseq[2], n=length(xseq))$y
    ctrls.integral <- sum(fhat.ctrls * x.stepsize)
    cases.integral <- sum(fhat.cases * x.stepsize)
    cat("reweighted density in ctrls normalizes to", ctrls.integral,
        ", flipped density to",  sum(exp(xseq) * fhat.ctrls * x.stepsize), "\n")
    cat("reweighted density in cases normalizes to", cases.integral,
        ", flipped density to",  sum(exp(-xseq) * fhat.cases * x.stepsize), "\n")
    if(abs(log(ctrls.integral)) > 0.05 | abs(log(cases.integral)) > 0.05) {
        stop("reweighted density failed to normalize: try value lower than ",
             x.stepsize, " for x.stepsize\n")
    }
    return(data.frame(x=xseq, f.ctrls.reweighted=fhat.ctrls,
                      f.cases.reweighted=fhat.cases))
}


Wdensities.adjusted <- function(y, W, range.xseq=c(-25, 25), x.stepsize=0.05,
                                adjust.bw=1) {# reweights raw observations density from raw W values, then averages direct and flipped
    densities.reweighted <- Wdensities.reweighted(y, W, range.xseq, x.stepsize,
                                                  adjust.bw=adjust.bw)
 
    n.ctrls <- length(W[y==0])
    n.cases <- length(W[y==1])
    xseq <- seq(range.xseq[1], range.xseq[2], by=x.stepsize)

    wts <- c(n.ctrls, n.cases) / (n.ctrls + n.cases)
    ## weights proportional to num controls, num cases
    wt.exp <- TRUE

    if(wt.exp) { # weighted average of direct and flipped densities
        wts.xseq <- cbind(exp(-xseq), 1)
        wts.xseq <- wts.xseq / rowSums(wts.xseq) # neg values of xseq weight col 1, pos values weight col 2
        ## reweight col1/col2 by ctrls/cases
        wts.xseq.cases <- t(wts * t(wts.xseq))
        ## for average density in controls and flipped density in cases, swap cols of wts.seq and reweight col1/col2 by ctrls/cases
        wts.xseq.swap <- cbind(wts.xseq[, 2], wts.xseq[, 1])
        wts.xseq.ctrls <- t(wts * t(wts.xseq.swap))
        wts.xseq.cases <- wts.xseq.cases / rowSums(wts.xseq.cases)
        wts.xseq.ctrls <- wts.xseq.ctrls / rowSums(wts.xseq.ctrls)
    } else { # unweighted average of direct and flipped densities
        wts.xseq.cases <- wts.xseq.ctrls <- matrix(rep(wts, each=length(xseq)), ncol=2)
    }   
        
    fhat.ctrls.fromcases <- exp(-xseq) * densities.reweighted$f.cases.reweighted
    fhat.cases.fromctrls <- exp(xseq) * densities.reweighted$f.ctrls.reweighted

    ## average these densities
    fhat.ctrls.average <- average.densities(densities.reweighted$f.ctrls.reweighted,
                                            fhat.ctrls.fromcases,
                                            wts.xseq.ctrls, x.stepsize)
    fhat.cases.average <- average.densities(fhat.cases.fromctrls,
                                            densities.reweighted$f.cases.reweighted,
                                            wts.xseq.cases, x.stepsize)
    ## check normalization and agreement
    ctrls.integral <- sum(exp(xseq) * fhat.ctrls.average * x.stepsize)
    cases.integral <- sum(exp(-xseq) * fhat.cases.average * x.stepsize)
    cat("flipped averaged density in ctrls normalizes to", ctrls.integral, "\n")
    cat("flipped averaged density in cases normalizes to", cases.integral, "\n")

    if(abs(log(ctrls.integral)) > 0.1 | abs(log(cases.integral)) > 0.1) {
        cat("Averaged densities would not normalize\n") 
        fhat.ctrls.adjusted <- fhat.ctrls
        fhat.cases.adjusted <- fhat.cases
    } else {
       fhat.ctrls.adjusted <- fhat.ctrls.average
       fhat.cases.adjusted <- fhat.cases.average
    }
    
    cat("KL divergence of density of e^-W p(W) in cases from p(W) in controls",
        kl(fhat.ctrls.adjusted, exp(-xseq) * fhat.cases.adjusted), "\n")
    return(data.frame(x=xseq, f.ctrls=fhat.ctrls.adjusted, f.cases=fhat.cases.adjusted))
}

########### these functions take densities as arguments  ###############################

error.dens2 <- function(par, f, range.xseq, x.stepsize) {# returns objective function for finding optimal weights
    theta <- par[1]
    mu <- par[2]
    xseq <- seq(range.xseq[1], range.xseq[2], by=x.stepsize)
    weights <- exp(-theta * (xseq - mu)^2)
    weights <- weights / sum(weights * f * x.stepsize) # normalize integral
    integral.f <- sum(weights * f * x.stepsize)
    integral.f.flipped <- sum(weights * exp(-xseq) * f * x.stepsize)
    mean.f <- sum(weights * xseq * f * x.stepsize)
    mean.unweighted <- sum(xseq * f * x.stepsize)
#    cat("theta", theta, "mu", mu, "integral", integral.f,
#        "integral", integral.f.flipped, "meandiff", mean.f - mean.unweighted, "\n")
    return(abs(log(integral.f.flipped)) + 0.1 * abs(mean.f - mean.unweighted))
}

weights2.flipdens <- function(f, range.xseq, x.stepsize) {#returns optimal weights so that flipped density integrates to 1
    xseq <- seq(range.xseq[1], range.xseq[2], by=x.stepsize)
    mean.unweighted <- sum(xseq * f * x.stepsize)
    ## reweight outliers so that smoothed exp(-x) * p(x) integrates to 1
    optim.result <- optim(par=c(0, mean.unweighted), fn=error.dens2, f=f,
                       range.xseq=range.xseq, x.stepsize=x.stepsize)
    theta <- optim.result$par[1]
    mu <- optim.result$par[2]
    weights <- exp(-theta * (xseq - mu)^2)
    weights <- weights / sum(weights * f * x.stepsize)
    return(weights)
}

fdensities.reweighted <- function(densities, range.xseq, x.stepsize) {# reweights densities so that flipped density integrates to 1
    xseq <- seq(range.xseq[1], range.xseq[2], by=x.stepsize)
    
    weights.cases <- weights2.flipdens(densities$f.cases, range.xseq, x.stepsize)
    fhat.cases <- weights.cases * densities$f.cases

    ## for controls, reverse sign of w to calculate weights
    weights.ctrls <- weights2.flipdens(rev(densities$f.ctrls), range.xseq, x.stepsize)
    fhat.ctrls <- weights.ctrls * rev(densities$f.ctrls)
    fhat.ctrls <- rev(fhat.ctrls)

    ctrls.integral <- sum(fhat.ctrls * x.stepsize)
    cases.integral <- sum(fhat.cases * x.stepsize)
    cat("reweighted density in ctrls normalizes to", ctrls.integral,
        ", flipped density to",  sum(exp(xseq) * fhat.ctrls * x.stepsize), "\n")
    cat("reweighted density in cases normalizes to", cases.integral,
        ", flipped density to",  sum(exp(-xseq) * fhat.cases * x.stepsize), "\n")
    if(abs(log(ctrls.integral)) > 0.05 | abs(log(cases.integral)) > 0.05) {
        stop("reweighted density failed to normalize\n")
    }
    return(data.frame(x=xseq, f.ctrls.reweighted=fhat.ctrls,
                      f.cases.reweighted=fhat.cases))
}

average.densities <- function(f1, f2, wts, x.stepsize) { # averages densities given weights
    f.average <- wts[, 1] * f1 + wts[, 2] * f2
    f.average <- f.average / sum(f.average * x.stepsize) #renormalize
    return(f.average)
}

Wdensities.reweighted.average <- function(f.cases, f.ctrls, n.ctrls, n.cases,
                                          range.xseq=c(-25, 25),
                                          x.stepsize=0.05) {# averages direct and flipped
    xseq <- seq(range.xseq[1], range.xseq[2], by=x.stepsize)

    wts <- c(n.ctrls, n.cases) / (n.ctrls + n.cases)
    ## weights proportional to num controls, num cases
    wt.exp <- TRUE

    if(wt.exp) { # weighted average of direct and flipped densities
        wts.xseq <- cbind(exp(-xseq), 1)
        wts.xseq <- wts.xseq / rowSums(wts.xseq) # neg values of xseq weight col 1, pos values weight col 2
        ## reweight col1/col2 by ctrls/cases
        wts.xseq.cases <- t(wts * t(wts.xseq))
        ## for average density in controls and flipped density in cases, swap cols of wts.seq and reweight col1/col2 by ctrls/cases
        wts.xseq.swap <- cbind(wts.xseq[, 2], wts.xseq[, 1])
        wts.xseq.ctrls <- t(wts * t(wts.xseq.swap))
        wts.xseq.cases <- wts.xseq.cases / rowSums(wts.xseq.cases)
        wts.xseq.ctrls <- wts.xseq.ctrls / rowSums(wts.xseq.ctrls)
    } else { # unweighted average of direct and flipped densities
        wts.xseq.cases <- wts.xseq.ctrls <- matrix(rep(wts, each=length(xseq)), ncol=2)
    }   
        
    fhat.ctrls.fromcases <- exp(-xseq) * f.cases
    fhat.cases.fromctrls <- exp(xseq) * f.ctrls

    ## average these densities
    fhat.ctrls.average <- average.densities(f.ctrls,
                                            fhat.ctrls.fromcases,
                                            wts.xseq.ctrls, x.stepsize)
    fhat.cases.average <- average.densities(fhat.cases.fromctrls,
                                            f.cases,
                                            wts.xseq.cases, x.stepsize)
    ## check normalization and agreement
    ctrls.integral <- sum(exp(xseq) * fhat.ctrls.average * x.stepsize)
    cases.integral <- sum(exp(-xseq) * fhat.cases.average * x.stepsize)
    cat("flipped averaged density in ctrls normalizes to", ctrls.integral, "\n")
    cat("flipped averaged density in cases normalizes to", cases.integral, "\n")

    if(abs(log(ctrls.integral)) > 0.1 | abs(log(cases.integral)) > 0.1) {
        cat("Averaged densities would not normalize\n") 
        fhat.ctrls.adjusted <- fhat.ctrls
        fhat.cases.adjusted <- fhat.cases
    } else {
       fhat.ctrls.adjusted <- fhat.ctrls.average
       fhat.cases.adjusted <- fhat.cases.average
    }
    
    cat("KL divergence of density of e^-W p(W) in cases from p(W) in controls",
        kl(fhat.ctrls.adjusted, exp(-xseq) * fhat.cases.adjusted), "\n")
    return(data.frame(x=xseq, f.ctrls=fhat.ctrls.adjusted, f.cases=fhat.cases.adjusted))
}

##################################################################

plotWdists <- function(Wdensities.unadj, Wdensities.adj,
                       distlabels=c("Adjusted", "Unadjusted")) {
    dists.data <- data.frame(W=Wdensities.unadj$x,
                             Controls=Wdensities.unadj$f.ctrls,
                             Cases=Wdensities.unadj$f.cases,
                             Controls.adj=Wdensities.adj$f.ctrls,
                             Cases.adj=Wdensities.adj$f.cases)
    dists.long <- melt(dists.data, id="W")
    names(dists.long)[2] <- "status"
    dists.long$adjusted <- ifelse(grepl(".adj", dists.long$status),
                                  distlabels[1], distlabels[2])
    dists.long$status <- gsub(".adj$", "", dists.long$status)
    p <- ggplot(dists.long, aes(x=log2(exp(1))*W, y=value,
                                linetype=adjusted, colour=status)) +
        geom_line(size=1.25) +
        scale_linetype_manual(values=c("dotted", "solid")) +
        scale_color_manual(values=c(Controls='#000000', Cases='#FF0000')) +
        scale_x_continuous(limit=2 * c(min(W), max(W))) + 
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

##############################################################################
