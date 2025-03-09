#' Causal inference via general estimands on binary outcomes: the weighted ATE, weighted ATT and weighted ATC
#' @param A vector of the treatment variable (binary, valued from 0 and 1)
#' @param Y vector of the binary outcome variable
#' @param X matrix of covariates/confounders that are included into the propensity score model
#' @param beta whether to include estimands from beta family weights; if so, the v parameter below is required to be specified
#' @param v model parameter for beta weights
#' @param trim whether to include trimming; if so, the trim.alpha parameter below needs to be specified
#' @param trim.alpha the trimming threshold; default is .05
#' @param trun whether to include truncation; if so, the trun.alpha parameter below needs to be specified
#' @param trun.alpha the truncation threshold; default is .05
#' @param boot whether to implement bootstrap for variance estimation; default is TRUE
#' @param n.boot number of bootstrap, if boot==TRUE; default is 500
#' @param return.psfig whether to plot and return the estimated propensity score distribution plots by treatment group
#' @param seed seed for initializing the random number generator, which is used in set.seed() function; default is 4399
#' @param conf.level level of confidence interval; default is .05 (for 0.95 confidence interval, using quantile method in bootstrap)
#' @return the function returns estimated causal estimands in the form of weighted risk difference (RD), risk ratio (RR) and odds ratio (OR);
#'         estimands defined by constant weights, overlap weights, matching weights and entropy weights will always be returned;
#'         whether to return estimands defined by trimming, truncation, and beta weights can be specified by user
ChiPS_bin <- function(A, Y, X,
                      ps.library="SL.glm",
                      beta=FALSE, v=NA,
                      trim=FALSE, trim.alpha=.05,
                      trun=FALSE, trun.alpha=.05,
                      boot=TRUE, n.boot=500,
                      seed=4399,
                      return.psfig=FALSE,
                      conf.level=.05) {

  if(sum(sort(unique(Y))!=c(0,1))!=0) { stop("The outcome must be binary in this function!") }

  set.seed(seed=seed)
  n <- length(A)
  X <- as.data.frame(X)
  e.h <- .propensity(A=A, Y=Y, X=X, X.pred=X, ps.library=ps.library)

  if(return.psfig) {
    df <- data.frame(A=A, e.h=e.h)
    ps.hist <- ggplot(df, aes(x=e.h, fill=factor(A))) +
      geom_histogram(position="identity", alpha=0.35, color="black", bins=35) +
      labs(x="Estimated propensity score", y="Count") +
      scale_fill_manual(values=c("darkblue", "coral"), name="Group") +
      theme_minimal()

    ps.density <- ggplot(df, aes(x=e.h, fill=factor(A))) +
      geom_density(alpha=0.35) +
      labs(x="Estimated propensity score", y="Count") +
      scale_fill_manual(values=c("darkblue", "coral"), name="Group") +
      theme_minimal()
  }

  result <- .point.est(A=A, Y=Y, X=X,
                       ps.library=ps.library,
                       beta=beta, v=v,
                       trim=trim, trim.alpha=trim.alpha,
                       trun=trun, trun.alpha=trun.alpha)

  if(!boot) { return(result) }
  if(boot) {
    watc.all.RD <- watt.all.RD <- wate.all.RD <- NULL
    watc.all.RR <- watt.all.RR <- wate.all.RR <- NULL
    watc.all.OR <- watt.all.OR <- wate.all.OR <- NULL
    for(i in 1:n.boot) {
      bt.samp <- sample(1:n, n, replace=T)
      A.bt <- A[bt.samp]
      Y.bt <- Y[bt.samp]
      X.bt <- X[bt.samp, ]

      all.result <- .point.est(A=A.bt, Y=Y.bt, X=X.bt,
                               ps.library=ps.library,
                               beta=beta, v=v,
                               trim=trim, trim.alpha=trim.alpha,
                               trun=trun, trun.alpha=trun.alpha)

      wate.all.RD <- rbind(all.result$RD.wate, wate.all.RD)
      watt.all.RD <- rbind(all.result$RD.watt, watt.all.RD)
      watc.all.RD <- rbind(all.result$RD.watc, watc.all.RD)

      wate.all.RR <- rbind(all.result$RR.wate, wate.all.RR)
      watt.all.RR <- rbind(all.result$RR.watt, watt.all.RR)
      watc.all.RR <- rbind(all.result$RR.watc, watc.all.RR)

      wate.all.OR <- rbind(all.result$OR.wate, wate.all.OR)
      watt.all.OR <- rbind(all.result$OR.watt, watt.all.OR)
      watc.all.OR <- rbind(all.result$OR.watc, watc.all.OR)

      if(i%%50==0) print(paste0("Bootstrap ", i, " is done."))
    }
    quant <- 1-conf.level/2

    upr.wate.RD <- apply(wate.all.RD, 2, quantile, quant)
    upr.watt.RD <- apply(watt.all.RD, 2, quantile, quant)
    upr.watc.RD <- apply(watc.all.RD, 2, quantile, quant)

    lwr.wate.RD <- apply(wate.all.RD, 2, quantile, 1-quant)
    lwr.watt.RD <- apply(watt.all.RD, 2, quantile, 1-quant)
    lwr.watc.RD <- apply(watc.all.RD, 2, quantile, 1-quant)

    upr.wate.RR <- apply(wate.all.RR, 2, quantile, quant)
    upr.watt.RR <- apply(watt.all.RR, 2, quantile, quant)
    upr.watc.RR <- apply(watc.all.RR, 2, quantile, quant)

    lwr.wate.RR <- apply(wate.all.RR, 2, quantile, 1-quant)
    lwr.watt.RR <- apply(watt.all.RR, 2, quantile, 1-quant)
    lwr.watc.RR <- apply(watc.all.RR, 2, quantile, 1-quant)

    upr.wate.OR <- apply(wate.all.OR, 2, quantile, quant)
    upr.watt.OR <- apply(watt.all.OR, 2, quantile, quant)
    upr.watc.OR <- apply(watc.all.OR, 2, quantile, quant)

    lwr.wate.OR <- apply(wate.all.OR, 2, quantile, 1-quant)
    lwr.watt.OR <- apply(watt.all.OR, 2, quantile, 1-quant)
    lwr.watc.OR <- apply(watc.all.OR, 2, quantile, 1-quant)

    return(list(df.WATE.RD=data.frame(Est=result$RD.wate, Upr=upr.wate.RD, Lwr=lwr.wate.RD),
                df.WATT.RD=data.frame(Est=result$RD.watt, Upr=upr.watt.RD, Lwr=lwr.watt.RD),
                df.WATC.RD=data.frame(Est=result$RD.watc, Upr=upr.watc.RD, Lwr=lwr.watc.RD),

                df.WATE.RR=data.frame(Est=result$RR.wate, Upr=upr.wate.RR, Lwr=lwr.wate.RR),
                df.WATT.RR=data.frame(Est=result$RR.watt, Upr=upr.watt.RR, Lwr=lwr.watt.RR),
                df.WATC.RR=data.frame(Est=result$RR.watc, Upr=upr.watc.RR, Lwr=lwr.watc.RR),

                df.WATE.OR=data.frame(Est=result$OR.wate, Upr=upr.wate.OR, Lwr=lwr.wate.OR),
                df.WATT.OR=data.frame(Est=result$OR.watt, Upr=upr.watt.OR, Lwr=lwr.watt.OR),
                df.WATC.OR=data.frame(Est=result$OR.watc, Upr=upr.watc.OR, Lwr=lwr.watc.OR)))
  }
}

.point.est <- function(A, Y, X,
                       ps.library="SL.glm",
                       beta=FALSE, v=NA,
                       trim=FALSE, trim.alpha=.05,
                       trun=FALSE, trun.alpha=.05) {

  # propensity score estimation
  n <- length(A)
  e.h <- .propensity(A=A, Y=Y, X=X, X.pred=X, ps.library=ps.library)

  # define the tilting functions
  weights <- c("overall", "treated", "control", "overlap", "matching", "entropy")
  tilt.wate <- data.frame(1, e.h, 1-e.h, e.h*(1-e.h), pmin(e.h, 1-e.h), -e.h*log(e.h)-(1-e.h)*log(1-e.h) )
  tilt.watt <- data.frame(1, NA, NA, e.h*(1-e.h)*(1-A) + A, pmin(e.h, 1-e.h)*(1-A) + A, (-e.h*log(e.h)-(1-e.h)*log(1-e.h))*(1-A) + A )
  tilt.watc <- data.frame(1, NA, NA, e.h*(1-e.h)*A + (1-A), pmin(e.h, 1-e.h)*A + (1-A), (-e.h*log(e.h)-(1-e.h)*log(1-e.h))*A + (1-A) )

  if(beta) {
    v <- v[!is.na(v)]
    if(length(v)>0) {
      weights <- c(weights, paste0("beta (v=", v, ")"))
      tilt.beta.e <- tilt.beta.t <- tilt.beta.c <- matrix(NA, nrow=n, ncol=length(v))
      v <- v-1
      for(i in 1:length(v)) {
        tilt.beta.e[,i] <- e.h^v[i]*(1-e.h)^v[i]
        tilt.beta.t[,i] <- e.h^v[i]*(1-e.h)^v[i]*(1-A) + A
        tilt.beta.c[,i] <- e.h^v[i]*(1-e.h)^v[i]*A + (1-A)
      }
      tilt.wate <- cbind(tilt.wate, tilt.beta.e)
      tilt.watt <- cbind(tilt.watt, tilt.beta.t)
      tilt.watc <- cbind(tilt.watc, tilt.beta.c)
    }
  }

  if(trim) {
    trim.alpha <- trim.alpha[!is.na(trim.alpha)]
    if(length(trim.alpha)>0) {
      weights <- c(weights, paste0("trimming (alpha=", trim.alpha, ")"))
      tilt.trim.wate <- tilt.trim.watt <- tilt.trim.watc <-
        matrix(NA, nrow=n, ncol=length(trim.alpha))
      for(i in 1:length(trim.alpha)) {
        tilt.trim.wate[,i] <- as.numeric(e.h<1-trim.alpha[i] & e.h>trim.alpha[i])
        tilt.trim.watt[,i] <- as.numeric(e.h<1-trim.alpha[i])*(1-A) + A
        tilt.trim.watc[,i] <- as.numeric(e.h>trim.alpha[i])*A + (1-A)
      }
      tilt.wate <- cbind(tilt.wate, tilt.trim.wate)
      tilt.watt <- cbind(tilt.watt, tilt.trim.watt)
      tilt.watc <- cbind(tilt.watc, tilt.trim.watc)
    }
  }

  if(trun) {
    trun.alpha <- trun.alpha[!is.na(trun.alpha)]
    if(length(trun.alpha)>0) {
      weights <- c(weights, paste0("truncation (alpha=", trun.alpha, ")"))
      tilt.trun.wate <- tilt.trun.watt <- tilt.trun.watc <-
        matrix(NA, nrow=n, ncol=length(trun.alpha))
      for(i in 1:length(trun.alpha)) {
        tilt.trun.wate[,i] <-
          as.numeric(e.h<1-trun.alpha[i] & e.h>trun.alpha[i]) +
          ( as.numeric(e.h<trun.alpha[i])/trun.alpha[i]*e.h + as.numeric(e.h<trun.alpha[i])/(1-trun.alpha[i])*(1-e.h) )
        tilt.trun.watt[,i] <-
          (as.numeric(e.h<1-trun.alpha[i]) + as.numeric(e.h>1-trun.alpha[i])*(trun.alpha[i]/(1-trun.alpha[i]) * (1-e.h)/e.h)) * (1-A) + A
        tilt.trun.watc[,i] <-
          (as.numeric(e.h>trun.alpha[i]) + as.numeric(e.h<trun.alpha[i])*((1-trun.alpha[i])/trun.alpha[i] * e.h/(1-e.h))) * A + (1-A)
      }
      tilt.wate <- cbind(tilt.wate, tilt.trun.wate)
      tilt.watt <- cbind(tilt.watt, tilt.trun.watt)
      tilt.watc <- cbind(tilt.watc, tilt.trun.watc)
    }
  }
  colnames(tilt.wate) <- colnames(tilt.watt) <- colnames(tilt.watc) <- weights

  # calculate the propensity score weighted estimators
  M <- ncol(tilt.wate)
  IPW.ate <- matrix(A/e.h - (1-A)/(1-e.h), ncol=M, nrow=n)
  IPW.att <- matrix(A - (1-A)*e.h/(1-e.h), ncol=M, nrow=n)
  IPW.atc <- matrix(A*(1-e.h)/e.h - (1-A), ncol=M, nrow=n)
  Y.mat <- matrix(Y, ncol=M, nrow=n)
  A.mat <- matrix(A, ncol=M, nrow=n)

  p1h <- apply(A.mat*IPW.ate*tilt.wate*Y.mat, 2, sum) / apply(A.mat*IPW.ate*tilt.wate, 2, sum)
  p0h <- apply((1-A.mat)*IPW.ate*tilt.wate*Y.mat, 2, sum) / apply((1-A.mat)*IPW.ate*tilt.wate, 2, sum)

  p1 <- apply(A.mat*IPW.att*tilt.watt*Y.mat, 2, sum) / apply(A.mat*IPW.att*tilt.watt, 2, sum)
  p0g <- apply((1-A.mat)*IPW.att*tilt.watt*Y.mat, 2, sum) / apply((1-A.mat)*IPW.att*tilt.watt, 2, sum)

  p1g <- apply(A.mat*IPW.atc*tilt.watc*Y.mat, 2, sum) / apply(A.mat*IPW.atc*tilt.watc, 2, sum)
  p0 <- apply((1-A.mat)*IPW.atc*tilt.watc*Y.mat, 2, sum) / apply((1-A.mat)*IPW.atc*tilt.watc, 2, sum)

  RD.wate <- p1h - p0h
  RR.wate <- p1h/p0h
  OR.wate <- p1h*(1-p0h) / (p0h*(1-p1h))

  RD.watt <- p1 - p0g
  RR.watt <- p1/p0g
  OR.watt <- p1*(1-p0g) / (p0g*(1-p1))

  RD.watc <- p1g - p0
  RR.watc <- p1g/p0
  OR.watc <- p1g*(1-p0) / (p0*(1-p1g))

  return(list(RD.wate=RD.wate, RR.wate=RR.wate, OR.wate=OR.wate,
              RD.watt=RD.watt[-c(2,3)], RR.watt=RR.watt[-c(2,3)], OR.watt=OR.watt[-c(2,3)],
              RD.watc=RD.watc[-c(2,3)], RR.watc=RR.watc[-c(2,3)], OR.watc=OR.watc[-c(2,3)]))
}

.propensity <- function(A, Y, X, X.pred, ps.library="SL.glm") {
  X <- as.data.frame(X)
  fit <- SuperLearner(Y=A, X=X, family=binomial(), SL.library=ps.library)
  e.h <- predict(fit, X.pred, type="response")$pred
  return(e.h)
}
