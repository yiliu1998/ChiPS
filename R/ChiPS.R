#' Causal inference via general estimands: the weighted ATE, weighted ATT and weighted ATC
#' @param A vector of the treatment variable (binary, valued from 0 and 1)
#' @param Y vector of the outcome variable (can be various types: continous, categorical, etc., but must be numerical)
#' @param X matrix of covariates/confounders that are included into the propensity score model
#' @param beta whether to include estimands from beta family weights; if so, the v parameter below is required to be specified
#' @param v model parameter for beta weights
#' @param trim whether to include trimming; if so, the trim.reEST and trim.alpha parameter below needs to be specified
#' @param trim.reEST whether to re-estimate the PS after trimming, if TRUE, both results by re-estimating the PS or not will be returned
#' @param trim.alpha the trimming threshold; default is .05
#' @param trun whether to include truncation; if so, the trun.alpha parameter below needs to be specified
#' @param trun.alpha the truncation threshold; default is .05
#' @param boot whether to implement bootstrap for variance estimation; default is TRUE
#' @param n.boot number of bootstrap, if boot==TRUE; default is 500
#' @param return.psfig whether to plot and return the estimated propensity score distribution plots by treatment group
#' @param seed seed for initializing the random number generator, which is used in set.seed() function; default is 4399
#' @param conf.level level of confidence interval; default is .05 (for 0.95 confidence interval, using normal approximation)
#' @return the function returns estimated causal estimands;
#'         estimands defined by constant weights, overlap weights, matching weights and entropy weights will always be returned;
#'         whether to return estimands defined by trimming, truncation, and beta weights can be specified by user
ChiPS <- function(A, Y, X,
                  ps.library="SL.glm",
                  beta=FALSE, v=NA,
                  trim=FALSE, trim.reEST=FALSE, trim.alpha=.05,
                  trun=FALSE, trun.alpha=.05,
                  boot=TRUE, n.boot=500,
                  return.psfig=FALSE,
                  seed=4399,
                  conf.level=.05) {

  set.seed(seed=seed)
  n <- length(A)
  X <- as.data.frame(X)
  e.h <- .propensity(A=A, X=X, X.pred=X, ps.library=ps.library)

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
    watc.all <- watt.all <- wate.all <- NULL
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

      wate.all <- rbind(all.result$WATE, wate.all)
      watt.all <- rbind(all.result$WATT, watt.all)
      watc.all <- rbind(all.result$WATC, watc.all)

      if(i%%50==0) print(paste0("Bootstrap ", i, " is done."))
    }
    sd.wate <- apply(wate.all, 2, sd, na.rm=T)
    sd.watt <- apply(watt.all, 2, sd, na.rm=T)
    sd.watc <- apply(watc.all, 2, sd, na.rm=T)

    quant <- qnorm(1-conf.level/2)

    upr.wate <- result$WATE + quant*sd.wate
    upr.watt <- result$WATT + quant*sd.watt
    upr.watc <- result$WATC + quant*sd.watc

    lwr.wate <- result$WATE - quant*sd.wate
    lwr.watt <- result$WATT - quant*sd.watt
    lwr.watc <- result$WATC - quant*sd.watc

    return(list(df.WATE=data.frame(Est=result$WATE, Std.Err=sd.wate, Upr=upr.wate, Lwr=lwr.wate),
                df.WATT=data.frame(Est=result$WATT, Std.Err=sd.watt, Upr=upr.watt, Lwr=lwr.watt),
                df.WATC=data.frame(Est=result$WATC, Std.Err=sd.watc, Upr=upr.watc, Lwr=lwr.watc)))
  }
}

.point.est <- function(A, Y, X,
                       ps.library="SL.glm",
                       beta=FALSE, v=NA,
                       trim=FALSE, trim.reEST=FALSE, trim.alpha=.05,
                       trun=FALSE, trun.alpha=.05) {

  # propensity score estimation
  n <- length(A)
  e.h <- .propensity(A=A, X=X, X.pred=X, ps.library=ps.library)

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
      tilt.trim.wate <- tilt.trim.watt <- tilt.trim.watc <- matrix(NA, nrow=n, ncol=length(trim.alpha))
      for(i in 1:length(trim.alpha)) {
        tilt.trim.wate[,i] <- as.numeric(e.h<1-trim.alpha[i] & e.h>trim.alpha[i])
        tilt.trim.watt[,i] <- as.numeric(e.h<1-trim.alpha[i])*(1-A) + A
        tilt.trim.watc[,i] <- as.numeric(e.h>trim.alpha[i])*A + (1-A)
      }
      tilt.wate <- cbind(tilt.wate, tilt.trim.wate)
      tilt.watt <- cbind(tilt.watt, tilt.trim.watt)
      tilt.watc <- cbind(tilt.watc, tilt.trim.watc)

      if(trim.reEST) {
        weights <- c(weights, paste0("trimming (re-est PS, alpha=", trim.alpha, ")"))
        tilt.trim.reEST.wate <- tilt.trim.reEST.watt <- tilt.trim.reEST.watc <- matrix(NA, nrow=n, ncol=length(trim.alpha))
        for(i in 1:length(trim.alpha)) {
          tilt.trim.reEST.wate[,i] <- as.numeric(e.h<1-trim.alpha[i] & e.h>trim.alpha[i])
          tilt.trim.reEST.watt[,i] <- as.numeric(e.h<1-trim.alpha[i])*(1-A) + A
          tilt.trim.reEST.watc[,i] <- as.numeric(e.h>trim.alpha[i])*A + (1-A)
        }
        tilt.wate <- cbind(tilt.wate, tilt.trim.reEST.wate)
        tilt.watt <- cbind(tilt.watt, tilt.trim.reEST.watt)
        tilt.watc <- cbind(tilt.watc, tilt.trim.reEST.watc)
      }
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

  wate.est <- apply(A.mat*IPW.ate*tilt.wate*Y.mat, 2, sum) / apply(A.mat*IPW.ate*tilt.wate, 2, sum) -
    apply((1-A.mat)*IPW.ate*tilt.wate*Y.mat, 2, sum) / apply((1-A.mat)*IPW.ate*tilt.wate, 2, sum)

  watt.est <- apply(A.mat*IPW.att*tilt.watt*Y.mat, 2, sum) / apply(A.mat*IPW.att*tilt.watt, 2, sum) -
    apply((1-A.mat)*IPW.att*tilt.watt*Y.mat, 2, sum) / apply((1-A.mat)*IPW.att*tilt.watt, 2, sum)

  watc.est <- apply(A.mat*IPW.atc*tilt.watc*Y.mat, 2, sum) / apply(A.mat*IPW.atc*tilt.watc, 2, sum) -
    apply((1-A.mat)*IPW.atc*tilt.watc*Y.mat, 2, sum) / apply((1-A.mat)*IPW.atc*tilt.watc, 2, sum)

  if(trim.reEST) {
    for(i in 1:length(trim.alpha)) {
      A.wate <- A[tilt.trim.reEST.wate[,i]==1]
      A.watt <- A[tilt.trim.reEST.watt[,i]==1]
      A.watc <- A[tilt.trim.reEST.watc[,i]==1]

      Y.wate <- Y[tilt.trim.reEST.wate[,i]==1]
      Y.watt <- Y[tilt.trim.reEST.watt[,i]==1]
      Y.watc <- Y[tilt.trim.reEST.watc[,i]==1]

      X.wate <- X[tilt.trim.reEST.wate[,i]==1,]
      X.watt <- X[tilt.trim.reEST.watt[,i]==1,]
      X.watc <- X[tilt.trim.reEST.watc[,i]==1,]

      e.h.wate <- .propensity(A=A.wate, X=X.wate, X.pred=X.wate, ps.library=ps.library)
      e.h.watt <- .propensity(A=A.watt, X=X.watt, X.pred=X.watt, ps.library=ps.library)
      e.h.watc <- .propensity(A=A.wate, X=X.watc, X.pred=X.watc, ps.library=ps.library)

      k <- length(weights)-length(trim.alpha)+i
      wate.est[k] <- sum(A.wate/e.h.wate * Y.wate)/sum(A.wate/e.h.wate) - sum((1-A.wate)/(1-e.h.wate) * Y.wate)/sum((1-A.wate)/(1-e.h.wate))
      watt.est[k] <- sum(A.watt * Y.watt)/sum(A.watt) - sum((1-A.watt)*e.h.watt/(1-e.h.watt) * Y.watt)/sum((1-A.watt)*e.h.watt/(1-e.h.watt))
      watc.est[k] <- sum(A.watc*(1-e.h.watc)/e.h.watc * Y.watc)/sum(A.watc*(1-e.h.watc)/e.h.watc) - sum((1-A.watc) * Y.watc)/sum(1-A.watc)
    }
  }

  return(list(WATE=wate.est, WATT=watt.est[-c(2,3)], WATC=watc.est[-c(2,3)]))
}

.propensity <- function(A, X, X.pred, ps.library="SL.glm") {
  X <- as.data.frame(X)
  fit <- SuperLearner(Y=A, X=X, family=binomial(), SL.library=ps.library)
  e.h <- predict(fit, X.pred, type="response")$pred
  return(e.h)
}
