#' Causal inference via general estimands: the weighted ATE, weighted ATT and weighted ATC
#' @param A vector of the treatment variable (binary, valued from 0 and 1)
#' @param Y vector of the outcome variable (can be various types: continous, categorical, etc., but must be numerical)
#' @param X matrix of covariates/confounders that are included into the propensity score model
#' @param ps.library methods used for fitting the propensity score model, adopted from the SuperLearner package
#' @param beta whether to include estimands from beta family weights; if so, the v parameter below is required to be specified
#' @param v model parameter for beta weights
#' @param trim whether to include trimming; if so, the trim.alpha parameter below needs to be specified
#' @param trim.alpha the trimming threshold; default is .05
#' @param trun whether to include truncation; if so, the trun.alpha parameter below needs to be specified
#' @param trun.alpha the truncation threshold; default is .05
#' @param n.folds number of data folds by each sample-splitting replication; default is 5
#' @param n.split number of sample-splitting replications; default is 10
#' @param seed seed for initializing the random number generator, which is used in set.seed() function; default is 4399
#' @param conf.level level of confidence interval; default is .05 (for 0.95 confidence interval, using normal approximation)
#' @return the function returns estimated causal estimands;
#'         estimands defined by constant weights, overlap weights, matching weights and entropy weights will always be returned;
#'         whether to return estimands defined by trimming, truncation, and beta weights can be specified by user
Chips_DML <- function(A, Y, X, ps.library="glm",
                      beta=FALSE, v=NA,
                      trim=FALSE, trim.alpha=.05,
                      trun=FALSE, trun.alpha=.05,
                      boot=TRUE, n.folds=5, n.split=10,
                      seed=4399, conf.level=.05) {

  n <- length(A)
  set.seed(seed=seed)
  seeds <- round(runif(n.split, min=0, max=10^3*n.split))
  Est.wate <- SE.wate <- Est.watt <- SE.watt <- Est.watc <- SE.watc <- NULL
  for(k in 1:n.split) {
    set.seed(seeds[k])
    pred.folds <- createFolds(1:n, k=n.folds, list=T)
    wate.est <- wate.se <- watt.est <- watt.se <- watc.est <- watc.se <- c()
    for(i in 1:n.folds) {
      pred.ind <- pred.folds[[i]]
      train.ind <- (1:n)[-pred.ind]

      e.h <- .propensity(A=A[train.ind],
                         Y=Y[train.ind],
                         X=X[train.ind,],
                         X.pred=X[pred.ind,],
                         ps.library=ps.library)

      result <- .point.est.DML(A=A[pred.ind], Y=Y[pred.ind],
                               e.h=e.h, beta=beta, v=v,
                               trim=trim, trun=trun,
                               trim.alpha=trim.alpha,
                               trun.alpha=trun.alpha)

      wate.est <- cbind(wate.est, result$WATE)
      wate.se <- cbind(wate.se, result$WATE.sd)
      watt.est <- cbind(watt.est, result$WATT)
      watt.se <- cbind(watt.se, result$WATT.sd)
      watc.est <- cbind(watc.est, result$WATC)
      watc.se <- cbind(watc.se, result$WATC.sd)
    }
    Est.wate <- rbind(Est.wate, apply(wate.est, 1, mean))
    SE.wate <- rbind(SE.wate, apply(wate.se, 1, mean)/sqrt(n.folds))
    Est.watt <- rbind(Est.watt, apply(watt.est, 1, mean))
    SE.watt <- rbind(SE.watt, apply(watt.se, 1, mean)/sqrt(n.folds))
    Est.watc <- rbind(Est.watc, apply(watc.est, 1, mean))
    SE.watc <- rbind(SE.watc, apply(watc.se, 1, mean)/sqrt(n.folds))

    if(k%%50==0) print(paste0("Sample-splitting ", k, " is done."))
  }

  quant <- qnorm(1-conf.level/2)

  wate <- apply(Est.wate, 2, mean, na.rm=T)
  sd.wate <- apply(SE.wate, 2, mean, na.rm=T)

  watt <- apply(Est.watt, 2, mean, na.rm=T)
  sd.watt <- apply(SE.watt, 2, mean, na.rm=T)

  watc <- apply(Est.watc, 2, mean, na.rm=T)
  sd.watc <- apply(SE.watc, 2, mean, na.rm=T)

  upr.wate <- wate + quant*sd.wate
  upr.watt <- watt + quant*sd.watt
  upr.watc <- watc + quant*sd.watc

  lwr.wate <- wate - quant*sd.wate
  lwr.watt <- watt - quant*sd.watt
  lwr.watc <- watc - quant*sd.watc

  return(list(df.WATE=data.frame(Est=wate, Std.Err=sd.wate, Upr=upr.wate, Lwr=lwr.wate),
              df.WATT=data.frame(Est=watt, Std.Err=sd.watt, Upr=upr.watt, Lwr=lwr.watt),
              df.WATC=data.frame(Est=watc, Std.Err=sd.watc, Upr=upr.watc, Lwr=lwr.watc)))
}

.propensity <- function(A, Y, X, X.pred, ps.library="glm") {
  if(ps.library=="glm") {
    X <- as.data.frame(X)
    colnames(X) <- paste0("X", 1:ncol(X))
    fit <- glm(A~.-A, data=data.frame(A, X), family=binomial(link="logit"))
    e.h <- predict(fit, X.pred, type="response")
  } else {
    X <- as.matrix(X)
    colnames(X) <- paste0("X", 1:ncol(X))
    fit <- SuperLearner(Y=A, X=X, family=binomial(), SL.library=ps.library)
    e.h <- predict(fit, X.pred, type="response")$pred
  }
  return(e.h)
}

.point.est.DML <- function(A, Y, X, e.h,
                           beta=FALSE, v=NA,
                           trim=FALSE, trim.alpha=.05,
                           trun=FALSE, trun.alpha=.05) {

  # define the tilting function
  n <- length(A)
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
          (as.numeric(e.h<trun.alpha[i])/trun.alpha[i] + as.numeric(e.h<trun.alpha[i])/(1-trun.alpha[i])) * e.h^A * (1-e.h)^(1-A)
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

  M <- ncol(tilt.wate)
  Y.mat <- matrix(Y, ncol=M, nrow=n)
  A.mat <- matrix(A, ncol=M, nrow=n)

  IPW.ate.1 <- matrix(A/e.h, ncol=M, nrow=n) / t(matrix(apply(tilt.wate, 2, mean), M, n))
  IPW.ate.0 <- matrix((1-A)/(1-e.h), ncol=M, nrow=n) / t(matrix(apply(tilt.wate, 2, mean), M, n))
  WATE.IF.1 <- IPW.ate.1*(Y.mat-mean(A/e.h*Y))*tilt.wate
  WATE.IF.0 <- IPW.ate.0*(Y.mat-mean((1-A)/(1-e.h)*Y))*tilt.wate
  WATE.IF   <- WATE.IF.1 - WATE.IF.0 + (mean(A/e.h*Y) - mean((1-A)/(1-e.h)*Y)) * tilt.wate / t(matrix(apply(tilt.wate, 2, mean), M, n))
  wate.est  <- apply(WATE.IF, 2, mean)
  wate.sd   <- sqrt(apply(mean(A)*WATE.IF.1, 2, var)/sum(A)+apply(mean(1-A)*WATE.IF.0, 2, var)/sum(1-A))

  IPW.att.1 <- matrix(A, ncol=M, nrow=n) / mean(A)
  IPW.att.0 <- matrix((1-A)*e.h/(1-e.h), ncol=M, nrow=n) / t(matrix(apply(tilt.watt*(1-A)*e.h/(1-e.h), 2, mean), M, n))
  WATT.IF.1 <- IPW.att.1*(Y.mat-sum(A*Y)/sum(A))
  WATT.IF.0 <- IPW.att.0*(Y.mat-sum((1-A)*e.h/(1-e.h)*Y)/sum(A))*tilt.watt
  WATT.IF   <- IPW.att.1*Y.mat - IPW.att.0*Y.mat*tilt.watt
  watt.est  <- apply(WATT.IF, 2, mean)
  watt.sd   <- sqrt(apply(mean(A)*WATT.IF.1, 2, var)/sum(A)+apply(mean(1-A)*WATT.IF.0, 2, var)/sum(1-A))

  IPW.atc.1 <- matrix(A*(1-e.h)/e.h, ncol=M, nrow=n) / t(matrix(apply(tilt.watc*A*(1-e.h)/e.h, 2, mean), M, n))
  IPW.atc.0 <- matrix(1-A, ncol=M, nrow=n) / mean(1-A)
  WATC.IF.1 <- IPW.atc.1*(Y.mat-sum((A*(1-e.h)/e.h*Y)/sum(1-A)))*tilt.watc
  WATC.IF.0 <- IPW.atc.0*(Y.mat-sum((1-A)*Y)/sum(1-A))
  WATC.IF   <- IPW.atc.1*Y.mat*tilt.watc - IPW.atc.0*Y.mat
  watc.est  <- apply(WATC.IF, 2, mean)
  watc.sd   <- sqrt(apply(mean(A)*WATC.IF.1, 2, var)/sum(A)+apply(mean(1-A)*WATC.IF.0, 2, var)/sum(1-A))

  return(list(WATE=wate.est, WATT=watt.est, WATC=watc.est, WATE.sd=wate.sd, WATT.sd=watt.sd, WATC.sd=watc.sd))
}
