library(dplyr)
library(glmnet)
library(hdi)
library(selectiveInference)
source("lasso_inference.R")
roundpval <- function(pvalue) {
  pvalue <- as.numeric(as.character(pvalue))
  pvareturn <- rep(NA, length(pvalue))

  if (!(all(pvalue <= 1 & pvalue >= 0))) {
    warning("At least one p-value is not in the 0, 1 interval")
  }

  for (i in 1:length(pvalue)) {
    pva <- pvalue[i]
    if (is.na(pva)) {
      pvareturn[i] <- NA
    } else {
      if (pva >= 0.06) {
        pvareturn[i] <- as.character(format(round(pva, digits = 2),
          nsmall = 2
        ))
      } else {
        if (pva < 0.06 & pva >= 0.04) {
          pvareturn[i] <- as.character(format(round(pva, digits = 3),
            nsmall = 3
          ))
        } else {
          if (pva < 0.04 & pva >= 0.01) {
            pvareturn[i] <- as.character(format(round(pva, digits = 2),
              nsmall = 2
            ))
          } else {
            if (pva < 0.01 & pva >= 0.001) {
              pvareturn[i] <-
                as.character(format(round(pva, digits = 3),
                  nsmall = 3
                ))
            } else {
              if (pva < 0.001 & pva >= 0.0001) {
                pvareturn[i] <-
                  as.character(format(round(pva, digits = 4),
                    nsmall = 4,
                    scientific = FALSE
                  ))
              } else {
                if (pva < 0.0001) {
                  pvareturn[i] <- "$<0.0001$"
                }
              }
            }
          }
        }
      }
    }
  }
  return(pvareturn)
}

debias.plmm <- function(simu, model, a = 1, level = 0.95, Z = NULL) {
  posi <- NULL
  data.lmm <- simu
  data.lmm$Y <- data.lmm$Y - model$Res.F$out.F$F.fit
  se <- model$se
  lambda <- model$hyper.parameters$lambda.grid

  X <- as.matrix(subset(data.lmm, select = -c(series, position, Y)))
  y <- data.lmm$Y

  y.debias <- scale(y - rep(model$U2$U2, as.vector(table(simu$series))))
  x.debias <- scale(X)

  x.debias.sd <- attr(x.debias, "scaled:scale")
  y.debias.sd <- attr(y.debias, "scaled:scale")

  gfit <- glmnet(x.debias,
    y.debias,
    alpha = 1, lambda = lambda,
    standardize = FALSE,
    intercept = TRUE, thresh = 1e-20
  )

  lasso.test <- as.vector(coef(gfit))

  SSLasso.out <- SSLasso(
    X = x.debias, y = y.debias,
    theta = lasso.test, alpha = 0.05, lambda = lambda,
    mu = NULL, resol = 1.3, maxiter = 50,
    threshold = 1e-2, verbose = F
  )


  SSLasso.out$low.lim <- SSLasso.out$low.lim * (y.debias.sd / x.debias.sd)
  SSLasso.out$up.lim <- SSLasso.out$up.lim * (y.debias.sd / x.debias.sd)


  posi$debias <- list(
    ci = cbind(SSLasso.out$low.lim, SSLasso.out$up.lim),
    pv = SSLasso.out$pvals
  )



  beta <- coef(gfit,
    x = x.debias, y = y.debias, s = lambda / length(y.debias),
    exact = TRUE
  )[-1]

  out <- fixedLassoInf(x.debias, y.debias, beta, lambda,
    sigma = sqrt(se),
    alpha = 0.05
  )

  nonzero <- names(x.debias.sd) %in% names(out$vars)

  out$ci[, 1] <- out$ci[, 1] * (y.debias.sd / x.debias.sd[nonzero])
  out$ci[, 2] <- out$ci[, 2] * (y.debias.sd / x.debias.sd[nonzero])

  ci.low <- rep(0, length(x.debias.sd))
  ci.low[nonzero] <- out$ci[, 1]

  ci.up <- rep(0, length(x.debias.sd))
  ci.up[nonzero] <- out$ci[, 2]

  SI.ci <- cbind(ci.low, ci.up)

  SI.pv <- rep(1, length(x.debias.sd))
  SI.pv[nonzero] <- out$pv
  posi$SI <- list(ci = SI.ci, pv = SI.pv)

  de.sparsified <- lasso.proj(x.debias, y.debias,
    sigma = sqrt(se),
    suppress.grouptesting = TRUE, do.ZnZ = T,
    betainit = lasso.test[-1]
  )
  de.sparsified.ci <- confint(de.sparsified, level = 0.95)
  de.sparsified.ci <- cbind(
    de.sparsified.ci[, 1] * (y.debias.sd / x.debias.sd),
    de.sparsified.ci[, 2] * (y.debias.sd / x.debias.sd)
  )

  posi$de.sparsified <- list(
    ci = de.sparsified.ci,
    pv = de.sparsified$pval
  )

  grp <- factor(data.lmm$series)

  z <- model.matrix(~ as.factor(series) - 1, data.lmm[, "series", drop = F])

  Sigma <- a * model$su * z %*% t(z) + model$se * diag(rep(1, length(y)))

  Sigma.svd <- svd(Sigma)
  Sig.a.inv.half <- Sigma.svd$u %*% diag(1 / sqrt(Sigma.svd$d)) %*% t(Sigma.svd$u)

  X.a <- Sig.a.inv.half %*% X
  y.a <- Sig.a.inv.half %*% y

  if (is.null(Z)) {
    de.sparsified.sl <- hdi::lasso.proj(X.a, y.a,
      suppress.grouptesting = TRUE, return.Z = T,
      do.ZnZ = T, betainit = "scaled lasso"
    )

    lp.Z <- de.sparsified.sl$Z
  } else {
    lp.Z <- Z
  }


  N <- length(y)
  n <- length(unique(grp))
  q <- ncol(z)
  p <- ncol(X)

  beta.db.mlm <- NULL
  beta.db.sd.mlm <- NULL

  beta.hat <- model$Res.F$theta[names(model$Res.F$theta) %in% colnames(X.a)]

  res <- y.a - X.a %*% beta.hat

  for (j in 1:length(beta.hat)) {
    col.j <- j

    wj.mlm <- lp.Z[, j]

    beta.db.mlm[j] <- beta.hat[col.j] + sum(wj.mlm * res) / sum(wj.mlm * X.a[, col.j])

    ok <- data.frame(wj.mlm, res, grp) %>%
      group_by(grp) %>%
      summarize(num = (sum(wj.mlm * res))^2)
    num <- sum(ok$num)

    beta.db.sd.mlm[j] <- sqrt(num) / sqrt((sum(wj.mlm * X.a[, col.j]))^2)
  }

  ci.ad <- cbind(
    beta.db.mlm - 1.96 * beta.db.sd.mlm,
    beta.db.mlm + 1.96 * beta.db.sd.mlm
  )
  bprojrescaled <- beta.db.mlm * (1 / beta.db.sd.mlm)
  pv.ad <- 2 * pnorm(abs(bprojrescaled), lower.tail = FALSE)

  posi$adapt.debias <- list(
    ci = ci.ad, pv = pv.ad, beta.hat = beta.hat,
    bhat = beta.db.mlm, bhat.sd = beta.db.sd.mlm
  )

  return(posi)
}


compute.cov <- function(x) {
  upper <- c(3, 2, 1, rep(0, 2)) < x[1:5, 2]
  lower <- x[1:5, 1] < c(3, 2, 1, rep(0, 2))
  return(upper & lower)
}


compute.length <- function(x) {
  return(x[1:5, 2] - x[1:5, 1])
}


cov.list <- function(list.ci) {
  coverage <- colMeans(do.call("rbind", lapply(lapply(list.ci, function(x) {
    x$ci
  }), compute.cov)))

  ci.length <- do.call("rbind", lapply(lapply(list.ci, function(x) {
    x$ci
  }), compute.length))

  length <- apply(ci.length, 2, function(x) mean(x[abs(x) != Inf]))

  pvalue <- apply(do.call("rbind", lapply(list.ci, function(x) {
    x$pv[1:5]
  })), 2, function(x) mean(na.omit(x)))

  out <- as.data.frame(cbind(coverage, length, pvalue))
  out[, 1] <- round(out[, 1], 3)
  out[, 2] <- round(out[, 2], 2)
  out[, 3] <- roundpval(out[, 3])

  return(out)
}


compute.cov.irr <- function(x) {
  upper <- c(3, 2, 1, 0, 0.8, 0, 1, -1.5, 0, -1, -1.2) < x[c(1:7, 102:105), 2]
  lower <- x[c(1:7, 102:105), 1] < c(3, 2, 1, 0, 0.8, 0, 1, -1.5, 0, -1, -1.2)

  return(upper & lower)
}

compute.length.irr <- function(x) {
  return(x[c(1:7, 102:105), 2] - x[c(1:7, 102:105), 1])
}

cov.list.irr <- function(list.ci) {
  coverage <- colMeans(do.call("rbind", lapply(lapply(list.ci, function(x) {
    x$ci
  }), compute.cov.irr)))

  ci.length <- do.call("rbind", lapply(lapply(list.ci, function(x) {
    x$ci
  }), compute.length.irr))

  length <- apply(ci.length, 2, function(x) mean(x[abs(x) != Inf]))

  pvalue <- apply(do.call("rbind", lapply(list.ci, function(x) {
    x$pv[c(1:7, 102:105)]
  })), 2, function(x) mean(na.omit(x)))


  out <- as.data.frame(cbind(coverage, length, pvalue))
  out[, 1] <- round(out[, 1], 3)
  out[, 2] <- round(out[, 2], 2)
  out[, 3] <- roundpval(out[, 3])

  return(out)
}
