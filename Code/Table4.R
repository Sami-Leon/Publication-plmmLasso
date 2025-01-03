source("CreateSim.R")
source("modelselection.R")
source("posi.R")
source("test.nonlinear.functions.R")
source("lasso_inference.R")


lambda.grid <- round(exp(seq(log(0.1), log(1 * 0.0001),
  length.out = 10
)), digits = 4)[c(6:7)]

gamma.grid <- c(
  0.000001, 0.0000001, 0.00000001, 0.000000001
)[2]

n.sim <- 2

compute_posi <- function(N, n.mvnorm, timepoints,
                         sample_from, seed,
                         gamma, lambda, cor) {
  sim_data <- simulate_irregular(
    N = N,
    n.mvnorm = n.mvnorm,
    grouped = TRUE,
    timepoints = timepoints,
    nonpara.inter = TRUE,
    sample_from = sample_from,
    cst_ni = FALSE,
    seed = seed, cor = cor
  )

  X <- subset(sim_data[[1]], select = -c(Y, series, position))
  Groupvar <- X$Group
  X <- X[, -1]
  X.inter <- matrix(nrow = nrow(X), ncol = ncol(X) * 2)
  for (i in 1:nrow(X.inter)) {
    if (Groupvar[i] == 1) {
      X.inter[i, ] <- c(unlist(X[i, ]), rep(0, ncol(X)))
    } else {
      X.inter[i, ] <- c(rep(0, ncol(X)), unlist(X[i, ]))
    }
  }
  X.inter <- as.data.frame(X.inter)
  colnames(X.inter) <- c(
    paste0("Group1_", colnames(X)),
    paste0("Group0_", colnames(X))
  )
  X.inter <- data.frame(Group = Groupvar, X.inter)

  data.inter <- cbind(sim_data[[1]][, 1:3], X.inter)

  fit_plmm <- select.plmm(
    data = data.inter,
    gamma = gamma,
    lambda = lambda,
    intercept = TRUE,
    timexgroup = TRUE
  )
  return(debias.plmm(data.inter, fit_plmm))
}


n100.p100.0 <- function(seed) {
  compute_posi(
    N = 100,
    n.mvnorm = 100,
    timepoints = 3:15,
    sample_from = seq(0, 52, 13),
    seed = seed,
    gamma = gamma.grid,
    lambda = lambda.grid, cor = 0.01
  )
}

list.n100.p100.0 <- lapply(1:n.sim, n100.p100.0)

fixed.lambda <- cov.list.irr(lapply(list.n100.p100.0, function(x) {
  x$SI
}))

de.sparsified <- cov.list.irr(lapply(list.n100.p100.0, function(x) {
  x$de.sparsified
}))

naive.debias <- cov.list.irr(lapply(list.n100.p100.0, function(x) {
  x$debias
}))

adaptive.debias <- cov.list.irr(lapply(list.n100.p100.0, function(x) {
  x$adapt.debias
}))

fixed.lambda
de.sparsified
naive.debias
adaptive.debias


n100.p100.0.5 <- function(seed) {
  compute_posi(
    N = 100,
    n.mvnorm = 100,
    timepoints = 3:15,
    sample_from = seq(0, 52, 13),
    seed = seed,
    gamma = gamma.grid,
    lambda = lambda.grid, cor = 0.5
  )
}

list.n100.p100.0.5 <- lapply(1:n.sim, n100.p100.0.5)

fixed.lambda <- cov.list.irr(lapply(list.n100.p100.0.5, function(x) {
  x$SI
}))

de.sparsified <- cov.list.irr(lapply(list.n100.p100.0.5, function(x) {
  x$de.sparsified
}))

naive.debias <- cov.list.irr(lapply(list.n100.p100.0.5, function(x) {
  x$debias
}))

adaptive.debias <- cov.list.irr(lapply(list.n100.p100.0.5, function(x) {
  x$adapt.debias
}))

fixed.lambda
de.sparsified
naive.debias
adaptive.debias


n100.p100.0.9 <- function(seed) {
  compute_posi(
    N = 100,
    n.mvnorm = 100,
    timepoints = 3:15,
    sample_from = seq(0, 52, 13),
    seed = seed,
    gamma = gamma.grid,
    lambda = lambda.grid, cor = 0.9
  )
}

list.n100.p100.0.9 <- lapply(1:n.sim, n100.p100.0.9)

fixed.lambda <- cov.list.irr(lapply(list.n100.p100.0.9, function(x) {
  x$SI
}))

de.sparsified <- cov.list.irr(lapply(list.n100.p100.0.9, function(x) {
  x$de.sparsified
}))

naive.debias <- cov.list.irr(lapply(list.n100.p100.0.9, function(x) {
  x$debias
}))

adaptive.debias <- cov.list.irr(lapply(list.n100.p100.0.9, function(x) {
  x$adapt.debias
}))

fixed.lambda
de.sparsified
naive.debias
adaptive.debias
