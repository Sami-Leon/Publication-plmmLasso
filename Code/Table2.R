source("CreateSim.R")
source("modelselection.R")
source("posi.R")
source("test.nonlinear.functions.R")
source("lasso_inference.R")


lambda.grid <- round(exp(seq(log(0.1), log(1 * 0.0001),
  length.out = 10
)), digits = 4)[5:6]

gamma.grid <- c(
  0.000001, 0.0000001, 0.00000001, 0.000000001
)[2:3]

compute_posi <- function(N = 100, n.mvnorm = 100, timepoints = 3:5,
                         sample_from = seq(0, 52, 13), seed = 1,
                         gamma, lambda) {
  sim_data <- simulate_group_inter(
    N = N,
    n.mvnorm = n.mvnorm,
    grouped = TRUE,
    timepoints = timepoints,
    nonpara.inter = TRUE,
    sample_from = sample_from,
    cst_ni = FALSE,
    seed = seed
  )

  fit_plmm <- select.plmm(
    data = sim_data[[1]],
    gamma = gamma,
    lambda = lambda,
    intercept = TRUE,
    timexgroup = TRUE
  )

  return(debias.plmm(sim_data[[1]], fit_plmm))
}


n.sim <- 2
n100.p100 <- function(seed) {
  compute_posi(
    N = 100,
    n.mvnorm = 100,
    timepoints = 3:5,
    sample_from = seq(0, 52, 13),
    seed = seed,
    gamma = gamma.grid,
    lambda = lambda.grid
  )
}

list.n100.p100 <- lapply(1:n.sim, n100.p100)

fixed.lambda <- cov.list(lapply(list.n100.p100, function(x) {
  x$SI
}))

de.sparsified <- cov.list(lapply(list.n100.p100, function(x) {
  x$de.sparsified
}))

naive.debias <- cov.list(lapply(list.n100.p100, function(x) {
  x$debias
}))

adaptive.debias <- cov.list(lapply(list.n100.p100, function(x) {
  x$adapt.debias
}))

fixed.lambda
de.sparsified
naive.debias
adaptive.debias


n100.p500 <- function(seed) {
  compute_posi(
    N = 100,
    n.mvnorm = 500,
    timepoints = 3:5,
    sample_from = seq(0, 52, 13),
    seed = seed,
    gamma = gamma.grid,
    lambda = lambda.grid
  )
}

list.n100.p500 <- lapply(1:n.sim, n100.p500)

fixed.lambda <- cov.list(lapply(list.n100.p500, function(x) {
  x$SI
}))

de.sparsified <- cov.list(lapply(list.n100.p500, function(x) {
  x$de.sparsified
}))

naive.debias <- cov.list(lapply(list.n100.p500, function(x) {
  x$debias
}))

adaptive.debias <- cov.list(lapply(list.n100.p500, function(x) {
  x$adapt.debias
}))

fixed.lambda
de.sparsified
naive.debias
adaptive.debias


n500.p100 <- function(seed) {
  compute_posi(
    N = 500,
    n.mvnorm = 100,
    timepoints = 3:5,
    sample_from = seq(0, 52, 13),
    seed = seed,
    gamma = gamma.grid,
    lambda = lambda.grid
  )
}

list.n500.p100 <- lapply(1:n.sim, n500.p100)

fixed.lambda <- cov.list(lapply(list.n500.p100, function(x) {
  x$SI
}))

de.sparsified <- cov.list(lapply(list.n500.p100, function(x) {
  x$de.sparsified
}))

naive.debias <- cov.list(lapply(list.n500.p100, function(x) {
  x$debias
}))

adaptive.debias <- cov.list(lapply(list.n500.p100, function(x) {
  x$adapt.debias
}))

fixed.lambda
de.sparsified
naive.debias
adaptive.debias


n100.p100.d <- function(seed) {
  compute_posi(
    N = 100,
    n.mvnorm = 100,
    timepoints = 10:14,
    sample_from = seq(0, 52, 4),
    seed = seed,
    gamma = gamma.grid,
    lambda = lambda.grid
  )
}

list.n100.p100.d <- lapply(1:n.sim, n100.p100.d)

fixed.lambda <- cov.list(lapply(list.n100.p100.d, function(x) {
  x$SI
}))

de.sparsified <- cov.list(lapply(list.n100.p100.d, function(x) {
  x$de.sparsified
}))

naive.debias <- cov.list(lapply(list.n100.p100.d, function(x) {
  x$debias
}))

adaptive.debias <- cov.list(lapply(list.n100.p100.d, function(x) {
  x$adapt.debias
}))

fixed.lambda
de.sparsified
naive.debias
adaptive.debias
