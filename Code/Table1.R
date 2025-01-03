source("CreateSim.R")
source("modelselection.R")
source("posi.R")
source("test.nonlinear.functions.R")

metric.table1 <- function(fit.plmm) {
  theta <- fit.plmm$Res.F$theta
  out <- c(
    round(mean((theta[-1] - c(3, 2, 1, rep(0, length(theta) - 4)))^2), 2),
    round(100 * sum(theta[c(2, 3, 4)] != 0) / 3),
    round(100 * (sum(theta[-1] != 0) / length(theta[-1]))),
    round(length(fit.plmm$Res.F$FunctSelect) / 2),
    round(fit.plmm$time, 1),
    round(theta[2:4], 2),
    round(sqrt(fit.plmm$su), 2)
  )
  return(out)
}

compute_metric <- function(N = 100, n.mvnorm = 100, timepoints = 3:5,
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

  return(metric.table1(fit_plmm))
}


lambda.grid <- round(exp(seq(log(0.1), log(1 * 0.0001),
  length.out = 10
)), digits = 4)[5:6]

gamma.grid <- c(
  0.000001, 0.0000001, 0.00000001, 0.000000001
)[2:3]

n.sim <- 2

n100.p100 <- function(seed) {
  compute_metric(
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


n100.p500 <- function(seed) {
  compute_metric(
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


n500.p100 <- function(seed) {
  compute_metric(
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


n100.p100.d <- function(seed) {
  compute_metric(
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

table1 <- cbind(
  rowMeans(do.call("cbind", list.n100.p100)),
  rowMeans(do.call("cbind", list.n100.p500)),
  rowMeans(do.call("cbind", list.n500.p100)),
  rowMeans(do.call("cbind", list.n100.p100.d))
)

colnames(table1) <- c(
  "n = 100, p = 100, n_i in [3,5]",
  "n = 100, p = 500, n_i in [3,5]",
  "n = 500, p = 100, n_i in [3,5]",
  "n = 100, p = 100, n_i in [10,14]"
)
rownames(table1) <- c(
  "MSE", "beta_true (%)", "beta_nonzero (%)", "psi_j(t)", "Time (sec)",
  "Group (3)", "beta_2 (2)", "Group x beta_2 (1)", "sigma_phi (0.7)"
)

table1.sd <- round(cbind(
  apply(do.call("cbind", list.n100.p100), 1, sd)[6:9],
  apply(do.call("cbind", list.n100.p500), 1, sd)[6:9],
  apply(do.call("cbind", list.n500.p100), 1, sd)[6:9],
  apply(do.call("cbind", list.n100.p100.d), 1, sd)[6:9]
), 2)

colnames(table1.sd) <- colnames(table1)
rownames(table1.sd) <- c("Group", "beta_2", "Group x beta_2", "sigma_phi")

table1

table1.sd
