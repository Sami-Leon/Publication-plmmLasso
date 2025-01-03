source("CreateSim.R")
source("modelselection.R")
source("posi.R")
source("test.nonlinear.functions.R")

metric.table3 <- function(fit.plmm) {
  theta <- fit.plmm$Res.F$theta
  out <- c(
    round(mean((theta[-1] - c(
      c(3, 2, 1, 0, 0.8, 0, 1, rep(0, 94)),
      c(-1.5, 0, -1, -1.2, rep(0, 96))
    ))^2), 2),
    round(100 * (sum(names(theta) %in% c(
      "Group", "Group1_X1", "Group1_X2", "Group1_X4",
      "Group1_X6", "Group0_X1", "Group0_X3",
      "Group0_X4"
    )) / 8), 2),
    round(100 * (sum(theta[-1] != 0) / length(theta[-1]))),
    round(length(fit.plmm$Res.F$FunctSelect) / 2),
    round(fit.plmm$time, 1),
    round(theta[names(theta) %in% c(
      "Group", "Group1_X1", "Group1_X2", "Group1_X4",
      "Group1_X6", "Group0_X1", "Group0_X3",
      "Group0_X4"
    )], 2),
    round(sqrt(fit.plmm$su), 2)
  )
  return(out)
}




compute_metric <- function(N, n.mvnorm, timepoints,
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
    timexgroup = TRUE, crit = "BIC.nonpara"
  )

  return(metric.table3(fit_plmm))
}


lambda.grid <- round(exp(seq(log(0.1), log(1 * 0.0001),
  length.out = 10
)), digits = 4)[c(6:7)]

gamma.grid <- c(
  0.000001, 0.0000001, 0.00000001, 0.000000001
)[2]

n.sim <- 2

cor.0 <- function(seed) {
  compute_metric(
    N = 100,
    n.mvnorm = 100,
    timepoints = 3:15,
    sample_from = seq(0, 52, 13),
    seed = seed,
    gamma = gamma.grid,
    lambda = lambda.grid, cor = 0.01
  )
}

list.cor.0 <- lapply(1:n.sim, cor.0)


cor.0.5 <- function(seed) {
  compute_metric(
    N = 100,
    n.mvnorm = 100,
    timepoints = 3:15,
    sample_from = seq(0, 52, 13),
    seed = seed,
    gamma = gamma.grid,
    lambda = lambda.grid, cor = 0.5
  )
}

list.cor.0.5 <- lapply(1:n.sim, cor.0.5)


cor.0.9 <- function(seed) {
  compute_metric(
    N = 100,
    n.mvnorm = 100,
    timepoints = 3:15,
    sample_from = seq(0, 52, 13),
    seed = seed,
    gamma = gamma.grid,
    lambda = lambda.grid, cor = 0.9
  )
}

list.cor.0.9 <- lapply(1:n.sim, cor.0.9)


table3 <- cbind(
  rowMeans(do.call("cbind", list.cor.0)),
  rowMeans(do.call("cbind", list.cor.0.5)),
  rowMeans(do.call("cbind", list.cor.0.9))
)

colnames(table3) <- c(
  "n = 100, p = 100, rho = 0",
  "n = 100, p = 100, rho = 0.5",
  "n = 500, p = 100, rho = 0.9"
)
rownames(table3) <- c(
  "MSE", "beta_true (%)", "beta_nonzero (%)", "psi_j(t)", "Time (sec)",
  "Group (3)", "beta_1,2 (2)", "beta_1,3 (1)", "beta_1,6 (0.8)",
  "beta_1,8 (1)", "beta_0,2 (-1.5)", "beta_0,4 (-1)",
  "beta_0,5 (-1.2)", "sigma_phi (0.7)"
)

table3.sd <- round(cbind(
  apply(do.call("cbind", list.cor.0), 1, sd)[6:14],
  apply(do.call("cbind", list.cor.0.5), 1, sd)[6:14],
  apply(do.call("cbind", list.cor.0.9), 1, sd)[6:14]
), 2)

colnames(table3.sd) <- colnames(table3)
rownames(table3.sd) <- c(
  "Group (3)", "beta_1,2 (2)", "beta_1,3 (1)", "beta_1,6 (0.8)",
  "beta_1,8 (1)", "beta_0,2 (-1.5)", "beta_0,4 (-1)",
  "beta_0,5 (-1.2)", "sigma_phi (0.7)"
)

table3

table3.sd
