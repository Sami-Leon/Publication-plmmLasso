source("CreateSim.R")
source("modelselection.R")

data.sim.a <- simulate_group_inter(
  N = 100, n.mvnorm = 100, grouped = F,
  seed = 12, timepoints = 3:5, nonpara.inter = F,
  sample_from = seq(0, 52, 13), cst_ni = F
)

lambda.grid <- round(exp(seq(log(0.1), log(1 * 0.0001),
  length.out = 10
)), digits = 4)[5:6]

gamma.grid <- c(
  0.000001, 0.0000001, 0.00000001, 0.000000001
)[2:3]

fit.plmm1.a <- select.plmm(
  data = data.sim.a[[1]], gamma = gamma.grid, lambda = lambda.grid,
  crit = "BIC", intercept = T, timexgroup = T
)

plot.a <- plot.fit(fit.plmm1.a, data.sim.a, same = T, grouped = F, plot.f = FALSE)

data.sim.b <- simulate_group_inter(
  N = 100, n.mvnorm = 100, grouped = T,
  seed = 12, timepoints = 3:5, nonpara.inter = F,
  sample_from = seq(0, 52, 13), cst_ni = F
)

fit.plmm.b <- select.plmm(
  data = data.sim.b[[1]], gamma = gamma.grid, lambda = lambda.grid,
  crit = "BIC", intercept = T, timexgroup = T
)

plot.b <- plot.fit(fit.plmm.b, data.sim.b, same = T, grouped = T, plot.f = FALSE)

data.sim.c <- simulate_group_inter(
  N = 100, n.mvnorm = 100, grouped = F,
  seed = 12, timepoints = 3:5, nonpara.inter = T,
  sample_from = seq(0, 52, 13), cst_ni = F
)

fit.plmm.c <- select.plmm(
  data = data.sim.c[[1]], gamma = gamma.grid, lambda = lambda.grid,
  crit = "BIC", intercept = T, timexgroup = T
)

plot.c <- plot.fit(fit.plmm.c, data.sim.c, same = F, grouped = F, plot.f = FALSE)

data.sim.d <- simulate_group_inter(
  N = 100, n.mvnorm = 100, grouped = T,
  seed = 12, timepoints = 3:5, nonpara.inter = T,
  sample_from = seq(0, 52, 13), cst_ni = F
)

fit.plmm.d <- select.plmm(
  data = data.sim.d[[1]], gamma = gamma.grid, lambda = lambda.grid,
  crit = "BIC", intercept = T, timexgroup = T
)

plot.d <- plot.fit(fit.plmm.d, data.sim.d, same = F, grouped = T, plot.f = FALSE)


ggarrange(plot.a, plot.b, plot.c, plot.d,
  ncol = 2, nrow = 2, common.legend = TRUE,
  legend = "bottom", labels = c("(a)", "(b)", "(c)", "(d)")
)
