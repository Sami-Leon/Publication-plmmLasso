source("CreateSim.R")
source("modelselection.R")

lambda.grid <- round(exp(seq(log(0.1), log(1 * 0.0001),
  length.out = 10
)), digits = 4)[5:6]

gamma.grid <- c(
  0.000001, 0.0000001, 0.00000001, 0.000000001
)[2:3]


sim.a <- simulate_group_inter(
  N = 100, n.mvnorm = 100, grouped = T,
  timepoints = 3:5, nonpara.inter = T,
  sample_from = seq(0, 52, 13), cst_ni = F
)

fit.plmm.a <- select.plmm(
  data = sim.a[[1]], gamma = gamma.grid, lambda = lambda.grid,
  intercept = T, timexgroup = T
)

plot.fit.a <- plot.fit(fit.plmm.a, sim.a)

sim.b <- simulate_group_inter(
  N = 100, n.mvnorm = 500, grouped = T,
  timepoints = 3:5, nonpara.inter = T,
  sample_from = seq(0, 52, 13), cst_ni = F
)

fit.plmm.b <- select.plmm(
  data = sim.b[[1]], gamma = gamma.grid, lambda = lambda.grid,
  intercept = T, timexgroup = T
)

plot.fit.b <- plot.fit(fit.plmm.b, sim.b)

sim.c <- simulate_group_inter(
  N = 500, n.mvnorm = 100, grouped = T,
  timepoints = 3:5, nonpara.inter = T,
  sample_from = seq(0, 52, 13), cst_ni = F
)

fit.plmm.c <- select.plmm(
  data = sim.c[[1]], gamma = gamma.grid, lambda = lambda.grid,
  intercept = T, timexgroup = T
)

plot.fit.c <- plot.fit(fit.plmm.c, sim.c)

sim.d <- simulate_group_inter(
  N = 100, n.mvnorm = 100, grouped = T,
  timepoints = 10:14, nonpara.inter = T,
  sample_from = seq(0, 52, 4), cst_ni = F
)

fit.plmm.d <- select.plmm(
  data = sim.d[[1]], gamma = gamma.grid, lambda = lambda.grid,
  intercept = T, timexgroup = T
)

plot.fit.d <- plot.fit(fit.plmm.d, sim.d)

ggarrange(plot.fit.a, plot.fit.b, plot.fit.c, plot.fit.d,
  ncol = 2, nrow = 2,
  labels = c("(a)", "(b)", "(c)", "(d)"),
  font.label = list(color = "black", size = 16)
)
