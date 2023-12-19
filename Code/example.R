

data.sim = simulate_group_inter(N = 100, n.mvnorm = 100, grouped = T,
                                seed = 12, timepoints = 3:5, nonpara.inter = T,
                                sample_from = seq(0,52,13), cst_ni = F)

lambda.grid <- round(exp(seq(log(0.1), log(1 * 0.0001),
                             length.out = 10
)), digits = 4)

gamma.grid <- c(
  0.000001, 0.0000001, 0.00000001, 0.000000001
)

fit.plmm <- select.plmm(
  data = Group1, gamma = gamma.grid, lambda = lambda.grid,
  crit = "BIC", intercept = T, timexgroup = T
)

posi = de.sparsified.PLMM(simu = Group1, model = test.c)

test.nonlinear.functions = f.test(Group1, test.c, n = 1000)

