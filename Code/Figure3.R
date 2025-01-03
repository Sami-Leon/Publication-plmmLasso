source("CreateSim.R")
source("modelselection.R")
source("test.nonlinear.functions.R")

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

test.nonlinear.functions.a <- f.test(sim.a[[1]], fit.plmm.a, n = 100)

f.CI.dis.a <- test.nonlinear.functions.a[[2]][test.nonlinear.functions.a[[2]]$Month %in% sort(unique(sim.a[[1]]$position)), ]

plot.CI.a <- test.nonlinear.functions.a[[3]] + geom_line(aes(Month, `Group diff.`),
  data = f.CI.dis.a, col = "red", linetype = "longdash"
) +
  geom_line(aes(Month, `CI lower`), data = f.CI.dis.a, col = "blue", linetype = "longdash") +
  geom_line(aes(Month, `CI upper`), data = f.CI.dis.a, col = "blue", linetype = "longdash") +
  theme(text = element_text(size = 16)) +
  geom_line(aes(Month, `Group diff.`),
    data = test.nonlinear.functions.a[[2]], col = "red"
  ) +
  ylim(c(-1.1, 1.1))


sim.b <- simulate_group_inter(
  N = 100, n.mvnorm = 500, grouped = T,
  timepoints = 3:5, nonpara.inter = T,
  sample_from = seq(0, 52, 13), cst_ni = F
)

fit.plmm.b <- select.plmm(
  data = sim.b[[1]], gamma = gamma.grid, lambda = lambda.grid,
  intercept = T, timexgroup = T
)

test.nonlinear.functions.b <- f.test(sim.b[[1]], fit.plmm.b, n = 100)

f.CI.dis.b <- test.nonlinear.functions.b[[2]][test.nonlinear.functions.b[[2]]$Month %in% sort(unique(sim.b[[1]]$position)), ]

plot.CI.b <- test.nonlinear.functions.b[[3]] + geom_line(aes(Month, `Group diff.`),
  data = f.CI.dis.b, col = "red", linetype = "longdash"
) +
  geom_line(aes(Month, `CI lower`), data = f.CI.dis.b, col = "blue", linetype = "longdash") +
  geom_line(aes(Month, `CI upper`), data = f.CI.dis.b, col = "blue", linetype = "longdash") +
  theme(text = element_text(size = 16)) +
  geom_line(aes(Month, `Group diff.`),
    data = test.nonlinear.functions.b[[2]], col = "red"
  ) +
  ylim(c(-1.1, 1.1))

sim.c <- simulate_group_inter(
  N = 500, n.mvnorm = 100, grouped = T,
  timepoints = 3:5, nonpara.inter = T,
  sample_from = seq(0, 52, 13), cst_ni = F
)

fit.plmm.c <- select.plmm(
  data = sim.c[[1]], gamma = gamma.grid, lambda = lambda.grid,
  intercept = T, timexgroup = T
)

test.nonlinear.functions.c <- f.test(sim.c[[1]], fit.plmm.c, n = 100)

f.CI.dis.c <- test.nonlinear.functions.c[[2]][test.nonlinear.functions.c[[2]]$Month %in% sort(unique(sim.c[[1]]$position)), ]

plot.CI.c <- test.nonlinear.functions.c[[3]] + geom_line(aes(Month, `Group diff.`),
  data = f.CI.dis.c, col = "red", linetype = "longdash"
) +
  geom_line(aes(Month, `CI lower`), data = f.CI.dis.c, col = "blue", linetype = "longdash") +
  geom_line(aes(Month, `CI upper`), data = f.CI.dis.c, col = "blue", linetype = "longdash") +
  theme(text = element_text(size = 16)) +
  geom_line(aes(Month, `Group diff.`),
    data = test.nonlinear.functions.c[[2]], col = "red"
  ) +
  ylim(c(-1.1, 1.1))


sim.d <- simulate_group_inter(
  N = 100, n.mvnorm = 100, grouped = T,
  timepoints = 10:14, nonpara.inter = T,
  sample_from = seq(0, 52, 4), cst_ni = F
)

fit.plmm.d <- select.plmm(
  data = sim.d[[1]], gamma = gamma.grid, lambda = lambda.grid,
  intercept = T, timexgroup = T
)

test.nonlinear.functions.d <- f.test(sim.d[[1]], fit.plmm.d, n = 100)

f.CI.dis.d <- test.nonlinear.functions.d[[2]][test.nonlinear.functions.d[[2]]$Month %in% sort(unique(sim.d[[1]]$position)), ]

plot.CI.d <- test.nonlinear.functions.d[[3]] + geom_line(aes(Month, `Group diff.`),
  data = f.CI.dis.d, col = "red", linetype = "longdash"
) +
  geom_line(aes(Month, `CI lower`), data = f.CI.dis.d, col = "blue", linetype = "longdash") +
  geom_line(aes(Month, `CI upper`), data = f.CI.dis.d, col = "blue", linetype = "longdash") +
  theme(text = element_text(size = 16)) +
  geom_line(aes(Month, `Group diff.`),
    data = test.nonlinear.functions.d[[2]], col = "red"
  ) +
  ylim(c(-1.1, 1.1))

sim.e <- simulate_group_inter(
  N = 100, n.mvnorm = 100, grouped = T,
  timepoints = 3:5, nonpara.inter = F,
  sample_from = seq(0, 52, 13), cst_ni = F
)

fit.plmm.e <- select.plmm(
  data = sim.e[[1]], gamma = gamma.grid, lambda = lambda.grid,
  intercept = T, timexgroup = T
)

test.nonlinear.functions.e <- f.test(sim.e[[1]], fit.plmm.e, n = 100)

f.CI.dis.e <- test.nonlinear.functions.e[[2]][test.nonlinear.functions.e[[2]]$Month %in% sort(unique(sim.e[[1]]$position)), ]

plot.CI.e <- test.nonlinear.functions.e[[3]] + geom_line(aes(Month, `Group diff.`),
  data = f.CI.dis.e, col = "red", linetype = "longdash"
) +
  geom_line(aes(Month, `CI lower`), data = f.CI.dis.e, col = "blue", linetype = "longdash") +
  geom_line(aes(Month, `CI upper`), data = f.CI.dis.e, col = "blue", linetype = "longdash") +
  theme(text = element_text(size = 16)) +
  geom_line(aes(Month, `Group diff.`),
    data = test.nonlinear.functions.e[[2]], col = "red"
  ) +
  ylim(c(-1.1, 1.1))


first_row <- ggarrange(plot.CI.a, plot.CI.b, ncol = 2, labels = c("(a)", "(b)"))
second_row <- ggarrange(plot.CI.c, plot.CI.d, ncol = 2, labels = c("(c)", "(d)"))
last_row <- ggarrange(NULL, plot.CI.e, NULL,
  ncol = 3, labels = c("", "(e)", ""),
  widths = c(1, 2.1, 1), legend = "bottom"
)

ggarrange(first_row, second_row, last_row, ncol = 1, heights = c(1, 1, 1.3))
