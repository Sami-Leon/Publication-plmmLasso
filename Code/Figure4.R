source("CreateSim.R")
source("modelselection.R")

lambda.grid <- round(exp(seq(log(0.1), log(1 * 0.0001),
  length.out = 10
)), digits = 4)[c(6:7)]

gamma.grid <- c(
  0.000001, 0.0000001, 0.00000001, 0.000000001
)[2]


data.sim.a <- simulate_irregular(
  N = 100, n.mvnorm = 100, grouped = T,
  seed = 1, timepoints = 3:15, nonpara.inter = T,
  sample_from = seq(0, 52, 13), cst_ni = F,
  cor = 0.01
)

X <- subset(data.sim.a[[1]], select = -c(Y, series, position))
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

data.a.inter <- cbind(data.sim.a[[1]][, 1:3], X.inter)

fit.plmm.a <- select.plmm(
  data = data.a.inter, gamma = gamma.grid, lambda = lambda.grid,
  crit = "BIC", intercept = T, timexgroup = T
)


plot.a <- plot.fit.irr(fit.plmm.a, data.sim.a)


data.sim.b <- simulate_irregular(
  N = 100, n.mvnorm = 100, grouped = T,
  seed = 1, timepoints = 3:15, nonpara.inter = T,
  sample_from = seq(0, 52, 13), cst_ni = F,
  cor = 0.5
)

X <- subset(data.sim.b[[1]], select = -c(Y, series, position))
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

data.b.inter <- cbind(data.sim.b[[1]][, 1:3], X.inter)

fit.plmm.b <- select.plmm(
  data = data.b.inter, gamma = gamma.grid, lambda = lambda.grid,
  crit = "BIC", intercept = T, timexgroup = T
)


plot.b <- plot.fit.irr(fit.plmm.b, data.sim.b)


data.sim.c <- simulate_irregular(
  N = 100, n.mvnorm = 100, grouped = T,
  seed = 1, timepoints = 3:15, nonpara.inter = T,
  sample_from = seq(0, 52, 13), cst_ni = F,
  cor = 0.9
)

X <- subset(data.sim.c[[1]], select = -c(Y, series, position))
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

data.c.inter <- cbind(data.sim.c[[1]][, 1:3], X.inter)

fit.plmm.c <- select.plmm(
  data = data.c.inter, gamma = gamma.grid, lambda = lambda.grid,
  crit = "BIC", intercept = T, timexgroup = T
)


plot.c <- plot.fit.irr(fit.plmm.c, data.sim.c)


first_row <- ggarrange(plot.a, plot.b, ncol = 2, labels = c("(a)", "(b)"))
last_row <- ggarrange(NULL, plot.c, NULL,
  ncol = 3, labels = c("", "(c)", ""),
  widths = c(1, 2, 1)
)

ggarrange(first_row, last_row, ncol = 1, heights = c(0.9, 1))
