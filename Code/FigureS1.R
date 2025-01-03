source("CreateSim.R")
source("set.zero.constraint.R")
source("CreateBases.R")

f2 <- function(t, A, omega) {
  out <- as.vector(A * sin(2 * pi * t / omega))

  out[out == 0] <- -10
  out <- out + 5
  return(out)
}

simulate_S1 <- function(N = 50, n.mvnorm = 100, grouped = T, seed,
                        timepoints = 3:5, nonpara.inter = T, sample_from, cst_ni) {
  set.seed(seed)
  if (nonpara.inter) {
    A <- c(1, 1.5)
    omega <- c(60, 110)
  } else {
    A <- c(1, 1)
    omega <- c(60, 60)
  }

  f0mean <- mean(f(sample_from, A[1], omega[1]))

  Y <- NULL
  out <- NULL
  phi <- rnorm(N, 0, 0.5)

  f.val <- NULL

  for (i in 1:N) {
    if (cst_ni) {
      ni <- timepoints
    } else {
      ni <- sample(timepoints, 1)
    }


    if (grouped) {
      theta <- c(3, 2, 0)
    } else {
      theta <- c(0, 2, 0)
    }

    group <- rep(sample(c(0, 1), 1), ni)

    X1 <- rep(rnorm(1, 1, 0.3), ni)

    eps <- rnorm(ni, 0, 0.05)

    t <- sort(sample(sample_from, ni, replace = F))

    if (group[1] == 0) {
      out <- rbind(out, cbind(
        rep(i, ni), t, phi[i] + f(t, A[1], omega[1]) - f0mean + eps,
        group, X1
      ))

      f.val <- c(f.val, f(t, A[1], omega[1]) - f0mean)
    } else {
      out <- rbind(out, cbind(
        rep(i, ni), t, phi[i] + f2(t, A[2], omega[2]) + eps,
        group, X1
      ))

      f.val <- c(f.val, f2(t, A[2], omega[2]))
    }
  }

  X <- MASS::mvrnorm(nrow(out), rep(0, n.mvnorm - 1), cst_cor(n.mvnorm - 1, 0))

  out <- cbind(out, X)

  colnames(out) <- c("series", "position", "Y", "Group", paste0("X", 1:(ncol(X) + 1)))

  Group1 <- out

  Group1[, "X2"] <- Group1[, "Group"] * Group1[, "X1"]
  Group1[, "Y"] <- Group1[, "Y"] + Group1[, c("Group", "X1", "X2"), drop = F] %*% theta

  Group1 <- as.data.frame(Group1)

  phi <- rep(phi, table(Group1$series))

  Group1 <- Group1[order(Group1$series, Group1$position), ]

  f.val <- f.val[order(Group1$series, Group1$position)]

  phi <- phi[order(Group1$series, Group1$position)]

  return(list(Group1, phi, f.val))
}

CreationBasesZero <- function(position, keep = NULL) {
  Ff <- c()
  numKnots <- 4
  n <- max(position)

  for (i in seq(1, 100, 1)) {
    Ff <- cbind(Ff, sin(2 * pi * i * position / n))
  }

  x <- position / n
  Fx <- c()
  rg <- seq(0.1, 2, 0.01)
  for (i in 1:length(rg)) {
    Fx <- cbind(Fx, x^rg[i])
  }

  F.tot <- cbind(Ff, Fx)

  if (!is.null(keep)) {
    F.tot <- F.tot[, keep]
  }

  F.Bases <- c()
  F.Bases.tot <- c()
  F.Bases.tot <- Bases.NonNulles(F.tot)
  Num.Bases.Pres <- F.Bases.tot$PresentBases
  F.Bases <- F.Bases.tot$Fh.P
  return(list(F.Bases = F.Bases, Num.Bases.Pres = Num.Bases.Pres))
}

Simu <- simulate_S1(
  N = 500, n.mvnorm = 2, grouped = T,
  seed = 1, timepoints = 3:5, nonpara.inter = T,
  sample_from = seq(0, 52, 13), cst_ni = F
)

Data <- Simu[[1]][, 1:6]
position <- Data$position

Dico.norm <- CreationBasesZero(position)
F.Bases.norm <- Dico.norm$F.Bases
Num.Bases.Pres.norm <- Dico.norm$Num.Bases.Pres


fit.plmm.S1 <- plmmlasso.zero(
  Y = Data$Y, series = Data$series, position = Data$position,
  X = Data[, c("Group", "X1", "X2")], F.Bases = F.Bases.norm,
  gam.cste = 0.1, intercept = TRUE, lambda.grid = 0.000000001,
  timexgroup = TRUE
)

Data$phi <- Simu[[2]]
Data$f <- Simu[[3]]

Data$X <- as.matrix(Data[, 4:6]) %*% cbind(c(3, 2, 0))

bluemean_X <- Data %>%
  group_by(Group) %>%
  mutate(mean = mean(X)) %>%
  ungroup()

bluemean_U2 <- Data %>%
  group_by(Group) %>%
  mutate(mean = mean(phi)) %>%
  ungroup()

Data$blueline <- Data$f + bluemean_U2$mean + bluemean_X$mean

mean_blueline <- Data %>%
  group_by(Group, position) %>%
  mutate(mean = mean(blueline)) %>%
  ungroup()

Data$blueline <- mean_blueline$mean

Data$F.fit <- fit.plmm.S1$Res.F$out.F$F.fit
Data$X.fit <- fit.plmm.S1$Res.F$X.fit
Data$U2 <- rep(fit.plmm.S1$U2$U2, table(Data$series))

mean_X <- Data %>%
  group_by(Group) %>%
  mutate(mean = mean(X.fit)) %>%
  ungroup()

mean_U2 <- Data %>%
  group_by(Group) %>%
  mutate(mean = mean(U2)) %>%
  ungroup()

redline <- fit.plmm.S1$Res.F$out.F$F.fit + mean_X$mean + mean_U2$mean
Data$redline <- redline

mean_redline <- Data %>%
  group_by(Group, position) %>%
  mutate(mean = mean(redline)) %>%
  ungroup()

Data$redline <- mean_redline$mean

group_label <- c("0" = "Group1", "1" = "Group2")
g1 <- ggplot(data = Data, aes(x = position, y = Y)) +
  geom_line(aes(x = position, y = Y, group = series)) +
  geom_line(aes(x = position, y = blueline, color = "truth"), data = Data, size = 1.3) +
  geom_line(aes(x = position, y = redline, color = "estimate"),
    data = Data, size = 1.3
  ) +
  xlab("Time") +
  scale_color_manual(name = "", values = c("truth" = "blue", "estimate" = "red")) +
  scale_x_continuous(breaks = c(0, 13, 26, 39, 52)) +
  theme(
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "bottom"
  ) +
  facet_grid(. ~ Group, labeller = as_labeller(group_label))

mean.0 <- attr(scale(unique(Data[Data$Group == 0, ]$F.fit),
  scale = F
), "scaled:center")
mean.1 <- attr(scale(unique(Data[Data$Group == 1, ]$F.fit),
  scale = F
), "scaled:center")

Data$F.fit.adj <- Data$F.fit
Data[Data$Group == 0, ]$F.fit.adj <- Data[Data$Group == 0, ]$F.fit - mean.0
Data[Data$Group == 1, ]$F.fit.adj <- Data[Data$Group == 1, ]$F.fit - mean.1


g2 <- ggplot(data = Data, aes(x = position, y = Y)) +
  geom_line(aes(x = position, y = F.fit.adj, color = "mean-zero", linetype = "mean-zero"), size = 1.3) +
  geom_line(aes(x = position, y = F.fit, color = "set-to-zero", linetype = "set-to-zero"), size = 1.3) +
  scale_color_manual(name = "", values = c("mean-zero" = "red", "set-to-zero" = "#32CD32")) +
  scale_linetype_manual(name = "", values = c("mean-zero" = "solid", "set-to-zero" = "dashed")) +
  facet_grid(. ~ Group, labeller = as_labeller(group_label)) +
  ylim(c(-10, 15)) +
  xlab("Time") +
  scale_x_continuous(breaks = c(0, 13, 26, 39, 52)) +
  theme(
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "bottom"
  )

ggarrange(g1, g2, ncol = 2, labels = c("(a)", "(b)"), font.label = list(size = 16))
