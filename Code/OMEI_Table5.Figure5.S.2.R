library(openxlsx)

source("CreateBases.R")
source("posi.R")
source("modelselection.R")
source("test.nonlinear.functions.R")
source("CreateSim.R")

OMEI <- read.xlsx("OMEI_synth.xlsx")

lambda.grid <- round(exp(seq(log(0.2), log(1 * 0.001),
  length.out = 20
)), digits = 3)[3:5]

gamma.grid <- c(
  0.0000001, 0.00000001, 0.000000001,
  0.0000000001, 0.00000000001
)[4:5]


best.BIC <- select.plmm(
  data = OMEI,
  gamma = gamma.grid,
  lambda = lambda.grid,
  intercept = TRUE,
  timexgroup = TRUE
)

redline <- best.BIC$Res.F$out.F$F.fit + rep(best.BIC$U2$U2, best.BIC$ni) +
  best.BIC$Res.F$X.fit

OMEI$redline <- redline
redline <- OMEI %>%
  group_by(Group, position) %>%
  mutate(mean = mean(redline)) %>%
  ungroup()

OMEI$redline <- redline$mean

OMEI$F.fit <- best.BIC$Res.F$out.F$F.fit

g1 <- ggplot(data = OMEI, aes(x = position, y = Y)) +
  geom_line(aes(x = position, y = Y, group = series, col = factor(Group)),
    linewidth = 0.3, linetype = "longdash"
  ) +
  geom_line(aes(x = position, y = redline, group = factor(Group), col = factor(Group)),
    data = OMEI, linewidth = 1.4, linetype = "solid"
  ) +
  theme(text = element_text(size = 12)) +
  xlab("time (Months)") +
  ylab("Saliva_Ca_cfu_Ln") +
  scale_color_manual(
    name = "", values = c("0" = "blue", "1" = "red"),
    labels = c("non-black", "black")
  ) +
  scale_x_continuous(breaks = c(1, 2, 4, 6, 12, 18, 24))

g2 <- ggplot(data = best.BIC$Res.F$out.F, aes(x = position, y = F.fit)) +
  geom_line(aes(x = position, y = F.fit, group = factor(Group), col = factor(Group)),
    data = best.BIC$Res.F$out.F, linewidth = 1.4, linetype = "solid"
  ) +
  theme(text = element_text(size = 12)) +
  ylim(c(-1.5, 15)) +
  xlab("time (Months)") +
  ylab("f(time)") +
  scale_color_manual(
    name = "", values = c("0" = "blue", "1" = "red"),
    labels = c("non-black", "Black")
  ) +
  scale_x_continuous(breaks = c(1, 2, 4, 6, 12, 18, 24))

OMEI$redline <- NULL
OMEI$F.fit <- NULL

test <- debias.plmm(OMEI, best.BIC)

test <- test$adapt.debias
df.omei <- data.frame(
  test$beta.hat, test$bhat,
  test$pv, test$ci
)

df.omei <- cbind(rownames(df.omei), df.omei)
rownames(df.omei) <- NULL
colnames(df.omei) <- c(
  "Variable", "Estimate", "Debiased", "p-value",
  "CI lower", "CI upper"
)

df.omei <- df.omei[test$pv < 0.05, ]

df.omei <- df.omei[order(abs(df.omei$Estimate), decreasing = T), ]

df.omei[, c(2:3, 5:6)] <- lapply(df.omei[, c(2:3, 5:6)], round, 2)
df.omei[, 4] <- roundpval(df.omei[, 4])

rownames(df.omei) <- NULL

df.omei

test.nonlinear.functions <- f.test(OMEI, best.BIC, n = 100, predicted = FALSE)

test.nonlinear.functions$df.CI

g3 <- test.nonlinear.functions$plot.CI


first_row <- ggarrange(g1, g2,
  ncol = 2, labels = c("(a)", "(b)"), legend = "top",
  common.legend = T
)

last_row <- ggarrange(NULL, g3, NULL,
  ncol = 3, labels = c("", "(c)", ""),
  widths = c(1, 2, 1)
)

ggarrange(first_row, last_row, ncol = 1, heights = c(1, 0.9))


test.nonlinear.functions$plot.boot
