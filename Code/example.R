source("CreateSim.R")
source("modelselection.R")
source("posi.R")
source("test.functions.R")

# Generating simulated dataset
data.sim = simulate_group_inter(N = 100, n.mvnorm = 100, grouped = T,
                                seed = 12, timepoints = 3:5, nonpara.inter = T,
                                sample_from = seq(0,52,13), cst_ni = F)

lambda.grid <- round(exp(seq(log(0.1), log(1 * 0.0001),
                             length.out = 10
)), digits = 4)

gamma.grid <- c(
  0.000001, 0.0000001, 0.00000001, 0.000000001
)

# Run model on a grid of hyperparameter and retrieve the model with the best BIC
fit.plmm <- select.plmm(
  data = data.sim[[1]], gamma = gamma.grid, lambda = lambda.grid,
  crit = "BIC", intercept = T, timexgroup = T
)

# Visualize overall fit and nonlinear functions 
plot.fit(fit.plmm, data.sim)

# Get debiased fixed-effects and pvalues
posi = debias.plmm(simu = data.sim[[1]], model = fit.plmm)

# Format post-selection inference output
df.posi = cbind(posi$beta.hat, posi$bhat, posi$ci, posi$pv)
colnames(df.posi) = c("Estimate", "Debiased", "CI lower", "CI upper", "p-value")
head(df.posi)

# Inference for the nonlinear functions
test.nonlinear.functions = f.test(data.sim[[1]], fit.plmm, n = 1000)

# Test of equality of the functions H0: f1 = f2
test.nonlinear.functions[[1]]

# Continuous joint confidence bands for the difference between f1 and f2
test.nonlinear.functions[[3]]
