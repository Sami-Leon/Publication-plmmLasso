source("EM.joint.U2.R") 
source("lasso_inference.r")

lambda.grid <- round(exp(seq(log(0.1), log(1 * 0.0001),
                             length.out = 10
)), digits = 4)

gamma.grid <- c(
  0.000001, 0.0000001, 0.00000001, 0.000000001
)


grid.param = expand.grid(lambda.grid, gamma.grid,  1:100)

# To minimize runtime, run this code in parallel
for (i in 1:nrow(grid.param)) {
  
  lambda = grid.param[i,][[1]]
  gam.cste = grid.param[i,][[2]]
  it = grid.param[i,][[3]]
  
  load(paste0("simu.n100.p100/simu_", it, ".RData"))
  Group1 = Group1_simulations[[1]]
  
  load(paste0("simu.n100.p100/F.Bases_", it, ".RData"))
  
  start_time <- Sys.time()
  
  list.MixteNonPara = EM.joint(Y = Group1$Y, series = Group1$series, position = Group1$position,
                               X = subset(Group1, select = -c(Y, series, position)),
                               F.Bases = F.Bases.norm, gam.cste = gam.cste, intercept = T,
                               lambda.grid = lambda, timexgroup = T, SimData = Group1_simulations)
  
  end_time <- Sys.time()
  
  list.MixteNonPara$time = difftime(end_time, start_time, units = "secs")
  list.MixteNonPara$Iter = it
  assign(paste0("list.MixteNonPara_", it, "_", gam.cste, "_", lambda),
         list.MixteNonPara)
  
  save(list = paste0("list.MixteNonPara_", it, "_", gam.cste, "_", lambda), 
       file = paste0("fit.n100.p100/list.MixteNonPara_", it,
                     "_", gam.cste, "_", lambda, ".RData"))
}

# To minimize runtime, run this code in parallel
for (i in 1:nrow(grid.param)) {
  
  lambda = grid.param[i,][[1]]
  gam.cste = grid.param[i,][[2]]
  it = grid.param[i,][[3]]
  
  load(paste0("simu.n500.p100/simu_", it, ".RData"))
  Group1 = Group1_simulations[[1]]
  
  load(paste0("simu.n500.p100/F.Bases_", it, ".RData"))
  
  start_time <- Sys.time()
  
  list.MixteNonPara = EM.joint(Y = Group1$Y, series = Group1$series, position = Group1$position,
                               X = subset(Group1, select = -c(Y, series, position)),
                               F.Bases = F.Bases.norm, gam.cste = gam.cste, intercept = T,
                               lambda.grid = lambda, timexgroup = T, SimData = Group1_simulations)
  
  end_time <- Sys.time()
  
  list.MixteNonPara$time = difftime(end_time, start_time, units = "secs")
  list.MixteNonPara$Iter = it
  assign(paste0("list.MixteNonPara_", it, "_", gam.cste, "_", lambda),
         list.MixteNonPara)
  
  save(list = paste0("list.MixteNonPara_", it, "_", gam.cste, "_", lambda), 
       file = paste0("fit.n500.p100/list.MixteNonPara_", it,
                     "_", gam.cste, "_", lambda, ".RData"))
}

# To minimize runtime, run this code in parallel
for (i in 1:nrow(grid.param)) {
  
  lambda = grid.param[i,][[1]]
  gam.cste = grid.param[i,][[2]]
  it = grid.param[i,][[3]]
  
  load(paste0("simu.n100.p500/simu_", it, ".RData"))
  Group1 = Group1_simulations[[1]]
  
  load(paste0("simu.n100.p500/F.Bases_", it, ".RData"))
  
  start_time <- Sys.time()
  
  list.MixteNonPara = EM.joint(Y = Group1$Y, series = Group1$series, position = Group1$position,
                               X = subset(Group1, select = -c(Y, series, position)),
                               F.Bases = F.Bases.norm, gam.cste = gam.cste, intercept = T,
                               lambda.grid = lambda, timexgroup = T, SimData = Group1_simulations)
  
  end_time <- Sys.time()
  
  list.MixteNonPara$time = difftime(end_time, start_time, units = "secs")
  list.MixteNonPara$Iter = it
  assign(paste0("list.MixteNonPara_", it, "_", gam.cste, "_", lambda),
         list.MixteNonPara)
  
  save(list = paste0("list.MixteNonPara_", it, "_", gam.cste, "_", lambda), 
       file = paste0("fit.n100.p500/list.MixteNonPara_", it,
                     "_", gam.cste, "_", lambda, ".RData"))
}
