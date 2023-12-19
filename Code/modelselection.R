source("EM.joint.U2.R") 

lambda.grid <- round(exp(seq(log(0.1), log(1 * 0.0001),
                             length.out = 10
)), digits = 4)

gamma.grid <- c(
  0.000001, 0.0000001, 0.00000001, 0.000000001
)


grid.param = expand.grid(lambda.grid, gamma.grid)

# To minimize runtime, run this code in parallel
for (i in 1:nrow(grid.param)) {
  
  lambda = grid.param[i,][[1]]
  gam.cste = grid.param[i,][[2]]

  start_time <- Sys.time()
  
  list.MixteNonPara = EM.joint(Y = Group1$Y, series = Group1$series, position = Group1$position,
                               X = subset(Group1, select = -c(Y, series, position)),
                               F.Bases = F.Bases.norm, gam.cste = gam.cste, intercept = T,
                               lambda.grid = lambda, timexgroup = T)
  
  end_time <- Sys.time()
  
  list.MixteNonPara$time = difftime(end_time, start_time, units = "secs")

}
