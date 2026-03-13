setwd("~/immunoparasite")
library(Rcpp)
library(parallel)
## this is the version of the model where initial Th1 and Th2 are fixed
## based on 'exploring_within-host_dynamics.R', I know that for a fixed Th2ness (e.g., Th1=700, Th2=500)
## changing dose and bp will alter the dynamics. In particular, there is only a narrow range of bp values
## for which dose has a variable affect on infection outcome. For low or high bp, infections are always chronic;
## for moderate bp they can be variable or always acute (which is kind of crazy, actually, and is 
## something I need to explore in much more detail with a bifurcation diagram.

## sourceCpp("nested_model_fixed_Th2.cpp")

## this is the version of the model where initial Th2ness varies from 400-700, with totalT fixed at 1200
sourceCpp("nested_model_vary_Th2.cpp")

for (th2 in seq(400,800,50)) {
  params = c(S1=1000, S2=1000, s1=2000, s2=2000,
             b1=0.1, b2=0.1, I12=10000, I21=10000,
             m=0.9, c1=50, c2=130, C1=50, C2=50,
             bp=8, Kp=300, a=0.004,
             b=1e-1, K=100,
             c=3e-3, v=1e-4, v0=1e-4, cv_v=0.5,
             tmax=400, S0=95, I0=5,
             minTh2=th2,maxTh2=th2,timestep=1)
  print(th2)
  tIn <- Sys.time()
  mclapply(1:50,
           function(x) out = nested_model_vary_Th2(params),
           mc.cores=15) -> out
  tOut <- Sys.time()
  print(tOut-tIn)
  saveRDS(out, file=paste0("Nested_model_variable_dose_fixed_Th2=",th2,".RDS"))
}

png(file="Epidemiological_dynamics_as_fixed_Th2_varies.png", height=9, width=9,  units='in', res=300)
par(mfrow=c(3,3), mar=c(3.5,3.5,0.5,0.5), oma=rep(0,4))
for (th2 in seq(400,800,50)) {
  out2 <- readRDS(file=paste0("Nested_model_variable_dose_fixed_Th2=",th2,".RDS"))
  ## Y axis limits
  ymax <- ceiling((lapply(out2, function(o) max(o[[2]][,3])) %>% unlist %>% max)/10)*10
  plot.new()
  plot.window(xlim=c(0,400), ylim=c(0,ymax))
  axis(1); axis(2); box('plot')
  for (i in 1:50) lines(out2[[i]][[2]][,c(1,3)], col=ifelse(out2[[i]][[2]][401,3]>0,1,2))
  mtext(side=1, line=2.5, "Time")
  mtext(side=2, line=2.5, "No. infected")
  legend(x='bottomright', legend=paste0("Th2=",th2), bty='n')
}
dev.off()

png(file="Evolutionary_dynamics_as_fixed_Th2_varies.png", height=9, width=9,  units='in', res=300)
par(mfrow=c(3,3), mar=c(3.5,3.5,0.5,0.5), oma=rep(0,4))
for (th2 in seq(400,800,50)) {
  out2 <- readRDS(file=paste0("Nested_model_variable_dose_fixed_Th2=",th2,".RDS"))
  plot.new()
  plot.window(xlim=c(0,400), ylim=c(0,1e-3))
  axis(1); axis(2); box('plot')
  for (i in 1:50) lines(out2[[i]][[2]][,c(1,6)], col=gray(0.7))
  lines(1:401, sapply(out2, function(o) o[[2]][,6]) %>% apply(., 1, mean), lwd=2)
  mtext(side=1, line=2.5, "Time")
  mtext(side=2, line=2.5, "Per-parasite virulence")
  legend(x='bottomright', legend=paste0("Th2=",th2), bty='n')
}
dev.off()


for (th2 in seq(400,550,50)) {
  params = c(S1=1000, S2=1000, s1=2000, s2=2000,
             b1=0.1, b2=0.1, I12=10000, I21=10000,
             m=0.9, c1=50, c2=130, C1=50, C2=50,
             bp=8, Kp=300, a=0.004,
             b=1e-1, K=100,
             c=3e-3, v=1e-4, v0=1e-4, cv_v=0.5,
             tmax=400, S0=95, I0=5,
             minTh2=th2,maxTh2=1200-th2,timestep=1)
  print(th2)
  tIn <- Sys.time()
  mclapply(1:50,
           function(x) out = nested_model_vary_Th2(params),
           mc.cores=12) -> out
  tOut <- Sys.time()
  print(tOut-tIn)
  saveRDS(out, file=paste0("Nested_model_variable_dose_5050_bias_variable_Th2=",th2,".RDS"))
}

