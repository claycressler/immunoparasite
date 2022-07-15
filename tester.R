setwd("~/immunoparasite")
library(Rcpp)
sourceCpp("nested_model.cpp")
source("nested_model.R")
source("within_host_model.R")

params = c(S1=1000, S2=1000, s1=2000, s2=2000, 
           b1=0.1, b2=0.1, I12=10000, I21=10000, 
           m=0.9, c1=50, c2=130, C1=50, C2=50, 
           bp=4, Kp=300, a=0.004, 
           b=1e-1, K=100, 
           c=1e-2, v=1e-3, v0=1e-3, cv_v=0.5,
           tmax=10, S0=10, I0=5,
           Th1=800,Th2=500)

## What are the deterministic dynamics of the within-host model if everyone starts at (T1=0,T2=0)?
out <- lapply(seq(10,100,10), function(p) ode(y=c(T1=0,T2=0,P=p), times=seq(0,50,0.1), within_host_model, params))
## Infections are always chronic
lapply(out, function(o) o[501,])

## okay, what if you assume a much higher initial values of T1 and T2?
out <- lapply(seq(10,100,10), function(p) ode(y=c(T1=800,T2=500,P=p), times=seq(0,50,0.1), within_host_model, params))
## Here you get a chronic infection if dose is 10 or 20, and an acute infection if dose is 30 or larger
lapply(out, function(o) o[501,])

## What about an initial dose of 30 exactly? What do the dynamics look like?
det = ode(y=c(T1=800,T2=500,P=30), times=seq(0,500,0.1), within_host_model, params)
#plot(det[,c(1,4)], type='l') ## should see an initial rise, and then squashed!

## what about here?
out <- lapply(seq(10,100,10), function(p) ode(y=c(T1=500,T2=800,P=p), times=seq(0,50,0.1), within_host_model, params))
## Infections are always acute
lapply(out, function(o) o[501,])

o1 = ode(y=c(T1=500,T2=800,P=100), times=seq(0,50,0.1), within_host_model, params)

## Set up the nested model such that we expect variability in infection duration based on the infecting dose 
if(!file.exists("Nested_model_out_variable_outcome.RDS")) {
  params = c(S1=1000, S2=1000, s1=2000, s2=2000, 
             b1=0.1, b2=0.1, I12=10000, I21=10000, 
             m=0.9, c1=50, c2=130, C1=50, C2=50, 
             bp=4, Kp=300, a=0.004, 
             b=1e-1, K=100, 
             c=3e-3, v=2e-4, v0=1e-5, cv_v=0.5,
             tmax=200, S0=45, I0=5,
             Th1=800,Th2=500)
  out <- vector(mode='list', length=20)
  for (i in 1:20) {
    print(i)
    tIn <- Sys.time()
    out[[i]] <- nested_modelC(params)
    tOut <- Sys.time()
    print(tOut-tIn)
  }
  saveRDS(out, file="Nested_model_out_variable_outcome.RDS")
}

if(!file.exists("Nested_model_out_chronic_outcome.RDS")) {
  print("chronic model")
  library(parallel)
  params = c(S1=1000, S2=1000, s1=2000, s2=2000, 
             b1=0.1, b2=0.1, I12=10000, I21=10000, 
             m=0.9, c1=50, c2=130, C1=50, C2=50, 
             bp=4, Kp=300, a=0.004, 
             b=1e-1, K=100, 
             c=3e-3, v=2e-4, v0=1e-5, cv_v=0.5,
             tmax=200, S0=45, I0=5,
             Th1=0,Th2=0)
   
  out <- vector(mode='list', length=20)
  for (i in 1:20) {
    print(i)
    tIn <- Sys.time()
    out[[i]] <- nested_modelC(params)
    tOut <- Sys.time()
    print(tOut-tIn)
  }
  saveRDS(out, file="Nested_model_out_chronic_outcome.RDS")
}


# if(!file.exists("Nested_model_out_acute_outcome.RDS")) {
#   print("acute model")
#   library(parallel)
#   params = c(S1=1000, S2=1000, s1=2000, s2=2000, 
#              b1=0.1, b2=0.1, I12=10000, I21=10000, 
#              m=0.9, c1=50, c2=130, C1=50, C2=50, 
#              bp=4, Kp=300, a=0.004, 
#              b=1e-1, K=100, 
#              c=3e-3, v=2e-4, v0=1e-5, cv_v=0.5,
#              tmax=10, S0=45, I0=5,
#              Th1=500,Th2=700)
#   mclapply(1:50, function(i) nested_modelC(params), mc.cores=10) -> o
#   
#   for (i in 1:50) o = nested_modelC(params)
#   
#   source("nested_model.R")
#   nested_modelR(params)
#   
#   out <- vector(mode='list', length=20)
#   for (i in 1:20) {
#     print(i)
#     tIn <- Sys.time()
#     out[[i]] <- nested_modelC(params)
#     tOut <- Sys.time()
#     print(tOut-tIn)
#   }
#   saveRDS(out, file="Nested_model_out_acute_outcome.RDS")
# }


# o = nested_modelC(params)
# 
# ## track the fates of each of the initially infected individuals - did some sustain an acute, and others a chronic, infection?
# lapply(10:14, function(i) unlist(lapply(o[[1]], function(o2) o2[o2[,5]==i,3])))
# ## In this run, one of the original infections was acute; the other 4 were chronic
# t(sapply(o[[1]], function(o2) o2[o2[,5]==12,]))[200:300,]
# ## if you look at all of the infected hosts at the final timestep, what are their states
# 
# 
