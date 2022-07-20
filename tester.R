setwd("~/immunoparasite")
library(Rcpp)
library(parallel)
## this is the version of the model where initial Th1 and Th2 are fixed
## based on 'exploring_within-host_dynamics.R', I know that for a fixed Th2ness (e.g., Th1=700, Th2=500)
## changing dose and bp will alter the dynamics. In particular, there is only a narrow range of bp values
## for which dose has a variable affect on infection outcome. For low or high bp, infections are always chronic;
## for moderate bp they can be variable or always acute (which is kind of crazy, actually, and is 
## something I need to explore in much more detail with a bifurcation diagram.
sourceCpp("nested_model_fixed_Th2.cpp")

## this is the version of the model where initial Th2ness varies from 400-700, with totalT fixed at 1200
sourceCpp("nested_model_vary_Th2.cpp")

## Based on 'exploring_within-host_dynamics.R', I know that all infections will become
## chronic at bp=8. Let's allow that to be the maximum possible bp, but set the initial 
## virulence so that the realized bp=4 (e.g., v=v0, since realized bp is bp*v/(v+v0))
params = c(S1=1000, S2=1000, s1=2000, s2=2000, 
           b1=0.1, b2=0.1, I12=10000, I21=10000, 
           m=0.9, c1=50, c2=130, C1=50, C2=50, 
           bp=8, Kp=300, a=0.004, 
           b=1e-1, K=100, 
           c=3e-3, v=1e-4, v0=1e-4, cv_v=0.5,
           tmax=10, S0=95, I0=5,
           Th1=700,Th2=500,timestep=1)
out <- nested_model_fixed_Th2(params)
out2 <- nested_model_vary_Th2(params)

params["tmax"] <- 400
tIn <- Sys.time()
mclapply(1:50, 
         function(x) out = nested_model_vary_Th2(params),
         mc.cores=10) -> out
tOut <- Sys.time()
print(tOut-tIn)
saveRDS(out, file="Nested_model_variable_dose_variable_Th2_variable_outcome.RDS")

tIn <- Sys.time()
mclapply(1:50, 
         function(x) out = nested_model_fixed_Th2(params),
         mc.cores=10) -> out
tOut <- Sys.time()
print(tOut-tIn)
saveRDS(out, file="Nested_model_variable_dose_fixed_Th2_variable_outcome.RDS")

## modify fixed initial Th2ness so that infections are always chronic
params = c(S1=1000, S2=1000, s1=2000, s2=2000, 
           b1=0.1, b2=0.1, I12=10000, I21=10000, 
           m=0.9, c1=50, c2=130, C1=50, C2=50, 
           bp=8, Kp=300, a=0.004, 
           b=1e-1, K=100, 
           c=3e-3, v=1e-4, v0=1e-4, cv_v=0.5,
           tmax=400, S0=95, I0=5,
           Th1=0,Th2=0,timestep=1)
tIn <- Sys.time()
mclapply(1:50, 
         function(x) out = nested_model_fixed_Th2(params),
         mc.cores=10) -> out
tOut <- Sys.time()
print(tOut-tIn)
saveRDS(out, file="Nested_model_variable_dose_fixed_Th2_chronic_outcome.RDS")

## modify fixed initial Th2ness so that infections are always acute (unless bp gets very large)
params = c(S1=1000, S2=1000, s1=2000, s2=2000, 
           b1=0.1, b2=0.1, I12=10000, I21=10000, 
           m=0.9, c1=50, c2=130, C1=50, C2=50, 
           bp=8, Kp=300, a=0.004, 
           b=1e-1, K=100, 
           c=3e-3, v=1e-4, v0=1e-4, cv_v=0.5,
           tmax=400, S0=95, I0=5,
           Th1=500,Th2=700,timestep=1)
tIn <- Sys.time()
mclapply(1:50, 
         function(x) out = nested_model_fixed_Th2(params),
         mc.cores=10) -> out
tOut <- Sys.time()
print(tOut-tIn)
saveRDS(out, file="Nested_model_variable_dose_fixed_Th2_acute_outcome.RDS")


