library(Rcpp)

sourceCpp("nested_model.cpp")
source("nested_model.R")

params = c(S1=1000, S2=1000, s1=2000, s2=2000, 
           b1=0, b2=0, I12=10000, I21=10000, 
           m=0.9, c1=50, c2=130, C1=50, C2=50, 
           bp=1,Kp=300, a=0.004, 
           c=1e-2, v=1e-2, v0=1e-3, cv_v=0.5,
           tmax=50, S0=95, I0=5)


system.time(nested_modelC(params) -> outC)
system.time(nested_modelR(params) -> outR)
