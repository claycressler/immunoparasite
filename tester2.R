setwd("~/immunoparasite")
library(Rcpp)
sourceCpp("nested_model2.cpp")

params = c(S1=1000, S2=1000, s1=2000, s2=2000, 
           b1=0.1, b2=0.1, I12=10000, I21=10000, 
           m=0.9, c1=50, c2=130, C1=50, C2=50, 
           bp=4, Kp=300, a=0.004, 
           b=1e-1, K=100, 
           c=3e-3, v=2e-4, v0=1e-5, cv_v=0.5,
           tmax=400, S0=45, I0=5,
           Th1=500,Th2=700,timestep=1)

out <- nested_modelC2(params)
