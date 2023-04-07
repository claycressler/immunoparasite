library(tidyverse)
library(deSolve)
library(rootSolve)
library(parallel)

simulate_model = function(pars) {
  
  dyn.load("dimensionless_model.so")
  
  dimensionless_modelJ = function(t, y, parms) {
    
    b1 = parms[1]
    b2 = parms[2]
    s1 = parms[3]
    s2 = parms[4]
    i12 = parms[5]
    i21 = parms[6]
    k1 = parms[7]
    k2 = parms[8]
    c1 = parms[9]
    c2 = parms[10]
    r = parms[11]
    a = parms[12]
    
    y1dot = b1 + c1*y[3]/(k1+y[3]) + s1*y[1]*y[1]/(1+y[1]*y[1]) * i12/(i12+y[2]) - y[1];
    y2dot = b2 + c2*y[3]/(k2+y[3]) + s2*y[2]*y[2]/(1+y[2]*y[2]) * i21/(i21+y[1]) - y[2];
    y3dot = r*y[3]*(1-y[3]) - a*y[3]*y[2];
    
    return(list(c(y1dot,y2dot,y3dot)))
  }
  
  iij = unname(pars["iij"])
  rel = unname(pars["rel"])
  totalI = unname(pars["totalI"])
  Th1.0 = unname(pars["Th1.0"])
  P.0 = unname(pars["P.0"])
  
  params = c(b1=0.000111111, b2=.000111111, 
             s1=2.22222, s2=2.22222, 
             i12=iij, i21=iij, 
             k1=1/6, k2=1/6, 
             c1=(1-rel)*0.111111, c2=(1+rel)*0.111111, 
             r=3.33333, a=3.33333)
  y0 = c(Th1=Th1.0, Th2=totalI-Th1.0, P=P.0)
  times = seq(0,100)
  out = ode(y0, times, func="derivs", parms=params, dllname="dimensionless_model", initfunc="initmod", atol=1e-14, rtol=1e-14)
  while((abs(diff(tail(out[,2],2))) > 1e-12 | abs(diff(tail(out[,3],2))) > 1e-12 | abs(diff(tail(out[,4],2))) > 1e-12) & tail(out[,1],1) < 5000) {
    #print(paste("more simulation required for", rel,totalI,Th1.0,P.0))
    times = seq(tail(out[,1],1),tail(out[,1],1)+100,1)
    y0 = c(Th1=tail(out[,"Th1"],1),Th2=tail(out[,"Th2"],1),P=tail(out[,"P"],1))
    out = ode(y0, times, func="derivs", parms=params, dllname="dimensionless_model", initfunc="initmod", atol=1e-14, rtol=1e-14)
    if (out[101,4] < 0) {
      #print("generating negative values"); print(rel); print(iij); print(y0); 
      out[101,4] = 0} # catch negatives
  }
  ## Classify the equilibrium reached (so long as the final t < 5000)
  if (tail(out[,1],1) < 5000) {
    ## Calculate the stability of the equilibrium
    equilStab = eigen(jacobian.full(y=as.numeric(tail(out[,2:4],1)),
                                    func=dimensionless_modelJ, 
                                    parms=params))$values
    if (all(Re(equilStab) < 0)) {## this equilibrium is stable, classify it
      if (tail(out[,"Th1"],1) < 0.1 & tail(out[,"Th2"],1) > 1 & tail(out[,"P"],1) < 1e-9)
        outcome = 1 ## Th2 polarization
      else if (tail(out[,"Th1"],1) > 0.8 & tail(out[,"Th2"],1) > 0.8 & tail(out[,"P"],1) < 0.2)
        outcome = 2 ## High Th1+Th2
      else if (tail(out[,"Th1"],1) > 0.9 & tail(out[,"Th2"],1) < 0.3 & tail(out[,"P"],1) > 0.6)
        outcome = 3 ## Th1 polarization
      else outcome = 4 ## Low Th1+Th2
    } else outcome = 5 ## Saddle
  } else outcome = 6 ## Cycling
  return(c(totalI,
           Th1.0,
           P.0,
           tail(out[,"time"],1),
           round(tail(out[,"Th1"],1),3),
           round(tail(out[,"Th2"],1),3),
           round(tail(out[,"P"],1),3),
           Th1.0/totalI,
           round((1+rel)/(1-rel),1),
           outcome))
}


for (iij in c(5,1,0.5,0.25)) {
  print(iij)
  parameters = vector(mode='list', length=6*10*101*101)
  i = 1
  for (rel in seq(0,0.5,0.1)) {
    for (totalI in seq(0.2,2,0.2)) {
      for (Th1.0 in seq(0,totalI,totalI/100)) {
        for (P.0 in seq(0.01,0.2,0.19/100)) {
          parameters[[i]] = c(iij=iij, rel=rel, totalI=totalI,Th1.0=Th1.0,P.0=P.0)
          i = i+1
        }
      }
    }
  }
  mclapply(parameters,
           function(p) simulate_model(p),
           mc.cores=12) -> results
  results2 = do.call("rbind.data.frame",results)
  colnames(results2) = c("initT","initTh1","initP","t","Th1","Th2","P","Th1bias","Th2bias","outcome") 
  mutate(results2,
         Th2bias=paste0("Th2=",Th2bias,"*Th1")) -> results2
  results2$Th2bias = factor(results2$Th2bias,
                           levels=c("Th2=1*Th1",
                                    "Th2=1.2*Th1",
                                    "Th2=1.5*Th1",
                                    "Th2=1.9*Th1", 
                                    "Th2=2.3*Th1", 
                                    "Th2=3*Th1"))
  saveRDS(results2, file=paste0("Dimensionless_model_results_across_parameters_initials_iij=",iij,".RDS"))
}


for (iij in c(20,10,5,1,0.5,0.25)) {
  print(iij)
  parameters = vector(mode='list', length=6*201*201*4)
  i = 1
  for (rel in seq(0,0.5,0.1)) {
    for (totalI in seq(0.2,2,1.8/200)) {
      for (Th1.0 in seq(0,totalI,totalI/200)) {
        for (P.0 in c(0.05,0.1,0.15,0.2)) {
          parameters[[i]] = c(iij=iij, rel=rel, totalI=totalI,Th1.0=Th1.0,P.0=P.0)
          i = i+1
        }
      }
    }
  }
  mclapply(parameters,
           function(p) simulate_model(p),
           mc.cores=12) -> results
  results2 = do.call("rbind.data.frame",results)
  colnames(results2) = c("initT","initTh1","initP","t","Th1","Th2","P","Th1bias","Th2bias","outcome") 
  mutate(results2,
         Th2bias=paste0("Th2=",Th2bias,"*Th1")) -> results2
  results2$Th2bias = factor(results2$Th2bias,
                            levels=c("Th2=1*Th1",
                                     "Th2=1.2*Th1",
                                     "Th2=1.5*Th1",
                                     "Th2=1.9*Th1", 
                                     "Th2=2.3*Th1", 
                                     "Th2=3*Th1"))
  
  saveRDS(results2, file=paste0("Dimensionless_model_results_across_parameters_initials_iij=",iij,"_fixed_P0.RDS"))
}



