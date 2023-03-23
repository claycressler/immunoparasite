
x1 = read.csv("Dimensionless_results_across_immune_parasite_initials_and_total_immune_equal_parasite_activation_weak_cross-inhibition.csv", header=FALSE)
x2 = read.csv("Dimensionless_results_across_immune_parasite_initials_and_total_immune_equal_parasite_activation_moderate_cross-inhibition.csv", header=FALSE)
x3 = read.csv("Dimensionless_results_across_immune_parasite_initials_and_total_immune_biased_parasite_activation_weak_cross-inhibition.csv", header=FALSE)
x4 = read.csv("Dimensionless_results_across_immune_parasite_initials_and_total_immune_biased_parasite_activation_moderate_cross-inhibition.csv", header=FALSE)
y1 = as.data.frame(matrix(unlist(x1), ncol=4, byrow=TRUE))
y2 = as.data.frame(matrix(unlist(x2), ncol=4, byrow=TRUE))
y3 = as.data.frame(matrix(unlist(x3), ncol=4, byrow=TRUE))
y4 = as.data.frame(matrix(unlist(x4), ncol=4, byrow=TRUE))
colnames(y1) <- colnames(y2) <- colnames(y3) <- colnames(y4) <- c("totalI","initialP","initialTh1bias","equilP")
y1$equilP <- round(y1$equil,3)
y2$equilP <- round(y2$equil,3)
y3$equilP <- round(y3$equil,3)
y4$equilP <- round(y4$equil,3)
y1$totalI <- round(y1$totalI,1)
y2$totalI <- round(y2$totalI,1)
y3$totalI <- round(y3$totalI,1)
y4$totalI <- round(y4$totalI,1)

mybreaks1 = c(sort(unique(y1$equilP)),1) # th2 polarized, coactivation, naive, th1 polarized
mybreaks2 = c(sort(unique(y2$equilP)),1) # th2 polarized, coactivation, naive, th1 polarized
mybreaks3 = c(sort(unique(y3$equilP)),1) # th2 polarized, th1 polarized
mybreaks4 = c(sort(unique(y4$equilP)),1) # th2 polarized, th1 polarized

library(ggplot2)
ggplot(y1, aes(initialTh1bias,initialP,z=equilP)) + 
  geom_contour_filled(breaks=mybreaks1) + 
  scale_fill_discrete(type=c("black","lightgrey","blue","red"), labels=c("Th2","Co","Low","Th1")) + 
  facet_wrap(~totalI) + theme_bw()

ggplot(y2, aes(initialTh1bias,initialP,z=equilP)) + 
  geom_contour_filled(breaks=mybreaks2) + 
  scale_fill_discrete(type=c("black","lightgrey","blue","red"), labels=c("Th2","Co","Low","Th1")) + 
  facet_wrap(~totalI) + theme_bw()

ggplot(y3, aes(initialTh1bias,initialP,z=equilP)) + 
  geom_contour_filled(breaks=mybreaks2) + 
  scale_fill_discrete(type=c("black","red"), labels=c("Th2","Th1")) + 
  facet_wrap(~totalI) + theme_bw()

ggplot(y4, aes(initialTh1bias,initialP,z=equilP)) + 
  geom_contour_filled(breaks=mybreaks4) + 
  scale_fill_discrete(type=c("black","red"), labels=c("Th2","Th1")) + 
  facet_wrap(~totalI) + theme_bw()


library(deSolve)
model <- function(t, y, pars) {
  b1 = unname(params["b1"])
  b2 = unname(params["b2"])
  c1 = unname(params["c1"])
  c2 = unname(params["c2"])
  k1 = unname(params["k1"])
  k2 = unname(params["k2"])
  s1 = unname(params["s1"])
  s2 = unname(params["s2"])
  i12 = unname(params['i12'])  
  i21 = unname(params['i21'])
  r = unname(params["r"])
  a = unname(params["a"])
  Th1 = y[1]
  Th2 = y[2]
  P = y[3]
  
  dTh1dt = b1 + c1*P/(k1+P) + s1*Th1^2/(1+Th1^2) * i12/(i12+Th2) - Th1
  dTh2dt = b2 + c2*P/(k2+P) + s2*Th2^2/(1+Th2^2) * i21/(i21+Th1) - Th2
  dPdt = r*P*(1-P) - a*P*Th2
  
  return(list(c(dTh1dt,dTh2dt,dPdt)))
  
}
 

params = c(b1=0.000111111, b2=.000111111, 
           s1=2.22222, s2=2.22222, 
           i12=10, i21=10, 
           k1=1/6, k2=1/6, 
           c1=0.111111, c2=0.111111, 
           r=3.33333, a=3.33333)
y0 = c(Th1=0.5, Th2=0.5, P=0.1)
times = seq(0,1000,1)
## Takes 24 seconds to simulate 1000 timesteps 1000 times using uncompiled code
system.time(sapply(1:1000, function(i) ode(y0, times, model, params)))
dyn.load("dimensionless_model.dll")
## Takes less than 1 second to do the same simulation, but using compiled code instead
system.time(sapply(1:1000, function(i) ode(y0, times, func="derivs", parms=params, dllname="dimensionless_model", initfunc="initmod")))

for (iij in c(20, 10, 5, 1, 0.5)) {
  tic <- Sys.time()
  print(iij)
  results = vector(mode='list', length=6*10*201*191)
  i = 1
  for (rel in seq(0,0.5,0.1)) {
    for (totalI in seq(0.2,2,0.2)) {
      for (Th1.0 in seq(0,totalI,totalI/200)) {
        for (P.0 in seq(0.01,0.2,0.001)) {
          params = c(b1=0.000111111, b2=.000111111, 
                     s1=2.22222, s2=2.22222, 
                     i12=iij, i21=iij, 
                     k1=1/6, k2=1/6, 
                     c1=(1-rel)*0.111111, c2=(1+rel)*0.111111, 
                     r=3.33333, a=3.33333)
          y0 = c(Th1=Th1.0, Th2=totalI-Th1.0, P=P.0)
          times = seq(0,100)
          out = ode(y0, times, func="derivs", parms=params, dllname="dimensionless_model", initfunc="initmod")
          
          while(abs(diff(tail(out[,4],2))) > 1e-12 & tail(out[,1],1) < 5000) {
            #print(paste("more simulation required for", rel,totalI,Th1.0,P.0))
            times = seq(tail(out[,1],1),tail(out[,1],1)+100,1)
            y0 = tail(out[,2:4],1)
            out = ode(y0, times, func="derivs", parms=params, dllname="dimensionless_model", initfunc="initmod")
          }
          results[[i]] <- c(rel,totalI,Th1.0,P.0,as.numeric(tail(out[,1:4],1)))
          i = i+1
        }
      }
    }
  }
  print(Sys.time()-tic)
  saveRDS(results, file=paste0("Dimensionless_model_results_across_parameters_initials_iij=",iij,".RDS"))
}

## for rel=0.5, there is the emergence of an unstable, cyclic equilibrium
## this is a high coactivation equilibrium, and the system cycles around it
## so biologically, it is entirely infeasible because the parasite burden gets
## very, very low (like 1e-4).
## But it's very cool that one of the equilibria undergoes a Hopf bifurcation!
##

results = readRDS(file="Dimensionless_model_results_across_parameters_initials_iij=10.RDS")

results = do.call("rbind.data.frame",results)
colnames(results) = c("Th2activation","initT","Th1Th2bias","initP","t","Th1","Th2","P")
results$P = round(results$P,3)
results$initT = round(results$initT,1)
results$outcome = rep(0,nrow(results))
## 1. Th2 polarization and parasite extinction
results$outcome[which(results$P==0)] = 1
## 2. High coactivation and near-extinction
results$outcome[which(results$t==5000)] = 2
results$outcome[which(results$P > 0 & results$P < 0.3 & results$Th1 > 0.9 & results$Th2 > 0.9)] = 2
## 3. Th1 polarization and chronic infection
results$outcome[which(results$P > 0.7 & results$Th1 > 1 & results$Th2 < 0.3)] = 3
## 4. Low coactivation and chronic infection
results$outcome[which(results$P > 0.7 & results$Th1 < 0.3 & results$Th2 < 0.3)] = 4

saveRDS(results, 


ggplot(subset(results, initT==0.8 & Th2activation < 0.6),aes(Th1Th2bias,initP,z=outcome)) + 
  geom_contour_filled(breaks=c(1,2,3,4,5)) + 
  scale_fill_discrete(type=c("black","lightgrey","red","blue"), name="", labels=c("High Th2","High Th1+2", "High Th1", "Low Th1+2")) + 
  facet_wrap(~Th2activation) + xlab("Th2 <--------> Th1") + ylab("Initial parasite") + theme_bw()



range(results[which(results[,5]==5000),"P"])
range(results[which(results[,5]==5000),"Th1"])
range(results[which(results[,5]==5000),"Th2"])

