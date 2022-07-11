library(phaseR, quietly=TRUE)
library(magrittr)
library(deSolve)
library(animation)

full_model = function(t, y, params) {
  s1 <- params["s1"]
  s2 <- params["s2"]
  b1 <- params["b1"]
  b2 <- params["b2"]
  I12 <- params["I12"]
  I21 <- params["I21"]
  S1 <- params["S1"]
  S2 <- params["S2"]
  m <- params["m"]
  c1 <- params["c1"]
  c2 <- params["c2"]
  C1 <- params["C1"]
  C2 <- params["C2"]
  bp <- params["bp"]
  Kp <- params["Kp"]
  a <- params["a"]
  
  T1 <- y[1]
  T2 <- y[2]
  P <- y[3]
  
  dT1 <- b1 + c1*P/(C1+P) + s1*T1^2/(S1^2+T1^2) * I12/(I12+T2) - m*T1
  dT2 <- b2 + c2*P/(C2+P) + s2*T2^2/(S2^2+T2^2) * I21/(I21+T1) - m*T2
  dP <- bp*P*(1-P/Kp) - a*T2*P
  
  list(c(dT1, dT2, dP))
}

T2P_cline_model = function(t, y, params) {
  s1 <- params["s1"]
  s2 <- params["s2"]
  b1 <- params["b1"]
  b2 <- params["b2"]
  I12 <- params["I12"]
  I21 <- params["I21"]
  S1 <- params["S1"]
  S2 <- params["S2"]
  m <- params["m"]
  c1 <- params["c1"]
  c2 <- params["c2"]
  C1 <- params["C1"]
  C2 <- params["C2"]
  bp <- params["bp"]
  Kp <- params["Kp"]
  a <- params["a"]
  T1 <- params["T1"]
  
  T2 <- y[1]
  P <- y[2]
  
  dT2 <- b2 + c2*P/(C2+P) + s2*T2^2/(S2^2+T2^2) * I21/(I21+T1) - m*T2
  dP <- bp*P*(1-P/Kp) - a*T2*P
  
  list(c(dT2, dP))
}

T1T2_cline_model = function(t, y, params) {
  s1 <- params["s1"]
  s2 <- params["s2"]
  b1 <- params["b1"]
  b2 <- params["b2"]
  I12 <- params["I12"]
  I21 <- params["I21"]
  S1 <- params["S1"]
  S2 <- params["S2"]
  m <- params["m"]
  c1 <- params["c1"]
  c2 <- params["c2"]
  C1 <- params["C1"]
  C2 <- params["C2"]
  bp <- params["bp"]
  Kp <- params["Kp"]
  a <- params["a"]
  P <- params["P"]
  
  T1 <- y[1]
  T2 <- y[2]
  
  dT1 <- b1 + c1*P/(C1+P) + s1*T1^2/(S1^2+T1^2) * I12/(I12+T2) - m*T1
  dT2 <- b2 + c2*P/(C2+P) + s2*T2^2/(S2^2+T2^2) * I21/(I21+T1) - m*T2
  
  list(c(dT1, dT2))
}

## Variation in parasite dose can generate differences in infection duration
## With this parameter set, low doses can become chronic, but high doses get cleared due to stronger chronicity-promoting feedbacks
params = c(S1=1000, S2=1000, s1=2000, s2=2000, b1=0.1, b2=0.1, I12=10000, I21=10000, m=0.9, c1=50, c2=130, C1=50, C2=50, bp=4,Kp=300, a=0.004)
## Start Th1 biased
out10 <- ode(y=c(T1=800,T2=500,P=10), times=seq(0,50,0.1), full_model, params)
out80 <- ode(y=c(T1=800,T2=500,P=80), times=seq(0,50,0.1), full_model, params)

ani.options(interval=0.2, ani.dev = "png", ani.type = "png", ani.height = 800, ani.width = 1200)
saveVideo({ 
  for (i in c(1:101,seq(103,201,2))) {
    params2 <- c(params, T1=unname(out10[i,2]))
    par(mfrow=c(1,2), mar=c(4.5,4.5,1,1), oma=rep(0.5,4))

    params2 <- c(params, P=unname(out10[i,4]))
    flows <- flowField(T1T2_cline_model, 
                       xlim=c(0,2100),
                       ylim=c(0,2100),
                       parameters=params2,
                       system="two.dim",
                       state.names=c("Th1 response","Th2 response"),
                       add=FALSE, 
                       cex.axis=1.5,
                       cex.lab=1.5)
    clines <- nullclines(T1T2_cline_model, 
                         xlim=c(0,2100),
                         ylim=c(0,2100),
                         parameters=params2,
                         col=c(1,2),
                         lwd=4,
                         state.names=c("Th1 response","Th2 response"),
                         add.legend=FALSE,
                         cex.lab=1.5)
    legend(x='topleft', legend=c("Th1 nullclines", "Th2 nullclines"), cex=1.4, lwd=3, col=c(1,2), bty='n')
    points(out10[1,2], out10[1,3], pch=21, bg='orange', cex=4)
    lines(out10[1:i,c(2,3)], col='orange', lwd=4)
    
    plot(out10[1:i,c(1,4)], type='l', xlab='Time', ylab='Parasite biomass', lwd=4, cex.axis=1.5, cex.lab=1.5, xlim=c(0,20), ylim=c(0,250))
    
  }}, video.name=paste0("movie1_for_Andrea_3-21-2022.mp4"), other.opts = "-pix_fmt yuv420p -b 300k")
  
ani.options(interval=0.2, ani.dev = "png", ani.type = "png", ani.height = 800, ani.width = 1200)
saveVideo({ 
  for (i in c(1:101,seq(103,201,2))) {
    params2 <- c(params, T1=unname(out10[i,2]))
    par(mfrow=c(1,2), mar=c(4.5,4.5,1,1), oma=rep(0.5,4))
    
    params2 <- c(params, P=unname(out10[i,4]))
    flows <- flowField(T1T2_cline_model, 
                       xlim=c(0,2100),
                       ylim=c(0,2100),
                       parameters=params2,
                       system="two.dim",
                       state.names=c("Th1 response","Th2 response"),
                       add=FALSE, 
                       cex.axis=1.5,
                       cex.lab=1.5)
    clines <- nullclines(T1T2_cline_model, 
                         xlim=c(0,2100),
                         ylim=c(0,2100),
                         parameters=params2,
                         col=c(1,2),
                         lwd=4,
                         state.names=c("Th1 response","Th2 response"),
                         add.legend=FALSE,
                         cex.lab=1.5)
    legend(x='topleft', legend=c("Th1 nullclines", "Th2 nullclines"), cex=1.4, lwd=3, col=c(1,2), bty='n')
    points(out80[1,2], out80[1,3], pch=21, bg='orange', cex=4)
    lines(out80[1:i,c(2,3)], col='orange', lwd=4)
    
    plot(out80[1:i,c(1,4)], type='l', xlab='Time', ylab='Parasite biomass', lwd=4, cex.axis=1.5, cex.lab=1.5, xlim=c(0,20), ylim=c(0,250))
    
  }}, video.name=paste0("movie2_for_Andrea_3-21-2022.mp4"), other.opts = "-pix_fmt yuv420p -b 300k")

    

    
