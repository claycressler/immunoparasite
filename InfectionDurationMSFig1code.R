library(phaseR, quietly=TRUE)
derivs = function(t, y, params) {
  sig1 <- params["sig1"]
  sig2 <- params["sig2"]
  beta1 <- params["beta1"]
  beta2 <- params["beta2"]
  i12 <- params["i12"]
  i21 <- params["i21"]
  
  t1 <- y[1]
  t2 <- y[2]
  
  dt1 <- beta1 + sig1*t1^2/(1+t1^2) * i12/(i12+t2) - t1
  dt2 <- beta2 + sig2*t2^2/(1+t2^2) * i21/(i21+t1) - t2
  
  list(c(dt1, dt2))
}

png(file="Fig1A.png", height=4, width=4, units='in', res=400)
par(mfrow=c(1,1), mar=c(2,2,0.5,0.5), oma=c(0,0,0,0))
params = c(sig1=2, sig2=2, beta1=0.05, beta2=0.05, i12=10, i21=10)
flows <- flowField(derivs, 
                   xlim=c(0,2),
                   ylim=c(0,2),
                   parameters=params,
                   system="two.dim",
                   state.names=c("",""),
                   add=FALSE)
clines <- nullclines(derivs, 
                     xlim=c(0,2),
                     ylim=c(0,2),
                     parameters=params,
                     col=c(1,2),
                     lwd=2,
                     state.names=c("Th1","Th2"),
                     add.legend=FALSE)
legend(x='topright', "A", bty='n', cex=1.25)
eq1 <- findEquilibrium(derivs, y0=c(0,0), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq2 <- findEquilibrium(derivs, y0=c(0,0.6), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq3 <- findEquilibrium(derivs, y0=c(0,1.2), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq4 <- findEquilibrium(derivs, y0=c(0.6,0), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq5 <- findEquilibrium(derivs, y0=c(1.2,0), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
dev.off()

png(file="Fig1B.png", height=4, width=4, units='in', res=400)
par(mfrow=c(1,1), mar=c(2,2,0.5,0.5), oma=c(0,0,0,0))
params = c(sig1=2, sig2=2, beta1=0.08, beta2=0.08, i12=10, i21=10)
flows <- flowField(derivs, 
                   xlim=c(0,2),
                   ylim=c(0,2),
                   parameters=params,
                   system="two.dim",
                   state.names=c("",""),
                   add=FALSE)
clines <- nullclines(derivs, 
                     xlim=c(0,2),
                     ylim=c(0,2),
                     parameters=params,
                     col=c(1,2),
                     lwd=2,
                     state.names=c("Th1","Th2"),
                     add.legend=FALSE)
legend(x='topright', "B", bty='n', cex=1.25)
eq1 <- findEquilibrium(derivs, y0=c(1,1), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq2 <- findEquilibrium(derivs, y0=c(0.7,0.7), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq3 <- findEquilibrium(derivs, y0=c(0.1,0.1), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq4 <- findEquilibrium(derivs, y0=c(0.1,0.6), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq5 <- findEquilibrium(derivs, y0=c(0.1,1.4), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq6 <- findEquilibrium(derivs, y0=c(0.6,0.1), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq7 <- findEquilibrium(derivs, y0=c(1.4,0.1), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
dev.off()

png(file="Fig1C.png", height=4, width=4, units='in', res=400)
par(mfrow=c(1,1), mar=c(2,2,0.5,0.5), oma=c(0,0,0,0))
params = c(sig1=2, sig2=2, beta1=0.11, beta2=0.11, i12=10, i21=10)
flows <- flowField(derivs, 
                   xlim=c(0,2),
                   ylim=c(0,2),
                   parameters=params,
                   system="two.dim",
                   state.names=c("",""),
                   add=FALSE)
clines <- nullclines(derivs, 
                     xlim=c(0,2),
                     ylim=c(0,2),
                     parameters=params,
                     col=c(1,2),
                     lwd=2,
                     state.names=c("Th1","Th2"),
                     add.legend=FALSE)
legend(x='topright', "C", bty='n', cex=1.25)
eq1 <- findEquilibrium(derivs, y0=c(1.2,1.2), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq2 <- findEquilibrium(derivs, y0=c(0.6,1.2), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq3 <- findEquilibrium(derivs, y0=c(1.2,0.6), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq4 <- findEquilibrium(derivs, y0=c(0.1,0.1), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq5 <- findEquilibrium(derivs, y0=c(0.1,0.4), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq6 <- findEquilibrium(derivs, y0=c(0.4,0.1), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq7 <- findEquilibrium(derivs, y0=c(0.1,1.5), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq8 <- findEquilibrium(derivs, y0=c(1.5,0.1), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq9 <- findEquilibrium(derivs, y0=c(0.4,0.4), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
dev.off()

png(file="Fig1D.png", height=4, width=4, units='in', res=400)
par(mfrow=c(1,1), mar=c(2,2,0.5,0.5), oma=c(0,0,0,0))
params = c(sig1=2, sig2=2, beta1=0.14, beta2=0.14, i12=10, i21=10)
flows <- flowField(derivs, 
                   xlim=c(0,2),
                   ylim=c(0,2),
                   parameters=params,
                   system="two.dim",
                   state.names=c("",""),
                   add=FALSE)
clines <- nullclines(derivs, 
                     xlim=c(0,2),
                     ylim=c(0,2),
                     parameters=params,
                     col=c(1,2),
                     lwd=2,
                     state.names=c("Th1","Th2"),
                     add.legend=TRUE)
legend(x='topright', "D", bty='n', cex=1.25)
eq1 <- findEquilibrium(derivs, y0=c(1.2,1.2), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq2 <- findEquilibrium(derivs, y0=c(0.5,0.5), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq3 <- findEquilibrium(derivs, y0=c(0.1,1.5), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq4 <- findEquilibrium(derivs, y0=c(0.5,1.5), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq5 <- findEquilibrium(derivs, y0=c(1.5,0.1), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq6 <- findEquilibrium(derivs, y0=c(1.5,0.5), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
dev.off()

png(file="Fig1E.png", height=4, width=4, units='in', res=400)
par(mfrow=c(1,1), mar=c(2,2,0.5,0.5), oma=c(0,0,0,0))
params = c(sig1=2, sig2=2, beta1=0.08, beta2=0.04, i12=10, i21=10)
flows <- flowField(derivs, 
                   xlim=c(0,2),
                   ylim=c(0,2),
                   parameters=params,
                   system="two.dim",
                   state.names=c("",""),
                   add=FALSE)
clines <- nullclines(derivs, 
                     xlim=c(0,2),
                     ylim=c(0,2),
                     parameters=params,
                     col=c(1,2),
                     lwd=2,
                     state.names=c("Th1","Th2"),
                     add.legend=FALSE)
legend(x='topright', "E", bty='n', cex=1.25)
eq1 <- findEquilibrium(derivs, y0=c(0,0), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq2 <- findEquilibrium(derivs, y0=c(0,0.6), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq3 <- findEquilibrium(derivs, y0=c(0,1.2), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq4 <- findEquilibrium(derivs, y0=c(0.6,0), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq5 <- findEquilibrium(derivs, y0=c(1.2,0), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
dev.off()

png(file="Fig1F.png", height=4, width=4, units='in', res=400)
par(mfrow=c(1,1), mar=c(2,2,0.5,0.5), oma=c(0,0,0,0))
params = c(sig1=2, sig2=2, beta1=0.12, beta2=0.06, i12=10, i21=10)
flows <- flowField(derivs, 
                   xlim=c(0,2),
                   ylim=c(0,2),
                   parameters=params,
                   system="two.dim",
                   state.names=c("",""),
                   add=FALSE)
clines <- nullclines(derivs, 
                     xlim=c(0,2),
                     ylim=c(0,2),
                     parameters=params,
                     col=c(1,2),
                     lwd=2,
                     state.names=c("Th1","Th2"),
                     add.legend=FALSE)
legend(x='topright', "F", bty='n', cex=1.25)
eq1 <- findEquilibrium(derivs, y0=c(0,0), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq2 <- findEquilibrium(derivs, y0=c(0,0.6), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq3 <- findEquilibrium(derivs, y0=c(0,1.2), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq4 <- findEquilibrium(derivs, y0=c(0.6,0), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq5 <- findEquilibrium(derivs, y0=c(1.2,0), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq6 <- findEquilibrium(derivs, y0=c(0.6,0.6), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq7 <- findEquilibrium(derivs, y0=c(0.7,1.0), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
dev.off()

png(file="Fig1G.png", height=4, width=4, units='in', res=400)
par(mfrow=c(1,1), mar=c(2,2,0.5,0.5), oma=c(0,0,0,0))
params = c(sig1=2, sig2=2, beta1=0.14, beta2=0.07, i12=10, i21=10)
flows <- flowField(derivs, 
                   xlim=c(0,2),
                   ylim=c(0,2),
                   parameters=params,
                   system="two.dim",
                   state.names=c("",""),
                   add=FALSE)
clines <- nullclines(derivs, 
                     xlim=c(0,2),
                     ylim=c(0,2),
                     parameters=params,
                     col=c(1,2),
                     lwd=2,
                     state.names=c("Th1","Th2"),
                     add.legend=FALSE)
legend(x='topright', "G", bty='n', cex=1.25)
eq1 <- findEquilibrium(derivs, y0=c(0.2,0.5), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq2 <- findEquilibrium(derivs, y0=c(0.4,0.6), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq3 <- findEquilibrium(derivs, y0=c(0.2,1.2), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq4 <- findEquilibrium(derivs, y0=c(0.4,1.0), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq5 <- findEquilibrium(derivs, y0=c(1.2,0), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
dev.off()


png(file="Fig1H.png", height=4, width=4, units='in', res=400)
par(mfrow=c(1,1), mar=c(2,2,0.5,0.5), oma=c(0,0,0,0))
params = c(sig1=2, sig2=2, beta1=0.15, beta2=0.075, i12=10, i21=10)
flows <- flowField(derivs, 
                   xlim=c(0,2),
                   ylim=c(0,2),
                   parameters=params,
                   system="two.dim",
                   state.names=c("",""),
                   add=FALSE)
clines <- nullclines(derivs, 
                     xlim=c(0,2),
                     ylim=c(0,2),
                     parameters=params,
                     col=c(1,2),
                     lwd=2,
                     state.names=c("Th1","Th2"),
                     add.legend=FALSE)
legend(x='topright', "H", bty='n', cex=1.25)
eq1 <- findEquilibrium(derivs, y0=c(0.3,1.5), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq2 <- findEquilibrium(derivs, y0=c(0.5,1.4), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
eq3 <- findEquilibrium(derivs, y0=c(1.2,0), parameters=params, system="two.dim", plot.it=TRUE, summary=FALSE)
dev.off()




