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

for (th2 in c(400,500,600,610,620,630,640,650,700,800)) {
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

png(file="Epidemiological_dynamics_as_fixed_Th2_varies_3-25.png", height=9, width=9,  units='in', res=300)
par(mfrow=c(4,3), mar=c(3.5,3.5,0.5,0.5), oma=rep(0,4))
for (th2 in c(400,500,600,610,620,630,640,650,700,800)) {
  out2 <- readRDS(file=paste0("Nested_model_variable_dose_fixed_Th2=",th2,".RDS"))
  ## Y axis limits
  ymax <- ceiling((lapply(out2, function(o) max(o[[2]][,3])) %>% unlist %>% max)/10)*10
  plot.new()
  plot.window(xlim=c(0,400), ylim=c(0,ymax))
  axis(1); axis(2); box('plot')
  for (i in 1:50) lines(out2[[i]][[2]][,c(1,3)], col=ifelse(out2[[i]][[2]][401,3]>0,1,2))
  mtext(side=1, line=2.5, "Time")
  mtext(side=2, line=2.5, "No. infected")
  legend(x='bottomright', legend=paste0("Init Th2=",th2), bty='n')
}
dev.off()

png(file="Mean_burden_dynamics_as_fixed_Th2_varies_3-25.png", height=9, width=9,  units='in', res=300)
par(mfrow=c(4,3), mar=c(3.5,3.5,0.5,0.5), oma=rep(0,4))
for (th2 in c(400,500,600,610,620,630,640,650,700,800)) {
  out2 <- readRDS(file=paste0("Nested_model_variable_dose_fixed_Th2=",th2,".RDS"))
  plot.new()
  plot.window(xlim=c(0,400), ylim=c(0,210))
  axis(1); axis(2); box('plot')
  for (i in 1:50) lines(1:401,lapply(out2[[i]][[1]], function(o) mean(o[which(o[,3]>0),3])) %>% unlist)
  mtext(side=1, line=2.5, "Time")
  mtext(side=2, line=2.5, "Avg. burden")
  legend(x='bottomright', legend=paste0("Init Th2=",th2), bty='n')
}
dev.off()

## Peak burden during each individual infection, and how that changes over time
## Calculate the peak burden and the time of the peak for each individual infection
par(mfrow=c(4,3), mar=c(3.5,3.5,0.5,0.5), oma=rep(0,4))
for (th2 in c(400,500,600,610,620,630,640,650,700,800)) {
  out2 <- readRDS(file=paste0("Nested_model_variable_dose_fixed_Th2=",th2,".RDS"))
  
  
  plot.new()
  plot.window(xlim=c(0,400), ylim=c(0,210))
  axis(1); axis(2); box('plot')
  for (i in 1:50) lines(1:401,lapply(out2[[i]][[1]], function(o) mean(o[which(o[,3]>0),3])) %>% unlist)
  mtext(side=1, line=2.5, "Time")
  mtext(side=2, line=2.5, "Avg. burden")
  legend(x='bottomright', legend=paste0("Init Th2=",th2), bty='n')
}
dev.off()


png(file="Evolutionary_dynamics_as_fixed_Th2_varies_3-25.png", height=9, width=9,  units='in', res=300)
par(mfrow=c(4,3), mar=c(3.5,3.5,0.5,0.5), oma=rep(0,4))
for (th2 in c(400,500,600,610,620,630,640,650,700,800)) {
  out2 <- readRDS(file=paste0("Nested_model_variable_dose_fixed_Th2=",th2,".RDS"))
  plot.new()
  plot.window(xlim=c(0,400), ylim=c(0,1e-3))
  axis(1); axis(2); box('plot')
  for (i in 1:50) {
    vseq = out2[[i]][[2]][,6]
    if (any(is.na(vseq))) {
      out2[[i]][[2]][which(is.na(vseq)):length(vseq),6] = NA
      lines(out2[[i]][[2]][,c(1,6)], col='pink')
    }
    else lines(out2[[i]][[2]][,c(1,6)], col=gray(0.7))
  }
  lines(1:401, sapply(out2, function(o) o[[2]][,6]) %>% apply(., 1, function(q) mean(q,na.rm=T)), lwd=2)
  mtext(side=1, line=2.5, "Time")
  mtext(side=2, line=2.5, "Per-parasite virulence")
  legend(x='bottomright', legend=paste0("Init Th2=",th2), bty='n')
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

png(file="Epidemiological_dynamics_as_min_Th2_varies_5050_bias.png", height=9, width=9,  units='in', res=300)
par(mfrow=c(2,2), mar=c(3.5,3.5,0.5,0.5), oma=rep(0,4))
for (th2 in seq(400,550,50)) {
  out2 <- readRDS(file=paste0("Nested_model_variable_dose_5050_bias_variable_Th2=",th2,".RDS"))
  ## Y axis limits
  ymax <- ceiling((lapply(out2, function(o) max(o[[2]][,3])) %>% unlist %>% max)/10)*10
  plot.new()
  plot.window(xlim=c(0,400), ylim=c(0,ymax))
  axis(1); axis(2); box('plot')
  for (i in 1:50) lines(out2[[i]][[2]][,c(1,3)], col=ifelse(out2[[i]][[2]][401,3]>0,1,2))
  mtext(side=1, line=2.5, "Time")
  mtext(side=2, line=2.5, "No. infected")
  legend(x='bottomright', legend=paste0("Init Th2=",th2,"-",1200-th2), bty='n')
}
dev.off()

png(file="Burden_dynamics_as_min_Th2_varies_5050_bias.png", height=9, width=9,  units='in', res=300)
par(mfrow=c(2,2), mar=c(3.5,3.5,0.5,0.5), oma=rep(0,4))
for (th2 in seq(400,550,50)) {
  out2 <- readRDS(file=paste0("Nested_model_variable_dose_5050_bias_variable_Th2=",th2,".RDS"))
  plot.new()
  plot.window(xlim=c(0,400), ylim=c(0,210))
  axis(1); axis(2); box('plot')
  for (i in 1:50) lines(1:401,lapply(out2[[i]][[1]], function(o) mean(o[which(o[,3]>0),3])) %>% unlist)
  mtext(side=1, line=2.5, "Time")
  mtext(side=2, line=2.5, "Avg. burden")
  legend(x='bottomright', legend=paste0("Init Th2=",th2,"-",1200-th2), bty='n')
}
dev.off()


png(file="Evolutionary_dynamics_as_min_Th2_varies_5050_bias.png", height=9, width=9,  units='in', res=300)
par(mfrow=c(2,2), mar=c(3.5,3.5,0.5,0.5), oma=rep(0,4))
for (th2 in seq(400,550,50)) {
  out2 <- readRDS(file=paste0("Nested_model_variable_dose_5050_bias_variable_Th2=",th2,".RDS"))
  plot.new()
  plot.window(xlim=c(0,400), ylim=c(0,5e-4))
  axis(1); axis(2); box('plot')
  for (i in 1:50) {
    vseq = out2[[i]][[2]][,6]
    if (any(is.na(vseq))) {
      out2[[i]][[2]][which(is.na(vseq)):length(vseq),6] = NA
      lines(out2[[i]][[2]][,c(1,6)], col='pink')
    }
    else lines(out2[[i]][[2]][,c(1,6)], col=gray(0.7))
  }
  lines(1:401, sapply(out2, function(o) o[[2]][,6]) %>% apply(., 1, function(q) mean(q,na.rm=T)), lwd=2)
  mtext(side=1, line=2.5, "Time")
  mtext(side=2, line=2.5, "Per-parasite virulence")
  legend(x='bottomright', legend=paste0("Init Th2=",th2,"-",1200-th2), bty='n')
}
dev.off()


for (th2 in seq(200,700,100)) {
  params = c(S1=1000, S2=1000, s1=2000, s2=2000,
             b1=0.1, b2=0.1, I12=10000, I21=10000,
             m=0.9, c1=50, c2=130, C1=50, C2=50,
             bp=8, Kp=300, a=0.004,
             b=1e-1, K=100,
             c=3e-3, v=1e-4, v0=1e-4, cv_v=0.5,
             tmax=400, S0=95, I0=5,
             minTh2=th2,maxTh2=800,timestep=1)
  print(th2)
  tIn <- Sys.time()
  mclapply(1:50,
           function(x) out = nested_model_vary_Th2(params),
           mc.cores=12) -> out
  tOut <- Sys.time()
  print(tOut-tIn)
  saveRDS(out, file=paste0("Nested_model_variable_dose_min_Th2=",th2,"_3-25.RDS"))
}

for (th2 in seq(200,750,50)) {
  params = c(S1=1000, S2=1000, s1=2000, s2=2000,
             b1=0.1, b2=0.1, I12=10000, I21=10000,
             m=0.9, c1=50, c2=130, C1=50, C2=50,
             bp=8, Kp=300, a=0.004,
             b=1e-1, K=100,
             c=3e-3, v=1e-4, v0=1e-4, cv_v=0.5,
             tmax=600, S0=95, I0=5,
             minTh2=th2,maxTh2=800,timestep=1)
  print(th2)
  tIn <- Sys.time()
  mclapply(1:50,
           function(x) out = nested_model_vary_Th2(params),
           mc.cores=12) -> out
  tOut <- Sys.time()
  print(tOut-tIn)
  saveRDS(out, file=paste0("Nested_model_variable_dose_min_Th2=",th2,"_3-25.RDS"))
}

for (th2 in seq(525,575,25)) {
  params = c(S1=1000, S2=1000, s1=2000, s2=2000,
             b1=0.1, b2=0.1, I12=10000, I21=10000,
             m=0.9, c1=50, c2=130, C1=50, C2=50,
             bp=8, Kp=300, a=0.004,
             b=1e-1, K=100,
             c=3e-3, v=1e-4, v0=1e-4, cv_v=0.5,
             tmax=600, S0=95, I0=5,
             minTh2=th2,maxTh2=800,timestep=1)
  print(th2)
  tIn <- Sys.time()
  mclapply(1:50,
           function(x) out = nested_model_vary_Th2(params),
           mc.cores=12) -> out
  tOut <- Sys.time()
  print(tOut-tIn)
  saveRDS(out, file=paste0("Nested_model_variable_dose_min_Th2=",th2,"_3-25.RDS"))
}



png(file="Epidemiological_dynamics_as_min_Th2_varies.png", height=9, width=9,  units='in', res=300)
par(mfrow=c(5,3), mar=c(3.5,3.5,0.5,0.5), oma=rep(0,4))
for (th2 in c(seq(200,500,50),seq(525,600,25),seq(650,750,50))) {
  out2 <- readRDS(file=paste0("Nested_model_variable_dose_min_Th2=",th2,"_3-25.RDS"))
  ## Y axis limits
  plot.new()
  plot.window(xlim=c(0,max(out2[[1]][[2]][,1])), ylim=c(0,110))
  axis(1); axis(2); box('plot')
  nextinct=0
  for (i in 1:50) {
    lines(out2[[i]][[2]][,c(1,3)], col=ifelse(as.numeric(tail(out2[[i]][[2]],1)[,3])>0,1,2))
    nextinct = nextinct + ifelse(as.numeric(tail(out2[[i]][[2]],1)[,3])>0,0,1)
  }
  mtext(side=1, line=2.5, "Time")
  mtext(side=2, line=2.5, "No. infected")
  legend(x='topleft', legend=paste0("Init Th2=",th2,"-",800), bty='n', text.col='blue')
  legend(x='topright', legend=paste0("Prob(fadeout)=",nextinct/50), bty='n', text.col='red')
}
dev.off()

png(file="Mean_burden_thru_time_as_min_Th2_varies.png", height=9, width=9,  units='in', res=300)
par(mfrow=c(5,3), mar=c(3.5,3.5,0.5,0.5), oma=rep(0,4))
for (th2 in c(seq(200,500,50),seq(525,600,25),seq(650,750,50))) {
  out2 <- readRDS(file=paste0("Nested_model_variable_dose_min_Th2=",th2,"_3-25.RDS"))
  plot.new()
  plot.window(xlim=c(0,max(out2[[1]][[2]][,1])), ylim=c(0,225))
  axis(1); axis(2); box('plot')
  for (i in 1:50) lines(out2[[i]][[2]][,1],lapply(out2[[i]][[1]], function(o) mean(o[which(o[,3]>0),3])) %>% unlist)
  mtext(side=1, line=2.5, "Time")
  mtext(side=2, line=2.5, "Avg. burden")
  legend(x='topleft', legend=paste0("Init Th2=",th2,"-",800), bty='n', text.col='blue')
  
}
dev.off()

png(file="Evolutionary_dynamics_as_min_Th2_varies.png", height=9, width=9,  units='in', res=300)
par(mfrow=c(5,3), mar=c(3.5,3.5,0.5,0.5), oma=rep(0,4))
for (th2 in c(seq(200,500,50),seq(525,600,25),seq(650,750,50))) {
  out2 <- readRDS(file=paste0("Nested_model_variable_dose_min_Th2=",th2,"_3-25.RDS"))
  plot.new()
  plot.window(xlim=c(0,max(out2[[1]][[2]][,1])), ylim=c(0,2e-3))
  axis(1); axis(2); box('plot')
  for (i in 1:50) {
    vseq = out2[[i]][[2]][,6]
    if (any(is.na(vseq))) {
      out2[[i]][[2]][which(is.na(vseq)):length(vseq),6] = NA
      lines(out2[[i]][[2]][,c(1,6)], col='red')
    }
    else lines(out2[[i]][[2]][,c(1,6)], col=gray(0.7))
  }
  lines(out2[[i]][[2]][,1], sapply(out2, function(o) o[[2]][,6]) %>% apply(., 1, function(q) mean(q,na.rm=T)), lwd=2)
  mtext(side=1, line=2.5, "Time")
  mtext(side=2, line=2.5, "Per-parasite virulence")
  legend(x='topleft', legend=paste0("Init Th2=",th2,"-",800), bty='n', text.col='blue')
  
}
dev.off()

png(file="Peak_probability_dynamics_as_min_Th2_varies.png", height=9, width=9,  units='in', res=300)
par(mfrow=c(5,3), mar=c(3.5,3.5,0.5,0.5), oma=rep(0,4))
for (th2 in c(seq(200,500,50),seq(525,600,25),seq(650,750,50))) {
  out2 <- readRDS(file=paste0("Nested_model_variable_dose_min_Th2=",th2,"_3-25.RDS"))
  ## For individuals whose infections started at different time intervals, what was the peak infection
  for (t in seq(0, round(max(duration_data$ends))-150,50)) {
    filter(duration_data, starts >= t, starts < t+50) %>% 
      mutate(event_code=ifelse(event=="clearance",1,ifelse(event=="death",2,0))) %>%
      with(., cuminc(ftime=durations, fstatus=event_code)) -> ci
    prob_clear = c(prob_clear, tail(ci[[1]]$est,1))
    var_clear = c(var_clear, tail(ci[[1]]$var,1))
  }
  plot(seq(0, round(max(duration_data$ends))-150,50), prob_clear, type='l', lwd=2, xaxt='n', xlab='', ylab='', ylim=c(0,1))
  mtext(side=1, line=2.5, "Infection interval")
  mtext(side=2, line=2.5, "Clearance probability")
  axis(1, 
       tick=TRUE, 
       at=seq(0, round(max(duration_data$ends))-150, 50), 
       labels=paste0("[",seq(0, round(max(duration_data$ends))-150, 50),",",seq(50, round(max(duration_data$ends))-100, 50),")"))
  legend(x='topleft', legend=paste0("Init Th2=",th2,"-",800), bty='n', text.col='blue')
}
dev.off()

## Compute information about the duration of each individual infection and the ultimate outcome of infection (peak infection, clearance, death, neither)
for (th2 in c(seq(200,500,50),seq(525,600,25),seq(650,750,50))) {
  print(th2)
  out2 <- readRDS(file=paste0("Nested_model_variable_dose_min_Th2=",th2,"_3-25.RDS"))
  out3 = vector(mode='list', length=length(out2))
  for (j in 1:length(out2)) {
    for (i in 1:length(out2[[j]][[1]])) {
      as.data.frame(out2[[j]][[1]][[i]], col.names=c("Th1","Th2","P","v","id")) %>%
        mutate(time=i) -> out2[[j]][[1]][[i]]
    }
    out3[[j]] = do.call("rbind.data.frame", out2[[j]][[1]])
  }
  ## Calculate infection duration for every individual infection
  duration_data = vector(mode='list', length=length(out3))
  for (j in 1:length(out3)) {
    inds = unique(out3[[j]]$V5)
    infection_data = vector(mode='list', length=length(inds))
    for (i in 1:length(inds)) {
      ind = inds[i]
      oo = filter(out3[[j]], V5==ind)
      ## Duration of any infections for this individual
      if (any(oo$V3 > 0)) {
        times = oo$time[which(oo$V3 > 0)]
        max_age = max(times)
        data.frame(durations = lengths(split(times, cumsum(c(1, diff(times) != 1)))),
                   starts = sapply(split(times, cumsum(c(1, diff(times) != 1))), min),
                   ends = sapply(split(times, cumsum(c(1, diff(times) != 1))), max)) %>%
          mutate(event=ifelse(ends==max_age, "death", 
                              ifelse(max_age<max(out3[[j]]$time), "clearance", "censored"))) -> infection_data[[i]]
      }
    }
    duration_data[[j]] = do.call("rbind.data.frame", infection_data)
  }
  duration_data = do.call("rbind.data.frame", duration_data)
  
  ## durations are calculated as if the timestep over which data was recorded was dt=1
  ## If that is not true, then every time in duration_data needs to be multiplied by dt
  dt = min(diff(out2[[1]][[2]][,1]))
  duration_data %<>% 
    mutate(durations = dt*durations,
           starts = dt*starts,
           ends= dt*ends)
  
  saveRDS(duration_data, paste0("Duration_data_min_Th2=",th2,"_3-25.RDS"))
}
  
## Compute the median duration of infection for infections starting at different times
## This analysis does not distinguish between infections ending due to death versus clearance
## (That I will tackle later.)
library(survival)
png(file="Infection_duration_dynamics_as_min_Th2_varies.png", height=9, width=9,  units='in', res=300)
par(mfrow=c(5,3), mar=c(3.5,3.5,0.5,0.5), oma=rep(0,4))
for (th2 in c(seq(200,500,50),seq(525,600,25),seq(650,750,50))) {
  duration_data = readRDS(paste0("Duration_data_min_Th2=",th2,"_3-25.RDS"))
  mean_duration = c()
  se_duration = c()
  for (t in seq(0, round(max(duration_data$ends))-150,50)) {
    filter(duration_data, starts >= t, starts < t+50) %>% 
      mutate(event_code = ifelse(event!="censored",1,0)) %>% 
      with(., survfit(Surv(time = durations, event = event_code)~1)) -> km_fit
    mean_duration = c(mean_duration, summary(km_fit)$table["rmean"])
    se_duration = c(se_duration, summary(km_fit)$table["se(rmean)"])
  }
  plot(seq(0, round(max(duration_data$ends))-150, 50), mean_duration+se_duration, lwd=2, lty=2, type='l', xaxt='n', ylim=c(min(mean_duration-se_duration),max(mean_duration+se_duration)), xlab="", ylab="")
  mtext(side=1, line=2.5, "Infection interval")
  mtext(side=2, line=2.5, "Mean infection duration")
  axis(1, 
       tick=TRUE, 
       at=seq(0, round(max(duration_data$ends))-150, 50), 
       labels=paste0("[",seq(0, round(max(duration_data$ends))-150, 50),",",seq(50, round(max(duration_data$ends))-100, 50),")"))
  lines(seq(0, round(max(duration_data$ends))-150, 50), mean_duration, lwd=2)
  lines(seq(0, round(max(duration_data$ends))-150, 50), mean_duration-se_duration, lty=2, lwd=2)
  legend(x='topleft', legend=paste0("Init Th2=",th2,"-",800), bty='n', text.col='blue')
}
dev.off()
 
## What if instead I want to look at probabilities of clearance or death by a particular age of infection
duration_data = readRDS(paste0("Duration_data_min_Th2=",th2,"_3-25.RDS"))
duration_data %>% 
  mutate(event_code=ifelse(event=="clearance",1,ifelse(event=="death",2,0))) %>%
  filter(., starts >=0, starts < 50) %>%
  with(., cuminc(ftime=durations, fstatus=event_code)) %>%
  plot(., 
     lwd=2,
     col=c('blue','red'),
     xlab="Time since infection",
     ylab="Cumulative incidence",
     curvlab=c("Clearance","Death"))

duration_data = readRDS(paste0("Duration_data_min_Th2=",th2,"_3-25.RDS"))
duration_data %>% 
  mutate(event_code=ifelse(event=="clearance",1,ifelse(event=="death",2,0))) %>%
  filter(., starts > 200, starts < 250) %>%
  with(., cuminc(ftime=durations, fstatus=event_code)) %>%
  plot(., 
       lwd=2,
       col=c('blue','red'),
       xlab="Time since infection",
       ylab="Cumulative incidence",
       curvlab=c("Clearance","Death"))

## Clearance probability always appears to asymptote rather quickly, whereas the probability of death asymptotes much more slowly, especially when immunity is th1-biased
## It does asymptote rather quickly when immunity is th2-biased
## Regardless, I think we can just compute the long-term probability of clearance, since death is really the only 
## other outcome that is possible, it just depends how long it takes to get there!
library(cmprsk)
png(file="Clearance_probability_dynamics_as_min_Th2_varies.png", height=9, width=9,  units='in', res=300)
par(mfrow=c(5,3), mar=c(3.5,3.5,0.5,0.5), oma=rep(0,4))
for (th2 in c(seq(200,500,50),seq(525,600,25),seq(650,750,50))) {
  duration_data = readRDS(paste0("Duration_data_min_Th2=",th2,"_3-25.RDS"))
  prob_clear = var_clear = c()
  for (t in seq(0, round(max(duration_data$ends))-150,50)) {
    filter(duration_data, starts >= t, starts < t+50) %>% 
      mutate(event_code=ifelse(event=="clearance",1,ifelse(event=="death",2,0))) %>%
      with(., cuminc(ftime=durations, fstatus=event_code)) -> ci
    prob_clear = c(prob_clear, tail(ci[[1]]$est,1))
    var_clear = c(var_clear, tail(ci[[1]]$var,1))
  }
  plot(seq(0, round(max(duration_data$ends))-150,50), prob_clear, type='l', lwd=2, xaxt='n', xlab='', ylab='', ylim=c(0,1))
  mtext(side=1, line=2.5, "Infection interval")
  mtext(side=2, line=2.5, "Clearance probability")
  axis(1, 
       tick=TRUE, 
       at=seq(0, round(max(duration_data$ends))-150, 50), 
       labels=paste0("[",seq(0, round(max(duration_data$ends))-150, 50),",",seq(50, round(max(duration_data$ends))-100, 50),")"))
  legend(x='topleft', legend=paste0("Init Th2=",th2,"-",800), bty='n', text.col='blue')
}
dev.off()
