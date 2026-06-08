setwd("~/immunoparasite")
library(Rcpp)
library(parallel)
library(patchwork)
library(tidyverse)
library(magrittr)
library(deSolve)
sourceCpp("nested_model_vary_Th2.cpp")

for (th2 in c(seq(200,500,50),seq(525,600,25),seq(650,750,50))) {
  params = c(S1=1000, S2=1000, s1=2000, s2=2000,
             b1=0.1, b2=0.1, I12=10000, I21=10000,
             m=0.9, c1=50, c2=130, C1=50, C2=50,
             bp=8, Kp=300, a=0.004,
             b=1e-1, K=100,
             c=3e-3, v=1e-4, v0=1e-4, cv_v=0.5,
             tmax=700, S0=95, I0=5,
             minTh2=th2,maxTh2=800,timestep=1)
  print(th2)
  tIn <- Sys.time()
  mclapply(1:50,
           function(x) nested_model_vary_Th2(params),
           mc.cores=12) -> out
  tOut <- Sys.time()
  print(tOut-tIn)
  saveRDS(out, file=paste0("Nested_model_variable_dose_min_Th2=",th2,"_4-9.RDS"))
}

th2=575
params = c(S1=1000, S2=1000, s1=2000, s2=2000,
           b1=0.1, b2=0.1, I12=10000, I21=10000,
           m=0.9, c1=50, c2=130, C1=50, C2=50,
           bp=8, Kp=300, a=0.004,
           b=1e-1, K=100,
           c=3e-3, v=1e-4, v0=1e-4, cv_v=0.5,
           tmax=1500, S0=95, I0=5,
           minTh2=th2,maxTh2=800,timestep=1)
tIn <- Sys.time()
mclapply(1:50,
         function(x) nested_model_vary_Th2(params),
         mc.cores=12) -> out
tOut <- Sys.time()
print(tOut-tIn)
saveRDS(out, file=paste0("Nested_model_variable_dose_min_Th2=",th2,"_5-11.RDS"))

within_host_model_det = function(t, y, params) {
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
  v <- params["v"]
  v0 <- params["v0"]
  
  T1 <- y[1]
  T2 <- y[2]
  P <- y[3]
  
  dT1 <- b1 + c1*P/(C1+P) + s1*T1^2/(S1^2+T1^2) * I12/(I12+T2) - m*T1
  dT2 <- b2 + c2*P/(C2+P) + s2*T2^2/(S2^2+T2^2) * I21/(I21+T1) - m*T2
  dP <- bp*v/(v0+v)*P*(1-P/Kp) - a*T2*P
  
  list(c(dT1, dT2, dP))
}

## A cool example showing how the high Th1, high Th2, low P deterministic equilibrium disappears with stochasticity
## Very nice extension of previous work!
th2 = 517; p = 40
parms=c(S1=1000, S2=1000, s1=2000, s2=2000,
        b1=0.1, b2=0.1, I12=10000, I21=10000,
        m=0.9, c1=50, c2=130, C1=50, C2=50,
        bp=8, Kp=300, a=0.004, v=1e-4, v0=1e-4)
x0 = c(T1=1200-th2, T2=th2, P=p)
ode(y=x0, 
    times=seq(0,200,0.1), 
    func=within_host_model_det, 
    parms=parms) -> out_det
plot(out_det)

sourceCpp("within-host_model.cpp")
parms=c(S1=1000, S2=1000, s1=2000, s2=2000,
        b1=0.1, b2=0.1, I12=10000, I21=10000,
        m=0.9, c1=50, c2=130, C1=50, C2=50,
        bp=8, Kp=300, a=0.004, v=1e-4, v0=1e-4,
        tmax=200, T1=1200-th2, T2=th2, P=p)
lapply(1:100, function(i) within_host_model(parms)) -> out_stoch

par(mfrow=c(2,2), mar=c(2,4,0,0), oma=c(2,0.5,0.5,0.5))
plot(out_det[,c(1,2)], type='l', ylim=c(0,1800), xlab="", ylab="")
for (i in 1:100) lines(out_stoch[[i]][,c(1,2)], col="gray")
lines(out_det[,c(1,2)], lwd=2)
mtext("Th1", side=2, line=2.5)

plot(out_det[,c(1,3)], type='l', ylim=c(0,1800), xlab="", ylab="")
for (i in 1:100) lines(out_stoch[[i]][,c(1,3)], col="gray")
lines(out_det[,c(1,3)], lwd=2)
mtext("Th2", side=2, line=2.5)

plot(out_det[,c(1,4)], type='l', ylim=c(0,260), xlab="", ylab="")
for (i in 1:100) lines(out_stoch[[i]][,c(1,4)], col="gray")
lines(out_det[,c(1,4)], lwd=2)
mtext("P", side=2, line=2.5)
mtext("Time", side=1, outer=T, line=1)
legend(x='bottomright', legend=paste0("Pr(clearance)=",sum((lapply(out_stoch, function(o) o[201,4]) %>% unlist)==0)/100), text.col='red', bty='n', adj=c(0,-1))

#########################################################################################
#########################################################################################
#########################################################################################

## Figure 1: within-host dynamics at baseline parameter values
th2_seq = seq(200,800)
p_seq = seq(1,40)
init_cond_list = split(expand.grid(Th2=th2_seq, P=p_seq), seq(length(th2_seq)*length(p_seq)))
mclapply(init_cond_list, 
         function(ic) 
           ode(y=c(T1=1200-as.numeric(ic["Th2"]), T2=as.numeric(ic["Th2"]), P=as.numeric(ic["P"])), 
               times=seq(0,200,0.1), 
               func=within_host_model_det, 
               parms=c(S1=1000, S2=1000, s1=2000, s2=2000,
                       b1=0.1, b2=0.1, I12=10000, I21=10000,
                       m=0.9, c1=50, c2=130, C1=50, C2=50,
                       bp=8, Kp=300, a=0.004, b=1e-1, K=100,
                       c=3e-3, v=1e-4, v0=1e-4))[2001,2:4],
         mc.cores=12) -> det_grid
saveRDS(det_grid, file="Deterministic_within-host_outcomes_4-16.RDS")

sourceCpp("within-host_model.cpp")
th2_seq = seq(200,800,by=2)
p_seq = seq(1,40)
init_cond = expand.grid(Th2=th2_seq, P=p_seq)
init_cond_list = lapply(1:nrow(init_cond),
                        function(i)
                          c(S1=1000, S2=1000, s1=2000, s2=2000,
                            b1=0.1, b2=0.1, I12=10000, I21=10000,
                            m=0.9, c1=50, c2=130, C1=50, C2=50,
                            bp=8, Kp=300, a=0.004, v=1e-4, v0=1e-4,
                            tmax=200, 
                            T1=1200-as.numeric(init_cond[i,"Th2"]), 
                            T2=as.numeric(init_cond[i,"Th2"]), 
                            P=as.numeric(init_cond[i,"P"])))
stoch_grid = vector(mode='list', length=length(init_cond_list))
for (i in 1:length(init_cond_list)) {
  print(i)
  mclapply(1:100, 
           function(j) within_host_model(init_cond_list[[i]])[201,2:4],
         mc.cores=14) %>%
    do.call("rbind",.) -> stoch_grid[[i]]
}
saveRDS(stoch_grid, file="Stochastic_within-host_outcomes_4-16.RDS")

det_grid = readRDS(file="Deterministic_within-host_outcomes_4-16.RDS")
th2_seq = seq(200,800)
p_seq = seq(1,40)
det_grid %>% 
  do.call("rbind",.) %>%
  as.data.frame() %>%
  mutate(T0=rep(th2_seq, length(p_seq)),
         P0=rep(p_seq, each=length(th2_seq)),
         outcome = case_when(
           T1 > 800 & T2 > 800  ~ 1,
           T1 > 800 & T2 <= 800 ~ 2,
           T1 <= 800 & T2 > 800 ~ 3,
           T1 <= 800 & T2 <= 800 ~ 4
         )) %>%
  ggplot(., aes(x=T0, y=P0, fill=factor(outcome))) + 
  geom_tile() + 
  annotate("text", x=-Inf, y=Inf, label="A", hjust=-2, vjust=2, color="black", size=4.5) +
  xlab("Initial Th2") + ylab("Initial P") + 
  scale_fill_manual(labels=c("High T1+2, Low P", "High T1, High P", "High T2, Zero P"), values=c("1"="gray", "2"="red", "3"="blue")) + 
  labs(fill="") +
  theme_bw() + 
  theme(legend.position="bottom",
        legend.key.size=unit(0.5,"cm"),
        legend.text=element_text(size=8),
        legend.title=element_text(size=10)) -> p1

stoch_grid = readRDS(file="Stochastic_within-host_outcomes_4-16.RDS")
## No outcomes where you get high Th1, high Th2 and non-zero P
lapply(stoch_grid, function(s) any(s[,1]>800 & s[,2] > 800 & s[,3] > 0)) %>% unlist %>% any
## 19 out of 1,204,000 simulations produce a high Th1, high Th2 outcome - rare enough to ignore
lapply(stoch_grid, function(s) sum(s[,1]>800 & s[,2] > 800)) %>% unlist %>% sum
## Compute extinction probability for each simulation run
th2_seq = seq(200,800,by=2)
p_seq = seq(1,40)
init_cond = expand.grid(Th2=th2_seq, P=p_seq)
mutate(init_cond, 
       p_extinct=lapply(stoch_grid, function(s) sum(s[,3]==0)/100) %>% unlist) -> stoch_outcomes
stoch_outcomes$p_extinct[which(stoch_outcomes$p_extinct>0.999)] = 0.999

stoch_outcomes %>%
  ggplot(., aes(x=Th2, y=P, z=p_extinct)) + 
  geom_contour_filled(aes(fill=after_stat(level_mid))) + 
  scale_fill_gradient(low="red", 
                      high="blue", 
                      limits=c(0,1), 
                      oob=scales::squish,
                      guide=guide_colorbar(
                        barwidth=unit(5,"cm"),
                        ticks=TRUE
                      )) +
  xlab("Initial Th2") + 
  ylab("Initial P") +
  labs(fill="Pr(High T2, Zero P)") +
  annotate("text", x=-Inf, y=Inf, label="B", hjust=-2, vjust=2, color="black", size=4.5) +
  theme_bw() + 
  theme(legend.position="bottom",
        legend.key.size=unit(0.5,"cm"),
        legend.text=element_text(size=8),
        legend.title=element_text(size=10),
        legend.spacing.y=unit(2,"mm"),
        legend.margin=margin(1,1,1,1)) -> p2

png(file="Fig1_Phil_Trans.png", width=5, height=7, units='in', res=400)
p1 / p2
dev.off()




#########################################################################################
#########################################################################################
#########################################################################################



## Range of initial parasite densities in simulations
P0 <- vector(mode='list', length=length(c(seq(200,500,50),seq(525,600,25),seq(650,750,50))))
iii = 1
for (th2 in c(seq(200,500,50),seq(525,600,25),seq(650,750,50))) {
  out2 <- readRDS(file=paste0("Nested_model_variable_dose_min_Th2=",th2,"_4-9.RDS"))
  P0[[iii]] = unlist(lapply(out2[[1]][[1]], function(o) max(o[,7])))
  iii = iii+1
}
max(unlist(P0)) ## 45 is largest dose any individual ever received; the median is 30

## Appendix figure showing the epidmiological dynamics for every min Th2 value
png(file="Epidemiological_dynamics_as_min_Th2_varies_4-9.png", height=9, width=9,  units='in', res=300)
par(mfrow=c(5,3), mar=c(3.5,3.5,0.5,0.5), oma=rep(0,4))
for (th2 in c(seq(200,500,50),seq(525,600,25),seq(650,750,50))) {
  out2 <- readRDS(file=paste0("Nested_model_variable_dose_min_Th2=",th2,"_4-9.RDS"))
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
  legend(x='bottomright', legend=paste0("Prob(fadeout)=",nextinct/50), bty='n', text.col='red')
}
dev.off()


#########################################################################################
#########################################################################################
#########################################################################################

## Figure 2

prev = viru = vector(mode='list', length=length(c(seq(200,500,50),seq(525,600,25),seq(650,750,50))))
iii = 1
for (th2 in c(seq(200,500,50),seq(525,600,25),seq(650,750,50))) {
  out2 <- readRDS(file=paste0("Nested_model_variable_dose_min_Th2=",th2,"_4-9.RDS"))
  prev[[iii]] = data.frame(prev=unlist(lapply(out2, function(o) tail(o[[2]][,3],1)))/100, th2=th2)
  viru[[iii]] = data.frame(prev=unlist(lapply(out2, function(o) tail(o[[2]][,6],1)))/100, th2=th2)
  iii = iii+1
}
prev %<>% do.call("rbind.data.frame",.)
viru %<>% do.call("rbind.data.frame",.)

th2 = 200
out2 <- readRDS(file=paste0("Nested_model_variable_dose_min_Th2=",th2,"_4-9.RDS"))
lapply(1:length(out2), 
       function(i) data.frame(time= out2[[i]][[2]][,1], 
                              infecteds=out2[[i]][[2]][,3]/100, 
                              v=out2[[i]][[2]][,6],
                              color=ifelse(any(out2[[i]][[2]][,3]==0),0,1),
                              ind=i)) %>%
  do.call("rbind.data.frame",.) -> o200

th2 = 575
out2 <- readRDS(file=paste0("Nested_model_variable_dose_min_Th2=",th2,"_4-9.RDS"))
lapply(1:length(out2), 
       function(i) data.frame(time= out2[[i]][[2]][,1], 
                              infecteds=out2[[i]][[2]][,3]/100, 
                              v=out2[[i]][[2]][,6], 
                              color=ifelse(any(out2[[i]][[2]][,3]==0),0,1),
                              ind=i)) %>%
  do.call("rbind.data.frame",.) -> o575

th2 = 700
out2 <- readRDS(file=paste0("Nested_model_variable_dose_min_Th2=",th2,"_4-9.RDS"))
lapply(1:length(out2), 
       function(i) data.frame(time= out2[[i]][[2]][,1], 
                              infecteds=out2[[i]][[2]][,3]/100, 
                              v=out2[[i]][[2]][,6], 
                              color=ifelse(any(out2[[i]][[2]][,3]==0),0,1),
                              ind=i)) %>%
  do.call("rbind.data.frame",.) -> o700


ggplot(o200, aes(x=time, y=infecteds, group=ind, color=as.factor(color))) + 
  geom_line() + 
  scale_color_manual(values=c("0"="red", "1"="gray")) + 
  geom_line(aes(x=time, y=meanI), linewidth=1, color="#0072B2", data=(o200 %>% filter(infecteds > 0) %>% group_by(time) %>% summarize(meanI=mean(infecteds))), inherit.aes=FALSE) + 
  xlab("Time") + ylab("Prevalence") + ylim(0,1) + 
  annotate("text", x=-Inf, y=Inf, label="Min Th2=200", hjust=-0.25, vjust=2, color="#0072B2", size=2.75) +
  annotate("text", x=Inf, y=-Inf, label=paste0("P(fadeout)=",1-sum(filter(o200, time==max(time))$color)/50), hjust=1.25, vjust=-1.5, color="red", size=2.75) +
  theme_bw() + 
  theme(legend.position="none")  -> p1


ggplot(o575, aes(x=time, y=infecteds, group=ind, color=as.factor(color))) + 
  geom_line() + 
  scale_color_manual(values=c("0"="red", "1"="gray")) + 
  geom_line(aes(x=time, y=meanI), linewidth=1, color="#E69F00", data=(o550 %>% filter(infecteds > 0) %>% group_by(time) %>% summarize(meanI=mean(infecteds))), inherit.aes=FALSE) + 
  xlab("Time") + ylab("Prevalence") + ylim(0,1) + 
  annotate("text", x=-Inf, y=Inf, label="Min Th2=575", hjust=-0.25, vjust=2, color="#E69F00", size=2.75) +
  annotate("text", x=Inf, y=-Inf, label=paste0("P(fadeout)=",1-sum(filter(o550, time==max(time))$color)/50), hjust=1.25, vjust=-1.5, color="red", size=2.75) +
  theme_bw() + 
  theme(legend.position="none")  -> p2

ggplot(o700, aes(x=time, y=infecteds, group=ind, color=as.factor(color))) + 
  geom_line() + 
  scale_color_manual(values=c("0"="red", "1"="gray")) + 
  geom_line(aes(x=time, y=meanI), linewidth=1, color="#009E73", data=(o700 %>% filter(infecteds > 0) %>% group_by(time) %>% summarize(meanI=mean(infecteds))), inherit.aes=FALSE) + 
  xlab("Time") + ylab("Prevalence") + ylim(0,1) + 
  annotate("text", x=-Inf, y=Inf, label="Min Th2=700", hjust=-0.25, vjust=2, color="#009E73", size=2.75) +
  annotate("text", x=Inf, y=-Inf, label=paste0("P(fadeout)=",1-sum(filter(o700, time==max(time))$color)/50), hjust=1.25, vjust=-1.5, color="red", size=2.75) +
  theme_bw() + 
  theme(legend.position="none")  -> p3

prev$color = "black"
prev$color[which(prev$th2==200)] = "#0072B2"
prev$color[which(prev$th2==575)] = "#E69F00"
prev$color[which(prev$th2==700)] = "#009E73"

ggplot(prev, aes(x=th2, y=prev, color=color)) + 
  geom_point() + 
  scale_color_manual(values=c("black"="black", "#0072B2"="#0072B2", "#E69F00"="#E69F00", "#009E73"="#009E73")) + 
  theme_bw() + 
  theme(legend.position='none') +
  xlab("Minimum Th2") + 
  ylab("Prevalence at final time") -> p4.1


png(file="Fig2_Epidemiological_dynamics_Phil_Trans.png", width=8, height=6, units='in', res=400)
left = p1 / p2 / p3
left | p4.1
dev.off()

###############################################################################
###############################################################################
###############################################################################

## Figure 3

ggplot(filter(o200, color==1), aes(x=time, y=v, group=ind, color=as.factor(color))) + 
  geom_line() + 
  scale_color_manual(values=c("0"="red", "1"="gray")) + 
  geom_line(aes(x=time, y=meanv), linewidth=1, color="#0072B2", data=(o200 %>% filter(infecteds > 0) %>% group_by(time) %>% summarize(meanv=mean(v))), inherit.aes=FALSE) + 
  xlab("Time") + ylab("Mean virulence") + ylim(0,0.001) + 
  annotate("text", x=-Inf, y=Inf, label="Min Th2=200", hjust=-0.25, vjust=2, color="#0072B2", size=2.75) +
  theme_bw() + 
  theme(legend.position="none")  -> p1


ggplot(filter(o575, color==1), aes(x=time, y=v, group=ind, color=as.factor(color))) + 
  geom_line() + 
  scale_color_manual(values=c("0"="red", "1"="gray")) + 
  geom_line(aes(x=time, y=meanv), linewidth=1, color="#E69F00", data=(o575 %>% filter(infecteds > 0) %>% group_by(time) %>% summarize(meanv=mean(v))), inherit.aes=FALSE) + 
  xlab("Time") + ylab("Mean virulence") + ylim(0,0.001) + 
  annotate("text", x=-Inf, y=Inf, label="Min Th2=575", hjust=-0.25, vjust=2, color="#E69F00", size=2.75) +
  theme_bw() + 
  theme(legend.position="none")  -> p2

ggplot(filter(o700, color==1), aes(x=time, y=v, group=ind, color=as.factor(color))) + 
  geom_line() + 
  scale_color_manual(values=c("0"="red", "1"="gray")) + 
  geom_line(aes(x=time, y=meanv), linewidth=1, color="#009E73", data=(o700 %>% filter(infecteds > 0) %>% group_by(time) %>% summarize(meanv=mean(v))), inherit.aes=FALSE) + 
  xlab("Time") + ylab("Mean virulence") + ylim(0,0.001) + 
  annotate("text", x=-Inf, y=Inf, label="Min Th2=700", hjust=-0.25, vjust=2, color="#009E73", size=2.75) +
  theme_bw() + 
  theme(legend.position="none")  -> p3

viru$color = "black"
viru$color[which(viru$th2==200)] = "#0072B2"
viru$color[which(viru$th2==575)] = "#E69F00"
viru$color[which(viru$th2==700)] = "#009E73"

ggplot(filter(viru, prev > 0), aes(x=th2, y=prev, color=color)) + 
  geom_point() + 
  scale_color_manual(values=c("black"="black", "#0072B2"="#0072B2", "#E69F00"="#E69F00", "#009E73"="#009E73")) + 
  theme_bw() + 
  theme(legend.position='none') +
  xlab("Minimum Th2") + 
  ylab("Virulence at final time") -> p4.1


png(file="Fig3_Evolutionary_dynamics_Phil_Trans.png", width=8, height=6, units='in', res=400)
left = p1 / p2 / p3
left | p4.1
dev.off()

###############################################################################
###############################################################################
###############################################################################

## Figure 4


## Compute information about the duration of each individual infection and the ultimate outcome of infection (peak infection, clearance, death, neither)
for (th2 in c(seq(200,500,50),seq(525,600,25),seq(650,750,50))) {
  print(th2)
  out2 <- readRDS(file=paste0("Nested_model_variable_dose_min_Th2=",th2,"_4-9.RDS"))
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
        # Identify runs of values > 0
        runs <- rle(oo$V3 > 0)
        # Split into runs
        p_groups <- split(oo$V3, rep(seq_along(runs$lengths), runs$lengths))
        t_groups <- split(oo$time, rep(seq_along(runs$lengths), runs$lengths))
        # Take max of positive runs
        max_p <- sapply(p_groups[runs$values], max)
        max_p_t <- sapply(1:length(max_p), function(l) t_groups[runs$values][[l]][which.max(p_groups[runs$values][[l]])]-min(t_groups[runs$values][[l]]))
        times = oo$time[which(oo$V3 > 0)]
        max_age = max(times)
        data.frame(peaks = max_p,
                   peak_time = max_p_t,
                   durations = lengths(split(times, cumsum(c(1, diff(times) != 1)))),
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
  
  saveRDS(duration_data, paste0("Duration_data_min_Th2=",th2,"_4-9.RDS"))
}

library(survival)

infection_durations = vector(mode='list', length=length(c(seq(200,500,50),seq(525,600,25),seq(650,750,50))))
i = 1
for (th2 in c(seq(200,500,50),seq(525,600,25),seq(650,750,50))) {
  duration_data = readRDS(paste0("Duration_data_min_Th2=",th2,"_4-9.RDS"))
  mean_duration = c()
  se_duration = c()
  for (t in seq(0, 450, 50)) {
    filter(duration_data, starts >= t, starts < t+50) %>% 
      mutate(event_code = ifelse(event!="censored",1,0)) %>% 
      with(., survfit(Surv(time = durations, event = event_code)~1)) -> km_fit
    mean_duration = c(mean_duration, summary(km_fit)$table["rmean"])
    se_duration = c(se_duration, summary(km_fit)$table["se(rmean)"])
  }
  infection_durations[[i]] = data.frame(tint=seq(0,450,50),
                                        mean_duration=mean_duration,
                                        se_duration=se_duration,
                                        minTh2=th2)
  i = i + 1
}
infection_durations %<>% do.call("rbind",.)

infection_peaks = vector(mode='list', length=length(c(seq(200,500,50),seq(525,600,25),seq(650,750,50))))
i = 1
for (th2 in c(seq(200,500,50),seq(525,600,25),seq(650,750,50))) {
  duration_data = readRDS(paste0("Duration_data_min_Th2=",th2,"_4-9.RDS"))
  mean_peak = c()
  se_peak = c()
  for (t in seq(0, 450, 50)) {
    filter(duration_data, starts >= t, starts < t+50) -> dat
    mean_peak = c(mean_peak, mean(dat$peaks))
    se_peak = c(se_peak, sd(dat$peaks)/sqrt(nrow(dat)))
  }
  infection_peaks[[i]] = data.frame(tint=seq(0,450,50),
                                        mean_peak=mean_peak,
                                        se_peak=se_peak,
                                        minTh2=th2)
  i = i + 1
}
infection_peaks %<>% do.call("rbind",.)

infection_clearance = vector(mode='list', length=length(c(seq(200,500,50),seq(525,600,25),seq(650,750,50))))
i = 1
for (th2 in c(seq(200,500,50),seq(525,600,25),seq(650,750,50))) {
  duration_data <- readRDS(file=paste0("Duration_data_min_Th2=",th2,"_4-9.RDS"))
  prob_clear = se_clear = c()
  for (t in seq(0, 450, 50)) {
    filter(duration_data, starts >= t, starts < t+50) %>% 
      mutate(event_code=ifelse(event=="clearance",1,ifelse(event=="death",2,0))) %>%
      with(., cuminc(ftime=durations, fstatus=event_code)) -> ci
    prob_clear = c(prob_clear, tail(ci[[1]]$est,1))
    se_clear = c(se_clear, sqrt(tail(ci[[1]]$var,1))/sqrt(nrow(filter(duration_data, starts >= t, starts < t+50))))
  }
  infection_clearance[[i]] = data.frame(tint=seq(0,450,50),
                                        mean_prob=prob_clear,
                                        se_prob=se_clear,
                                        minTh2=th2)
  i = i + 1
}
infection_clearance %<>% do.call("rbind",.)


merge(merge(infection_clearance, infection_durations), infection_peaks) %>% 
  pivot_longer(., cols=c(3,5,7), names_to="Component", values_to="Values") %>%
  filter(., minTh2%in%c(200,300,400,500,550,575,600,650,700,750)) -> dat
dat$Component = factor(dat$Component, levels=c("mean_prob", "mean_duration", "mean_peak"))

dat %>%
  ggplot(., aes(x=tint, y=Values, group=minTh2, color=factor(minTh2))) + 
  geom_line() + 
  scale_color_manual(values = colorRampPalette(c("red", "blue"))(10)) +
  facet_wrap(.~Component, scales="free", labeller=labeller(Component = c("mean_duration"="B. Infection duration", "mean_peak"="C. Peak burden", "mean_prob"="A. Clearance probability"))) + 
  theme_bw() + 
  xlab("Infection start time") + 
  ylab("Fitness component") + 
  labs(color=NULL) + 
  theme(legend.position="bottom") -> p_base
leg <- get_legend(p_base)

legend_row <- plot_grid(
  ggdraw() + draw_label("Th1-biased", hjust = -1, x = 0, size = 10, color="red"),
  leg,
  ggdraw() + draw_label("Th2-biased", hjust = 2, x = 1, size = 10, color="blue"),
  nrow = 1,
  rel_widths = c(0.24, 0.52, 0.24)
)

final_plot <- plot_grid(
  p_base + theme(legend.position = "none"),
  legend_row,
  ncol = 1,
  rel_heights = c(1, 0.2)
)



png(file="Fig4_Fitness_components_Phil_Trans.png", width=6, height=4, units='in', res=400)
final_plot
dev.off()


#########################################################################

out2 <- readRDS(file="Nested_model_variable_dose_min_Th2=575_5-11.RDS")
lapply(1:length(out2), 
       function(i) data.frame(time= out2[[i]][[2]][,1], 
                              infecteds=out2[[i]][[2]][,3]/100, 
                              v=out2[[i]][[2]][,6], 
                              color=ifelse(any(out2[[i]][[2]][,3]==0),0,1),
                              ind=i)) %>%
  do.call("rbind.data.frame",.) -> o575

ggplot(o575, aes(x=time, y=infecteds, group=ind, color=as.factor(color))) + 
  geom_line() + 
  scale_color_manual(values=c("0"="red", "1"="gray")) + 
  geom_line(aes(x=time, y=meanI), linewidth=1, color="#0072B2", data=(o575 %>% filter(infecteds > 0) %>% group_by(time) %>% summarize(meanI=mean(infecteds))), inherit.aes=FALSE) + 
  xlab("Time") + ylab("Prevalence") + ylim(0,1) + 
  annotate("text", x=-Inf, y=Inf, label="Min Th2=575", hjust=-0.25, vjust=2, color="#0072B2", size=2.75) +
  annotate("text", x=Inf, y=-Inf, label=paste0("P(fadeout)=",1-sum(filter(o575, time==max(time))$color)/50), hjust=1.25, vjust=-1.5, color="red", size=2.75) +
  theme_bw() + 
  theme(legend.position="none")  -> p1

ggplot(filter(o575, color==1), aes(x=time, y=v, group=ind, color=as.factor(color))) + 
  geom_line() + 
  scale_color_manual(values=c("0"="red", "1"="gray")) + 
  geom_line(aes(x=time, y=meanv), linewidth=1, color="#E69F00", data=(o575 %>% filter(infecteds > 0) %>% group_by(time) %>% summarize(meanv=mean(v))), inherit.aes=FALSE) + 
  xlab("Time") + ylab("Mean virulence") + ylim(0,0.001) + 
  annotate("text", x=-Inf, y=Inf, label="Min Th2=575", hjust=-0.25, vjust=2, color="#E69F00", size=2.75) +
  theme_bw() + 
  theme(legend.position="none")  -> p2

png(file="FigS1_Extended_timestep_runs.png", width=5, height=3, units='in', res=400)
p1 + p2
dev.off()

###############################################################################
###############################################################################
###############################################################################

## Figure 5

ESS5 = read.csv("Fig5_ESS_5.csv")
ESS10 = read.csv("Fig5_ESS_10.csv")
ESS20 = read.csv("Fig5_ESS_20.csv")
ESSb = read.csv("Fig5_ESS_baseline.csv")

png(filename="Fig5_ESS_result_Phil_Trans.png", height=5, width=6, units='in', res=400)
par(mfrow=c(1,2), mar=c(4,2,0.5,0.5), oma=c(0,2,0,0))
plot(ESS5, xlab=expression("Clearance probability"~p), ylab="", type='l', lwd=2, xlim=c(0.6,1), col="#0072B2")
lines(ESS10, lwd=2, col= "#009E73")
lines(ESS20, lwd=2, col="#E69F00")
mtext("ESS virulence", side=3, line=0, outer=T)
legend(x=0.6, y=0.01, xjust=0.15, yjust=0.5, fill=c("#0072B2","#009E73","#E69F00"), legend=c(expression(gamma==0.2),expression(gamma==0.1),expression(gamma==0.05)), bty='n')
legend(x=0.6, y=0.044, xjust=0.5, yjust=0.5, legend="A", bty="n")

plot(ESSb, xlab=expression("Clearance rate"~gamma), ylab="", type='l', lwd=2, col=1)
legend(x=0, y=0.044, legend="B", xjust=0.5, yjust=0.5, bty='n')
dev.off()
