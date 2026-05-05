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
  scale_fill_manual(labels=c("High T1+2, Low P", "High T1, High P", "High T2, Zero P"), values=c("1"="gray", "2"="red", "3"="black")) + 
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
                      high="black", 
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
  out2 <- readRDS(file=paste0("Nested_model_variable_dose_min_Th2=",th2,"_3-25.RDS"))
  prev[[iii]] = data.frame(prev=unlist(lapply(out2, function(o) tail(o[[2]][,3],1)))/100, th2=th2)
  viru[[iii]] = data.frame(prev=unlist(lapply(out2, function(o) tail(o[[2]][,6],1)))/100, th2=th2)
  iii = iii+1
}
prev %<>% do.call("rbind.data.frame",.)
viru %<>% do.call("rbind.data.frame",.)

th2 = 200
out2 <- readRDS(file=paste0("Nested_model_variable_dose_min_Th2=",th2,"_3-25.RDS"))
lapply(1:length(out2), 
       function(i) data.frame(time= out2[[i]][[2]][,1], 
                              infecteds=out2[[i]][[2]][,3]/100, 
                              v=out2[[i]][[2]][,6],
                              color=ifelse(any(out2[[i]][[2]][,3]==0),0,1),
                              ind=i)) %>%
  do.call("rbind.data.frame",.) -> o200

th2 = 550
out2 <- readRDS(file=paste0("Nested_model_variable_dose_min_Th2=",th2,"_3-25.RDS"))
lapply(1:length(out2), 
       function(i) data.frame(time= out2[[i]][[2]][,1], 
                              infecteds=out2[[i]][[2]][,3]/100, 
                              v=out2[[i]][[2]][,6], 
                              color=ifelse(any(out2[[i]][[2]][,3]==0),0,1),
                              ind=i)) %>%
  do.call("rbind.data.frame",.) -> o550

th2 = 700
out2 <- readRDS(file=paste0("Nested_model_variable_dose_min_Th2=",th2,"_3-25.RDS"))
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


ggplot(o550, aes(x=time, y=infecteds, group=ind, color=as.factor(color))) + 
  geom_line() + 
  scale_color_manual(values=c("0"="red", "1"="gray")) + 
  geom_line(aes(x=time, y=meanI), linewidth=1, color="#E69F00", data=(o550 %>% filter(infecteds > 0) %>% group_by(time) %>% summarize(meanI=mean(infecteds))), inherit.aes=FALSE) + 
  xlab("Time") + ylab("Prevalence") + ylim(0,1) + 
  annotate("text", x=-Inf, y=Inf, label="Min Th2=550", hjust=-0.25, vjust=2, color="#E69F00", size=2.75) +
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
prev$color[which(prev$th2==550)] = "#E69F00"
prev$color[which(prev$th2==700)] = "#009E73"

ggplot(prev, aes(x=th2, y=prev, color=color)) + 
  geom_point() + 
  scale_color_manual(values=c("black"="black", "#0072B2"="#0072B2", "#E69F00"="#E69F00", "#009E73"="#009E73")) + 
  theme_bw() + 
  scale_x_continuous(breaks=c(200, 400, 550, 600, 700, 800)) + 
  theme(axis.text.x = element_text(face=c("bold","plain","bold","plain","bold","plain"),
                                   colour=c("#0072B2","black","#E69F00","black","#009E73","black")),
        legend.position='none') +
  xlab("Minimum Th2") + 
  ylab("Prevalence at final time") -> p4.1


ggplot(prev, aes(x=as.factor(th2), y=prev, color=color)) + 
  geom_boxplot() + 
  scale_color_manual(values=c("black"="black", "#0072B2"="#0072B2", "#E69F00"="#E69F00", "#009E73"="#009E73")) + 
  theme_bw() + 
  theme(axis.text.x = element_text(face=c("bold",rep("plain",7),"bold",rep("plain",3),"bold","plain"),
                                   colour=c("#0072B2",rep("black",7),"#E69F00",rep("black",3),"#009E73","black")),
        legend.position='none') +
  xlab("Minimum Th2ness") + 
  ylab("Evo-eco equilibrium prevalence") -> p4.2

png(file="Epidemiological_dynamics_summary_fig.png", width=8, height=6, units='in', res=400)
left = p1 / p2 / p3
left | p4.1
dev.off()

png(file="Epidemiological_dynamics_summary_fig_2.png", width=8, height=6, units='in', res=400)
left = p1 / p2 / p3
left | p4.2
dev.off()
