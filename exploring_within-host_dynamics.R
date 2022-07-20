source("within_host_model.R")
out <- lapply(seq(10,100,10), function(p) ode(y=c(T1=800,T2=500,P=p), times=seq(0,50,0.1), within_host_model, params))
## Infections are always chronic - the system evolves to a lower bp value because it makes all infections chronic
lapply(out, function(o) o[501,])

runSilent <- function(BP, t1, t2, p) {
  options(warn = -1)
  on.exit(options(warn = 0))
  params = c(S1=1000, S2=1000, s1=2000, s2=2000, b1=0.1, b2=0.1, I12=10000, I21=10000, m=0.9, c1=50, c2=130, C1=50, C2=50, bp=BP, Kp=300, a=0.004)
  capture.output(res <- try(lsoda(y=c(T1=t1,T2=t2,P=p), times=0:100,
                                  within_host_model,params)))
  if (!inherits(res, "try-error") && nrow(res)==101)
    res2 <- res[101,4]
  else res2 <- "error"
  res2
}

bpRates <- seq(1,8,0.2) %>% as.list
mclapply(bpRates, 
         function(BP) {
           params = c(S1=1000, S2=1000, s1=2000, s2=2000, b1=0.1, b2=0.1, I12=10000, I21=10000, m=0.9, c1=50, c2=130, C1=50, C2=50, bp=BP, Kp=300, a=0.004)
           results <- expand.grid(dose=seq(10,100,5), totalT=seq(200,2000,100), T2=0)
           for (i in 1:nrow(results)) {
             #print(i)
             litT2 <- 0
             bigT2 <- results$totalT[i]
             midT2 <- floor(bigT2/2)
             P1 = runSilent(BP=BP, t1=results$totalT[i]-litT2, t2=litT2, p=results$dose[i])
             P2 = runSilent(BP=BP, t1=results$totalT[i]-midT2, t2=midT2, p=results$dose[i])
             P3 = runSilent(BP=BP, t1=results$totalT[i]-bigT2, t2=bigT2, p=results$dose[i])
             if(sum(round(c(P1,P2,P3)))==0) 
               results[i,"T2"] <- 0
             else if(all(c(P1,P2,P3) > 20))
               results[i,"T2"] <- results$totalT[i]
             else {
               while(bigT2-midT2 > 1 && midT2-litT2 > 1 && P2!="error") {
                 if(round(abs(P1-P2)) > round(abs(P2-P3))) { ## acute to chronic switch happens in the first interval
                   ## litT2 stays the same and bigT2 becomes midT2 =>
                   ## P1 remains the same and P3 becomes P2
                   bigT2 <- midT2
                   P3 <- P2
                   ## midT2 is halfway between litT2 and bigT2
                   midT2 <- floor((bigT2+litT2)/2)
                   P2 = runSilent(BP=BP, t1=results$totalT[i]-midT2, t2=midT2, p=results$dose[i])
                   if (P2=="error") {
                     Perr <- lapply(seq(litT2,bigT2), function(midT2) runSilent(BP=BP, t1=results$totalT[i]-midT2, t2=midT2, p=results$dose[i]))
                     midT2 <- seq(litT2,bigT2)[(lapply(Perr, function(p) p < 1) %>% unlist %>% which %>% min)]
                   }
                 } else {
                   ## bigT2 stays the same and litT2 becomes midT2 =>
                   ## P3 remains the same and P1 becomes P2
                   litT2 <- midT2
                   P1 <- P2
                   ## midT2 is halfway between litT2 and bigT2
                   midT2 <- floor((bigT2+litT2)/2)
                   P2 = runSilent(BP=BP, t1=results$totalT[i]-midT2, t2=midT2, p=results$dose[i])
                   if (P2=="error") {
                     Perr <- lapply(seq(litT2,bigT2), function(midT2) runSilent(BP=BP, t1=results$totalT[i]-midT2, t2=midT2, p=results$dose[i]))
                     midT2 <- seq(litT2,bigT2)[(lapply(Perr, function(p) p < 1) %>% unlist %>% which %>% min)]
                   }
                 }
               }
               results[i,"T2"] <- midT2
             }
           }
           return(results)
         },
         mc.cores=12) -> output
saveRDS(output, file="Th2_thresholds_between_acute_and_chronic_across_doses_totalT-cells_and_parasite_replication.RDS")

## Can I find a set of initial conditions where the outcome is always variable, regardless of the bp and the dose value?
## Also need one such that the outcome shifts at the same value of T2, or at least in similar range of T2 values, and I may need to consider incorporating some variation in initial T2ness

## There is no totalT value where there is variation across every value of dose and across every value of bp
sapply(seq(200,2000,100), function(tot) all(unlist(lapply(output, function(o) all(with(subset(o, totalT==tot), totalT-T2)>0,na.rm=T)))))

## that was probably a big ask, honestly
## For any totalT > 700, there are always 31/36 bp values that produces variable outcomes
sapply(seq(200,2000,100), function(tot) sum(unlist(lapply(output, function(o) all(with(subset(o, totalT==tot), totalT-T2)>0,na.rm=T)))))

## And it's always for the lower values of bp
lapply(output, function(o) with(subset(o, totalT==800), totalT-T2))
lapply(output, function(o) with(subset(o, totalT==2000), totalT-T2))

## what is happening at large totalT?
## extinction for every dose
sapply(seq(10,60,2), function(d) runSilent(BP=2, t1=2000-1000, t2=1000, p=d))

## chronic for every dose
sapply(seq(10,60,2), function(d) runSilent(BP=2, t1=2000-500, t2=500, p=d))

sapply(seq(100,1900,100), function(t) all(sapply(seq(10,60,2), function(d) runSilent(BP=2, t1=2000-t, t2=t, p=d)) < 1))

## so there may be a switchpoint in here, but it's hard to find
sapply(seq(10,60,2), function(d) runSilent(BP=2, t1=2000-900, t2=900, p=d))

## But if you look at these lower values of bp, where we see lots of variation 
lapply(output, function(o) with(subset(o, totalT==1200), T2))
## as bp gets smaller, the Th2 breakpoints get smaller
runSilent(BP=1, t1=1200-664, t2=664, p=10)
runSilent(BP=1, t1=1200-665, t2=665, p=10)
## to this...
runSilent(BP=5, t1=1200-477, t2=477, p=10)
runSilent(BP=5, t1=1200-478, t2=478, p=10)

## smaller doses for the same bp and Th2ness are more likely to produce an acute infection 
runSilent(BP=1, t1=1200-664, t2=664, p=10)
runSilent(BP=1, t1=1200-664, t2=664, p=12)

runSilent(BP=5, t1=1200-477, t2=477, p=10)
runSilent(BP=5, t1=1200-477, t2=477, p=12)

## at large enough bp values, infections are ALWAYS chronic, regardless of dose or initial Th2ness
sapply(1:10, function(d) runSilent(BP=8, t1=1200-1199, t2=1199, p=d))

## That suggests the potential that an initial bp that was high enough could cause virulence to evolve upwards rather than downwards (evolutionary bistability)


## what that means is that we need to allow for variation in initial Th2ness and initial dose
range(unlist(lapply(output[1:31], function(o) with(subset(o, totalT==1200), T2))))
## draw initial Th2ness from a uniform distribution with min=400 and max=700 and fix totalT=1200

## For this fixed Th2ness, low doses give you chronic and high doses give you acute
sapply(1:100, function(d) runSilent(BP=4, t1=1200-500, t2=500, p=d))
## so reducing bp evolutionarily can be advantageous because it makes the likely initial dose smaller

## interestingly, there is a weird trough, though. 
## Notice: 
## for bp = 1,2,3, infections are chronic regardless of the initial dose
params = c(S1=1000, S2=1000, s1=2000, s2=2000, b1=0.1, b2=0.1, I12=10000, I21=10000, m=0.9, c1=50, c2=130, C1=50, C2=50, bp=1, Kp=300, a=0.004)
sapply(1:50, function(d) lsoda(y=c(T1=700,T2=500,P=d), times=0:200, within_host_model,params)[201,4])
## for bp = 4,5 infections are chronic for low dose and acute for high dose
params = c(S1=1000, S2=1000, s1=2000, s2=2000, b1=0.1, b2=0.1, I12=10000, I21=10000, m=0.9, c1=50, c2=130, C1=50, C2=50, bp=4, Kp=300, a=0.004)
sapply(1:50, function(d) lsoda(y=c(T1=700,T2=500,P=d), times=0:200, within_host_model,params)[201,4])
params = c(S1=1000, S2=1000, s1=2000, s2=2000, b1=0.1, b2=0.1, I12=10000, I21=10000, m=0.9, c1=50, c2=130, C1=50, C2=50, bp=5, Kp=300, a=0.004)
sapply(1:50, function(d) lsoda(y=c(T1=700,T2=500,P=d), times=0:200, within_host_model,params)[201,4])
## for bp = 6 infections are acute for all doses
params = c(S1=1000, S2=1000, s1=2000, s2=2000, b1=0.1, b2=0.1, I12=10000, I21=10000, m=0.9, c1=50, c2=130, C1=50, C2=50, bp=7, Kp=300, a=0.004)
sapply(1:50, function(d) lsoda(y=c(T1=700,T2=500,P=d), times=0:200, within_host_model,params)[201,4])
## for bp = 7,8 infections are chronic for all doses
params = c(S1=1000, S2=1000, s1=2000, s2=2000, b1=0.1, b2=0.1, I12=10000, I21=10000, m=0.9, c1=50, c2=130, C1=50, C2=50, bp=8, Kp=300, a=0.004)
sapply(1:50, function(d) lsoda(y=c(T1=700,T2=500,P=d), times=0:200, within_host_model,params)[201,4])

## even a moderate initial Th2 bias produces chronic infections only at very high bp
params = c(S1=1000, S2=1000, s1=2000, s2=2000, b1=0.1, b2=0.1, I12=10000, I21=10000, m=0.9, c1=50, c2=130, C1=50, C2=50, bp=7, Kp=300, a=0.004)
sapply(1:50, function(d) lsoda(y=c(T1=500,T2=700,P=d), times=0:200, within_host_model,params)[201,4])

## how long do infections last for different initial doses in this state?
params = c(S1=1000, S2=1000, s1=2000, s2=2000, b1=0.1, b2=0.1, I12=10000, I21=10000, m=0.9, c1=50, c2=130, C1=50, C2=50, bp=4, Kp=300, a=0.004)
(sapply(1:50, function(d) {
  o = lsoda(y=c(T1=500,T2=700,P=d), times=0:200, within_host_model,params)
  c(o[min(which(o[,4]<1)),1], max(o[,4]))}) %>% t) -> try1
## larger doses produce larger maxes, but shorter durations of the infection

## reducing bp actually just makes everything worse for the parasite
params = c(S1=1000, S2=1000, s1=2000, s2=2000, b1=0.1, b2=0.1, I12=10000, I21=10000, m=0.9, c1=50, c2=130, C1=50, C2=50, bp=3, Kp=300, a=0.004)
(sapply(1:50, function(d) {
  o = lsoda(y=c(T1=500,T2=700,P=d), times=0:200, within_host_model,params)
  c(o[min(which(o[,4]<1)),1], max(o[,4]))}) %>% t) -> try2

## increasing it is better. so the issue is whether the population can persist long enough to evolve towards higher bp, as we would expect it eventually would
params = c(S1=1000, S2=1000, s1=2000, s2=2000, b1=0.1, b2=0.1, I12=10000, I21=10000, m=0.9, c1=50, c2=130, C1=50, C2=50, bp=5, Kp=300, a=0.004)
(sapply(1:50, function(d) {
  o = lsoda(y=c(T1=500,T2=700,P=d), times=0:200, within_host_model,params)
  c(o[min(which(o[,4]<1)),1], max(o[,4]))}) %>% t) -> try3
