library(Rcpp)
sourceCpp("nested_model.cpp")

params = c(S1=1000, S2=1000, s1=2000, s2=2000, 
           b1=0, b2=0, I12=10000, I21=10000, 
           m=0.9, c1=50, c2=130, C1=50, C2=50, 
           bp=4,Kp=300, a=0.004, 
           c=1e-3, v=1e-2, v0=1e-3, cv_v=0.5,
           tmax=20, S0=10, I0=5)

system.time(nested_model(params) -> out)

system.time(nested_modelR(params) -> out)

nested_modelR <- function(params) {
  ## Between-host parameters
  c <- params["c"] ## contact rate between hosts
  v <- params["v"] ## virulence (per-capita mortality rate for infected hosts - I assume that parasite load does not determine virulence, but that virulence does determine parasite replication rate)
  v0 <- params["v0"] ## half-saturation constant scaling virulence into parasite replication rate
  cv_v <- params["cv_v"] ## coefficient of variation in virulence (for evolving simulations)
  
  ## Within-host parameters
  s1 <- params["s1"] ## self-stimulation of Th1
  s2 <- params["s2"] ## self-stimulation of Th2 
  b1 <- params["b1"] ## baseline production of Th1
  b2 <- params["b2"] ## baseline production of Th2
  I12 <- params["I12"] ## inhibition of Th1 by Th2
  I21 <- params["I21"] ## inhibition of Th2 by Th1
  S1 <- params["S1"] ## Th1 self-stimulation half-saturation constant
  S2 <- params["S2"] ## Th2 self-stimulation half-saturation constant 
  m <- params["m"] ## cytokine decay/downregulation
  c1 <- params["c1"] ## Th1 production due to parasites
  c2 <- params["c2"] ## Th2 production due to parasites
  C1 <- params["C1"] ## Th1 production half-saturation constant
  C2 <- params["C2"] ## Th2 production half-saturation constant
  bp <- params["bp"] ## parasite maximum birth rate parameter
  Kp <- params["Kp"] ## parasite carrying capacity
  a <- params["a"] ## immune killing rate
  
  ## Every individual in the host population has values of Th1, Th2, P, and v (if infected [P>0]).
  ## These values determine when/if the host recovers from its infection.
  ## They also determine the immune state of the host when it gets infected.
  S <- params["S0"]
  I <- params["I0"]
  R <- 0 ## number of recoveries
  Host <- vector(mode='list', length=S+I)
  for(i in 1:S) Host[[i]] <- c(0,0,0,0) ## susceptible host initial state 
  for(i in (S+1):(S+I)) {
    ## set virulences of initially infected hosts
    sd_v <- cv_v * v
    mu <- log(v^2 / sqrt(sd_v^2+v^2))
    sigma <- sqrt(log(cv_v^2/v^2 + 1))
    this_v <- rlnorm(1, meanlog=mu, sdlog=sigma)
    ## set initial infected doses
    if (fixedDose) ## if dose is fixed, it is 10% of the carrying capacity
      this_dose <- Kp/10
    else ## assume that dose is Poisson distributed (so variance=mean)
      this_dose <- rpois(1, lambda=Kp/10)
    Host[[i]] <- c(0,0,this_dose,this_v) ## infected host initial state
  }
  
  ## time 
  t <- 0
  
  ## set up storage for susceptible hosts, infected hosts, no. recoveries, and mean virulence
  ## store every 0.1
  out <- array(0, dim=c(length(seq(0,tmax,0.1))-1, 5))
  i <- 1
   ##while (t < tmax) {
  for (j in 1:100000) {
    
    ## extract the current T1, T2, P, and v values for the host population
    T1 <- lapply(Host, function(i) i[1]) %>% unlist
    T2 <- lapply(Host, function(i) i[2]) %>% unlist
    P <- lapply(Host, function(i) i[3]) %>% unlist
    v <- lapply(Host, function(i) i[4]) %>% unlist
    
    ## extract the current S and I state for the population (R is tracked below)
    S <- sum(P==0)
    I <- sum(P>0)
    
    ## compute within-host event rates
    ## to allow the immune system to "reset" after clearance, self-stimulation/cross-inhibition turns "off" when the parasite goes extinct, reflecting the action of Treg cells
    prod1 <- ifelse(P > 0, b1 + c1*P/(C1+P) + s1*T1^2/(S1^2+T1^2) * I12/(I12+T2), b1)
    prod2 <- ifelse(P > 0, b2 + c2*P/(C2+P) + s2*T2^2/(S2^2+T2^2) * I21/(I21+T1), b2)
    death1 <- m*T1
    death2 <- m*T2
    birthP <- bp*v/(v0+v)*P*(1-P/Kp)
    deathP <- a*T2*P
    deathI <- v
    contact <- c*S*I
    rates <- c(prod1, prod2, death1, death2, birthP, deathP, deathI, contact)
  }
  return(rates)
}
