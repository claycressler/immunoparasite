library(pryr)
library(tidyverse)

workspace.size <- function() {
  ws <- sum(sapply(ls(envir=globalenv()), function(x)object_size(get(x))))
  class(ws) <- "object_size"
  ws
}


nested_modelR <- function(params, seed=NULL, fixedDose=FALSE) {
  if (!is.null(seed))
    set.seed(seed)
  
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
  
  ## simulation parameters
  S <- params["S0"]
  I <- params["I0"]
  R <- 0 ## number of recoveries
  tmax <- params["tmax"]
  
  ## Every individual in the host population has values of Th1, Th2, P, and v (if infected [P>0]).
  ## These values determine when/if the host recovers from its infection.
  ## They also determine the immune state of the host when it gets infected.
  Host <- vector(mode='list', length=S+I)
  for(i in 1:S) Host[[i]] <- c(0,0,0,0,i) ## susceptible host initial state 
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
    Host[[i]] <- c(0,0,this_dose,this_v,i) ## infected host initial state
  }
  
  ## time 
  t <- 0
  
  ## set up storage for susceptible hosts, infected hosts, no. recoveries, and mean virulence
  ## store every 0.1
  Popn <- array(0, dim=c(length(seq(0,tmax,0.1))-1, 5))
  HostState <- vector(mode='list', length=nrow(Popn))
  i <- 1
  while (t < tmax) {
    
    ## extract the current T1, T2, P, and v values for the host population
    T1 <- lapply(Host, function(i) i[1]) %>% unlist
    T2 <- lapply(Host, function(i) i[2]) %>% unlist
    P <- lapply(Host, function(i) i[3]) %>% unlist
    v <- lapply(Host, function(i) i[4]) %>% unlist
    
    ## extract the current S and I state for the population (R is tracked below)
    S <- sum(P==0)
    I <- sum(P>0)
    
    ## if it's time, store the current system state (t, S, I, R mean virulence) 
    if(t >= seq(0,tmax,0.1)[i]) {
      ## check for memory "leaks" (dramatic increases in memory usage over runtime)
      ##print(workspace.size())
      Popn[i,] <- c(t, S, I, R, mean(v[v>0]))
      HostState[[i]] <- do.call("rbind.data.frame", Host)
      colnames(HostState[[i]]) <- c("T1","T2","P","v")
      i <- i+1
    }
    
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
    
    ## what time does the event happen?
    dt <- rexp(1, rate=sum(rates))
    
    ## update t 
    t <- t + dt
    
    ## "wheel of fortune"
    wheel <- cumsum(rates)/sum(rates)
    
    ## which event happens? Draw a random uniform to determine
    rand <- runif(1)
    
    event <- 1 + sum(rand > wheel)
    
    ## If event is 1:(S+I), a new Th1 cell is produced
    if (event%in%seq(1,S+I)) {
      ind <- event
      Host[[ind]][1] <- Host[[ind]][1]+1
    }
    ## If event is (S+I+1):(2*(S+I)), a new Th1 cell is produced
    if (event%in%seq(S+I+1,2*(S+I))) {
      ind <- event-(S+I)
      Host[[ind]][2] <- Host[[ind]][2]+1
    }
    ## If event is (2*(S+I)+1):(3*(S+I)), a Th1 cell is lost
    if (event%in%seq(2*(S+I)+1,3*(S+I))) {
      ind <- event-2*(S+I)
      Host[[ind]][1] <- Host[[ind]][1]-1
    }  
    ## If event is (3*(S+I)+1):(4*(S+I)), a Th2 cell is lost
    if (event%in%seq(3*(S+I)+1,4*(S+I))) {
      ind <- event-3*(S+I)
      Host[[ind]][2] <- Host[[ind]][2]-1
    }  
    ## If event is (4*(S+I)+1):(5*(S+I)), a parasite is born
    if (event%in%seq(4*(S+I)+1,5*(S+I))) {
      ind <- event-4*(S+I)
      Host[[ind]][3] <- Host[[ind]][3]+1
    }
    ## If event is (5*(S+I)+1):(6*(S+I)), a parasite is killed
    if (event%in%seq(5*(S+I)+1,6*(S+I))) {
      ind <- event-5*(S+I)
      Host[[ind]][3] <- Host[[ind]][3]-1
      ## if P = 0 now, set v = 0 so host can no longer die of infection and increment the no. of recoveries
      if(Host[[ind]][3]==0) {
        Host[[ind]][4] <- 0
        R <- R + 1
      }
    }
    ## If event is (6*(S+I)+1):(7*(S+I)), a host dies of infection
    if (event%in%seq(6*(S+I)+1,7*(S+I))) {
      ind <- event-6*(S+I)
      Host <- Host[-ind]
    }
    ## If event is 7*(S+I)+1, there was a contact between a susceptible and an infected host
    if (event==(7*(S+I)+1)) {
      ## choose the susceptible and infected host at random
      Sinds <- which(unlist(lapply(Host, function(h) h[3]))==0)
      Iinds <- which(unlist(lapply(Host, function(h) h[3]))>0)
      whichS <- Sinds[sample(1:length(Sinds), 1)]
      whichI <- Iinds[sample(1:length(Iinds), 1)]
      ## What is the dose and virulence of the new infection?
      ## infectious dose can be fixed or based on current infection load
      if (fixedDose)
        Host[[whichS]][3] <- Kp/10
      else ## dose is Poisson distributed based on the current infection load of the infecting individual
        Host[[whichS]][3] <- rpois(1, lambda=Host[[whichI]][3]/10)
      ## virulence also depends on the virulence of the infecting individual
      this_v <- Host[[whichI]][4]
      sd_v <- cv_v * this_v
      mu <- log(this_v^2 / sqrt(sd_v^2+this_v^2))
      sigma <- sqrt(log(cv_v^2/this_v^2 + 1))
      Host[[whichS]][4] <- rlnorm(1, meanlog=mu, sdlog=sigma)
    }
  }
  return(list(HostState,Popn))
}