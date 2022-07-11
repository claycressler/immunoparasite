gillespie_sim <- function(tmax, y0, params, fixedDose=FALSE, seed=NULL) {
  if (!is.null(seed))
    set.seed(seed)
  
  ## Between-host parameters
  c <- params["c"] ## contact rate between hosts
  v <- params["v"] ## virulence (per-capita mortality rate for infected hosts - I assume that parasite load does not determine virulence, but that virulence does determine parasite replication rate)
  v0 <- params["v0"] ## half-saturation constant scaling virulence into parasite replication rate
  
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
  bp <- params["bp"] * v/(v0+v) ## parasite maximum birth rate parameter (depends on virulence)
  Kp <- params["Kp"] ## parasite carrying capacity
  a <- params["a"] ## immune killing rate
  dose <- params["dose"] ## initial dose (if dose is fixed)
  
  ## Every individual in the host population has values of Th1, Th2, and P.
  ## These values determine when/if the host recovers from its infection.
  ## They also determine the immune state of the host when it gets infected.
  S <- y0[1]
  I <- y0[2]
  Host <- vector(mode='list', length=S+I)
  for(i in 1:S) Host[[i]] <- c(0,0,0) ## susceptible host initial state 
  for(i in (S+1):(S+I)) Host[[i]] <- c(0,0,dose) ## infected host initial state
 
  ## time 
  t <- 0
  
  ## set up storage for susceptible hosts, infected hosts, and mean and variance in virulence
  out <- array(0, dim=c(1e6, 5))
  colnames(out) <- c("Time", "S", "I", "meanV", "varV")
  out[1,] <- c(t, S, I, v, 0)
  i <- 2
  
  while (t < tmax) {

    ## extract the current T1, T2, and P values for the host population
    T1 <- lapply(Host, function(i) i[1]) %>% unlist
    T2 <- lapply(Host, function(i) i[2]) %>% unlist
    P <- lapply(Host, function(i) i[3]) %>% unlist
    
    ## compute within-host event rates
    ## to prevent the system from getting "stuck", self-stimulation/cross-inhibition turns "off" when the parasite goes extinct, reflecting the action of Treg cells
    prod1 <- sapply(P, function(p) ifelse(p > 0, b1 + c1*p/(C1+p) + s1*T1^2/(S1^2+T1^2) * I12/(I12+T2), b1))
    prod2 <- sapply(P, function(p) ifelse(p > 0, b1 + c1*p/(C1+p) + s1*T1^2/(S1^2+T1^2) * I12/(I12+T2), b1))
    death1 <- m*T1
    death2 <- m*T2
    birthP <- bp*P*(1-P/Kp)
    deathP <- a*T2*P
    contact <- c*S*I
    deathI <- v*I
    rates <- c(prod1, prod2, death1, death2, birthP, deathP, contact, deathI)
    
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
      Host[[ind]][1] <- Host[[event]][1]+1
    }
    ## If event is (S+I+1):(2*(S+I)), a new Th1 cell is produced
    if (event%in%seq(S+I+1,2*(S+I))) {
      ind <- event-(S+I)
      Host[[ind]][2] <- Host[[event]][2]+1
    }
    ## If event is (2*(S+I)+1):(3*(S+I)), a Th1 cell is lost
    if (event%in%seq(2*(S+I)+1,3*(S+I))) {
      ind <- event-2*(S+I)
      Host[[ind]][1] <- Host[[event]][1]-1
    }  
    ## If event is (3*(S+I)+1):(4*(S+I)), a Th2 cell is lost
    if (event%in%seq(3*(S+I)+1,4*(S+I))) {
      ind <- event-3*(S+I)
      Host[[ind]][2] <- Host[[event]][2]-1
    }  
    ## If event is (4*(S+I)+1):(5*(S+I)), a parasite is born
    if (event%in%seq(4*(S+I)+1,5*(S+I))) {
      ind <- event-4*(S+I)
      Host[[ind]][3] <- Host[[event]][3]+1
    }
    ## If event is (5*(S+I)+1):(6*(S+I)), a parasite is killed
    if (event%in%seq(5*(S+I)+1,6*(S+I))) {
      ind <- event-5*(S+I)
      Host[[ind]][3] <- Host[[event]][3]-1
    }
    ## If event is 6*(S+I)+1, there was a contact between a susceptible and an infected host
    if (event==(6*(S+I)+1)) {
      ## choose the susceptible and infected host at random
      Sinds <- which(unlist(lapply(Host, function(h) h[3]))==0)
      Iinds <- which(unlist(lapply(Host, function(h) h[3]))>0)
      whichS <- Sinds[sample(1:length(Sinds), 1)]
      whichI <- Iinds[sample(1:length(Iinds), 1)]
      
      ## infectious dose can be fixed or based on current infection load
      if (fixedDose)
        Host[[whichS]][3] <- dose
      else {
        thisP <- Host[[whichI]][3]/10
        ## assume that the mean and sd are the same (big assumption!)
        mu <- log(thisP^2 / sqrt(thisP^2+thisP^2))
        sigma <- sqrt(log(2))
        Host[[whichS]][3] <- ceiling(rlnorm(1, meanlog=mu, sdlog=sigma)))
      }
    }
    ## If event is 6*(S+I)+2, an infected host dies
    if (event==(6*(S+I)+2)) {
      ## choose the infected host at random
      Iinds <- which(unlist(lapply(Host, function(h) h[3]))>0)
      whichI <- Iinds[sample(1:length(Iinds), 1)]
      Host <- Host[-whichI]
    }
    
     out[i,] <- c(t, T1, T2, P)
    
    i <- i + 1
    #print(out[i,])
    if (i > nrow(out)) ## add more rows
      out <- rbind(out, array(0, dim=c(1e6, 4)))
  }
  out <- out[1:(i-1),]
  return(out)
}

params = c(S1=1000, S2=1000, s1=2000, s2=2000, 
           b1=0, b2=0, I12=10000, I21=10000, 
           m=0.9, c1=50, c2=130, C1=50, C2=50, 
           bp=4,Kp=300, a=0.004, dose=30,
           c=1e-3, v=1e-4, v0=1e-5)

