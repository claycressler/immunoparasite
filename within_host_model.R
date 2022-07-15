library(deSolve)

within_host_model = function(t, y, params) {
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
