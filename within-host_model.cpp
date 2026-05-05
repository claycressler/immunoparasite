#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix within_host_model(NumericVector params) {
  //int nested_model(NumericVector params) {
  // Extract all relevant model parameters and algorithm parameters from params
  // Note that it is critical that parameters are specified in EXACTLY the order they are extracted in
  double S1 = params[0]; // Th1 self-stimulation half-saturation constant
  double S2 = params[1]; // Th2 self-stimulation half-saturation constant
  double s1 = params[2]; // self-stimulation of Th1
  double s2 = params[3]; // self-stimulation of Th2
  double b1 = params[4]; // baseline production of Th1
  double b2 = params[5]; // baseline production of Th2
  double I12 = params[6]; // inhibition of Th1 by Th2
  double I21 = params[7]; // inhibition of Th2 by Th1
  double m = params[8]; // cytokine decay/downregulation
  double c1 = params[9]; // Th1 production due to parasites
  double c2 = params[10]; // Th2 production due to parasites
  double C1 = params[11]; // Th1 production half-saturation constant
  double C2 = params[12]; // Th1 production half-saturation constant
  double bp = params[13]; // parasite maximum birth rate parameter
  double Kp = params[14]; // parasite carrying capacity
  double a = params[15]; // immune killing rate
  double v = params[16]; // virulence (per-capita mortality rate for infected hosts - I assume that parasite load does not determine virulence, but that virulence does determine parasite replication rate)
  double v0 = params[17]; // half-saturation constant scaling virulence into parasite replication rate
  double tmax = params[18]; // length of time to run simulation
  int T1 = params[19];
  int T2 = params[20];
  int P = params[21];
  
  // initialize time
  double t = 0.0; 
  
  // initialize storage for population statistics 
  NumericMatrix Popn(tmax + 1, 4);
  NumericVector ones = rep(1.0,Popn.nrow());
  ones(0) = 0.0;
  NumericVector times(ones.length());
  std::partial_sum(ones.begin(), ones.end(), times.begin());
  Popn(_,0) = times;
  int tIter = 0;
  
  // create storage for the events so they aren't created every time through the loop
  double prod1,prod2,decay1,decay2,birthP,deathP, totalRate, rand;
  int event;
  IntegerVector whichWheel,whichShosts,whichIhosts; 
  
  while(t < tmax & P > 0) {
    // store population information?
    if (t >= times(tIter)) {
      Popn(tIter, 1) = T1;
      Popn(tIter, 2) = T2;
      Popn(tIter, 3) = P;
      tIter += 1;
    }
    
    // compute the rates of all model processes
    // speeds up computation considerably to assume no production or decay of Th1 and Th2 in uninfected hosts
    // biologically this is unreasonable but makes sense in the context of setting the initial Th1 and Th2 of
    // every newly infected host
    // basically just means we have less extraneous events occurring so time will advance faster
    prod1 = b1 + c1*P/(C1+P) + s1*pow(T1,2.0)/(pow(S1,2.0)+pow(T1,2.0)) * I12/(I12+T2);
    prod2 = b2 + c2*P/(C2+P) + s2*pow(T2,2.0)/(pow(S2,2.0)+pow(T2,2.0)) * I21/(I21+T1);
    decay1 = m*T1;
    decay2 = m*T2;
    birthP = bp*v/(v0+v)*P*(1-P/Kp);
    deathP = a*T2*P;
    
    // combine all rates into a single vector
    NumericVector rates(6);
    rates(0) = prod1;
    rates(1) = prod2;
    rates(2) = decay1;
    rates(3) = decay2;
    rates(4) = birthP;
    rates(5) = deathP;
    // Compute the total rate of events
    totalRate = std::accumulate(rates.begin(), rates.end(), 0.0);
    // divide each rate by the total rate
    NumericVector partialRates = rates / totalRate;
    // cumulative sum the rates to set up the "wheel of fortune"
    NumericVector wheel(rates.length());
    std::partial_sum(partialRates.begin(), partialRates.end(), wheel.begin());
    // generate a random uniform
    rand = runif(1)[0]; 
    // identify which event is happening
    whichWheel = ifelse(rand > wheel, 1, 0);
    event = std::accumulate(whichWheel.begin(), whichWheel.end(), 0);
    
    // increment time
    t += rexp(1, totalRate)[0];
    
    // Production of Th1 
    if (event==0) 
      T1 += 1;
    else if (event==1)
      T2 += 1;
    else if (event==2)
      T1 -= 1;
    else if (event==3)
      T2 -= 1;
    else if (event==4)
      P += 1;
    else 
      P -= 1;
    
  }
  // Get the final system state
  Popn(tIter, 1) = T1;
  Popn(tIter, 2) = T2;
  Popn(tIter, 3) = P;
  
  return Popn;
}