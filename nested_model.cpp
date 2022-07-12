#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
//NumericMatrix nested_model(NumericVector params) {
NumericVector nested_model(NumericVector params) {
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
  double c = params[16]; // contact rate between hosts
  double v = params[17]; // virulence (per-capita mortality rate for infected hosts - I assume that parasite load does not determine virulence, but that virulence does determine parasite replication rate)
  double v0 = params[18]; // half-saturation constant scaling virulence into parasite replication rate
  double vCV = params[19]; // coefficient of variation in virulence (for evolving simulations)
  double tmax = params[20]; // length of time to run simulation
  int S0 = params[21]; // initial number of susceptible hosts
  int I0 = params[22]; // initial number of infected hosts
  int R = 0; // initial number of recoveries
  
  // Set up the host population as a matrix, since the total number of hosts cannot grow
  NumericMatrix Hosts(S0+I0, 4);
  for (int i=0; i < S0; i++) { // set the initial state of the susceptible hosts
    Hosts(i,0) = 0; // initial Th1 
    Hosts(i,1) = 0; // initial Th2
    Hosts(i,2) = 0; // initial P
    Hosts(i,3) = 0; // initial v
  }
  for (int i=S0; i < S0+I0; i++) { // set the initial state of the susceptible hosts
    Hosts(i,0) = 0; // initial Th1 
    Hosts(i,1) = 0; // initial Th2
    Hosts(i,2) = rpois(1, Kp/10)[0]; // initial P is drawn from a Poisson distribution
    // initial virulence is drawn from a lognormal distribution
    double vSD = vCV * v; 
    double mu = log(pow(v,2.0) / sqrt(pow(vSD,2.0) + pow(v,2.0)));
    double sigma = sqrt(log(pow(vSD,2.0) / pow(v,2.0) + 1));
    Hosts(i,3) = rlnorm(1, mu, sigma)[0];
  }
  
  // initialize time
  double t = 0.0; 
  
  // initialize storage for population statistics (no. susceptible, infected, recovered, mean virulence)
  NumericMatrix Popn(tmax * 10, 4); 
  
  while(t < tmax) {
    // get the current states for all hosts
    NumericVector T1 = Hosts(_,0);
    NumericVector T2 = Hosts(_,1);
    NumericVector P = Hosts(_,2);
    NumericVector V = Hosts(_,3);
    
    // get the current population state
    int S = std::count(P.begin(), P.end(), 0);
    int I = Hosts.nrow()-S;
    
    NumericVector prod1 = ifelse(P > 0, b1 + c1*P/(C1+P) + s1*pow(T1,2.0)/(pow(S1,2.0)+pow(T1,2.0)) * I12/(I12+T2), b1);
    NumericVector prod2 = ifelse(P > 0, b2 + c2*P/(C2+P) + s2*pow(T2,2.0)/(pow(S2,2.0)+pow(T2,2.0)) * I21/(I21+T1), b2);
    NumericVector decay1 = m*T1;
    NumericVector decay2 = m*T2;
    NumericVector birthP = bp*V/(v0+V)*P*(1-P/Kp);
    NumericVector deathP = a*T2*P;
    NumericVector deathI = V;
    double contact = c * S * I;
    
    int nrates = prod1.length() + prod2.length() + decay1.length() + decay2.length() + birthP.length() + deathP.length() + deathI.length() + 1; 
    NumericVector rates(nrates);
    for (int i=0; i < prod1.length(); i++) {
      rates(i) = prod1(i);
      rates(i+prod1.length()) = prod2(i);
      rates(i+prod1.length()+prod2.length()) = decay1(i);
      rates(i+prod1.length()+prod2.length()+decay1.length()) = decay2(i);
      rates(i+prod1.length()+prod2.length()+decay1.length()+decay2.length()) = birthP(i);
      rates(i+prod1.length()+prod2.length()+decay1.length()+decay2.length()+birthP.length()) = deathP(i);
      rates(i+prod1.length()+prod2.length()+decay1.length()+decay2.length()+birthP.length()+deathP.length()) = deathI(i);
    }
    rates(prod1.length()+prod2.length()+decay1.length()+decay2.length()+birthP.length()+deathP.length()+deathI.length()) = contact;
    
 
  }
  return rates;
}