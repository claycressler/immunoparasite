#include <Rcpp.h>
using namespace Rcpp;

// a utility function for removing a row from a matrix
// this works equivalently to R's x[-i,]
NumericMatrix row_erase (NumericMatrix& x, int& rowID) {
  // a similar function would exist for removing a column.
  NumericMatrix x2(Dimension(x.nrow()-1, x.ncol()));
  int iter = 0; // possibly make this a pointer?
  for (int i = 0; i < x.nrow(); i++) {
    if (i != rowID) {
      x2.row(iter) = x.row(i);
      iter++;
    }
  }
  return x2;
}

// [[Rcpp::export]]
List nested_modelC(NumericVector params) {
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
  double c = params[16]; // contact rate between hosts
  double v = params[17]; // virulence (per-capita mortality rate for infected hosts - I assume that parasite load does not determine virulence, but that virulence does determine parasite replication rate)
  double v0 = params[18]; // half-saturation constant scaling virulence into parasite replication rate
  double vCV = params[19]; // coefficient of variation in virulence (for evolving simulations)
  double tmax = params[20]; // length of time to run simulation
  int S0 = params[21]; // initial number of susceptible hosts
  int I0 = params[22]; // initial number of infected hosts
  int R = 0; // initial number of recoveries
  
  // Set up the host population as a matrix, since the total number of hosts cannot grow
  NumericMatrix Hosts(S0+I0, 5);
  for (int i=0; i < S0; i++) { // set the initial state of the susceptible hosts
    Hosts(i,0) = 0; // initial Th1 
    Hosts(i,1) = 0; // initial Th2
    Hosts(i,2) = 0; // initial P
    Hosts(i,3) = 0; // initial v
    Hosts(i,4) = i; // individual ID for tracking through time
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
    Hosts(i,4) = i; // individual ID for tracking through time
  }
  
  // initialize time
  double t = 0.0; 
  
  // initialize storage for population statistics (no. susceptible, infected, recovered, mean virulence)
  NumericMatrix Popn(tmax * 10 + 1, 5);
  NumericVector ones = rep(1.0,Popn.nrow());
  ones(0) = 0.0;
  NumericVector times(ones.length());
  std::partial_sum(ones.begin(), ones.end(), times.begin());
  times = times/10;
  Popn(_,0) = times;
  int tIter = 0;
  // initialize storage for the states of the individual hosts
  List HostState(Popn.nrow());
  
  while(t < tmax) {
    // get the current states for all hosts
    NumericVector T1 = Hosts(_,0);
    NumericVector T2 = Hosts(_,1);
    NumericVector P = Hosts(_,2);
    NumericVector V = Hosts(_,3);
    
    // get the current population state
    double S = std::count(P.begin(), P.end(), 0.0);
    double I = Hosts.nrow()-S;
    
    // store population information?
    if (t >= times(tIter)) {
      Popn(tIter, 1) = S;
      Popn(tIter, 2) = I;
      Popn(tIter, 3) = R;
      Popn(tIter, 4) = mean(V);
      HostState(tIter) = clone(Hosts);
      tIter += 1;
    }
    
    // compute the rates of all model processes
    NumericVector prod1 = ifelse(P > 0.0, b1 + c1*P/(C1+P) + s1*pow(T1,2.0)/(pow(S1,2.0)+pow(T1,2.0)) * I12/(I12+T2), b1);
    NumericVector prod2 = ifelse(P > 0.0, b2 + c2*P/(C2+P) + s2*pow(T2,2.0)/(pow(S2,2.0)+pow(T2,2.0)) * I21/(I21+T1), b2);
    NumericVector decay1 = m*T1;
    NumericVector decay2 = m*T2;
    NumericVector birthP = bp*V/(v0+V)*P*(1-P/Kp);
    NumericVector deathP = a*T2*P;
    NumericVector deathI = V;
    double contact = c * S * I;
    
    // combine all rates into a single vector
    // this looks very tedious but is 30x faster than similar R code!
    // compute the total number of rates
    int nrates = prod1.length() + prod2.length() + decay1.length() + decay2.length() + birthP.length() + deathP.length() + deathI.length() + 1; 
    // set up the storage vector
    NumericVector rates(nrates);
    // add each set of rate processes to the storage vector
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
    double totalRate = std::accumulate(rates.begin(), rates.end(), 0.0);
    
    //Rcout << rates << "\n";
    
    // divide each rate by the total rate
    NumericVector partialRates = rates / totalRate;
    // cumulative sum the rates to set up the "wheel of fortune"
    NumericVector wheel(rates.length());
    std::partial_sum(partialRates.begin(), partialRates.end(), wheel.begin());
    
    // generate a random uniform
    double rand = runif(1)[0]; 
    // identify which event is happening
    IntegerVector whichWheel = ifelse(rand > wheel, 1, 0);
    int event = std::accumulate(whichWheel.begin(), whichWheel.end(), 0);

    // increment time
    t += rexp(1, totalRate)[0];

    // Production of Th1 
    if (event < prod1.length()) {
      int ind = event; 
      Hosts(ind,0) += 1.0;
    }
    // Production of Th2 
    else if (event >= prod1.length() && event < (prod1.length() + prod2.length())) {
      int ind = event-prod1.length(); 
      Hosts(ind,1) += 1.0;
    }
    // Decay of Th1
    else if (event >= (prod1.length() + prod2.length()) && event < (prod1.length() + prod2.length() + decay1.length())) {
      int ind = event-prod1.length()-prod2.length(); 
      Hosts(ind,0) -= 1.0;
    }
    // Decay of Th2
    else if (event >= (prod1.length() + prod2.length() + decay1.length()) && event < (prod1.length() + prod2.length() + decay1.length() + decay2.length())) {
      int ind = event-prod1.length()-prod2.length()-decay1.length(); 
      Hosts(ind,1) -= 1.0;
    }
    // Parasite birth
    else if (event >= (prod1.length() + prod2.length() + decay1.length() + decay2.length()) && event < (prod1.length() + prod2.length() + decay1.length() + decay2.length() + birthP.length())) {
      int ind = event-prod1.length()-prod2.length()-decay1.length()-decay2.length(); 
      Hosts(ind,2) += 1.0;
    }
    // Parasite death
    else if (event >= (prod1.length() + prod2.length() + decay1.length() + decay2.length() + birthP.length()) && event < (prod1.length() + prod2.length() + decay1.length() + decay2.length() + birthP.length() + deathP.length())) {
      int ind = event-prod1.length()-prod2.length()-decay1.length()-decay2.length()-birthP.length(); 
      Hosts(ind,2) -= 1.0;
      // if P = 0 now, set virulence = 0 and update the number of recoveries
      if (Hosts(ind,2) < 1.0) {
        Hosts(ind,3) = 0.0;
        R += 1.0;
      }
    }
    // Host death
    else if (event >= (prod1.length() + prod2.length() + decay1.length() + decay2.length() + birthP.length() + deathP.length()) && event < (prod1.length() + prod2.length() + decay1.length() + decay2.length() + birthP.length() + deathP.length() + deathI.length())) {
      int ind = event-prod1.length()-prod2.length()-decay1.length()-decay2.length()-birthP.length()-deathP.length(); 
      // Remove this host from the Hosts NumericMatrix
      Hosts = row_erase(Hosts, ind);
    }
    // Infection
    else {
      // choose a susceptible host and an infected host at random
      // this can be done in the same way as we chose the event
      // first identify all of the susceptible hosts
      NumericVector Shosts = ifelse(Hosts(_,2)>0.0, 0.0, 1.0);
      // divide by the total number of susceptible hosts
      Shosts = Shosts / (std::accumulate(Shosts.begin(), Shosts.end(), 0.0));
      // cumulative sum to set up the "wheel of fortune" for susceptible hosts
      NumericVector Swheel(Shosts.length());
      std::partial_sum(Shosts.begin(), Shosts.end(), Swheel.begin());
      // choose a random uniform
      double randS = runif(1)[0];
      // identify the index of the randomly chosen S
      IntegerVector whichShosts = ifelse(randS > Swheel, 1, 0);
      int Sind = std::accumulate(whichShosts.begin(), whichShosts.end(), 0);

      // first identify all of the infected hosts
      NumericVector Ihosts = ifelse(Hosts(_,2)>0.0, 1.0, 0.0);
      // divide by the total number of infected hosts
      Ihosts = Ihosts / (std::accumulate(Ihosts.begin(), Ihosts.end(), 0.0));
      // cumulative sum to set up the "wheel of fortune" for infected hosts
      NumericVector Iwheel(Ihosts.length());
      std::partial_sum(Ihosts.begin(), Ihosts.end(), Iwheel.begin());
      // choose a random uniform
      double randI = runif(1)[0];
      // identify the index of the randomly chosen I
      IntegerVector whichIhosts = ifelse(randI > Iwheel, 1, 0);
      int Iind = std::accumulate(whichIhosts.begin(), whichIhosts.end(), 0);
      
      // set the initial dose and virulence of the newly infected individual on the basis of the state of the infecting host
      Hosts(Sind,2) = rpois(1, Hosts(Iind,2)/10)[0]; // initial P is drawn from a Poisson distribution
      // if the initial dose is 0, don't set the virulence - this individual escaped without an infection!
      if (Hosts(Sind,2) > 0.0) {
        // initial virulence is drawn from a lognormal distribution
        double vSD = vCV * Hosts(Iind,3); 
        double mu = log(pow(Hosts(Iind,3),2.0) / sqrt(pow(vSD,2.0) + pow(Hosts(Iind,3),2.0)));
        double sigma = sqrt(log(pow(vSD,2.0) / pow(Hosts(Iind,3),2.0) + 1));
        Hosts(Sind,3) = rlnorm(1, mu, sigma)[0];
      }
    }
  }
  // Get the final system state
  NumericVector P = Hosts(_,2);
  NumericVector V = Hosts(_,3);
  Popn(tIter, 1) = std::count(P.begin(), P.end(), 0.0);
  Popn(tIter, 2) = Hosts.nrow()-Popn(tIter, 1);
  Popn(tIter, 3) = R;
  Popn(tIter, 4) = mean(V);
  //HostState[tIter] = Hosts;
  
  // return a list
  return List::create(HostState,Popn);
}