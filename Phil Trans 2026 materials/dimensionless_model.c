#include <R.h>
static double parms[12];
#define b1 parms[0]
#define b2 parms[1]
#define s1 parms[2]
#define s2 parms[3]
#define i12 parms[4]
#define i21 parms[5]
#define k1 parms[6]
#define k2 parms[7]
#define c1 parms[8]
#define c2 parms[9]
#define r parms[10]
#define a parms[11]

/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
  int N = 12;
  odeparms(&N, parms);
}

/* Derivatives */
void derivs (int *neq, double *t, double *y, double *ydot,
             double *yout, int *ip)
{
  ydot[0] = b1 + c1*y[2]/(k1+y[2]) + s1*y[0]*y[0]/(1+y[0]*y[0]) * i12/(i12+y[1]) - y[0];
  ydot[1] = b2 + c2*y[2]/(k2+y[2]) + s2*y[1]*y[1]/(1+y[1]*y[1]) * i21/(i21+y[0]) - y[1];
  ydot[2] = r*y[2]*(1-y[2]) - a*y[2]*y[1];
  
}
