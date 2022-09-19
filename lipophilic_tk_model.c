/* lipophilic_tk_model.c for R deSolve package
   ___________________________________________________

   Model File:  lipophilic_tk.model

   Date:  Mon Sep 12 16:28:18 2022

   Created by:  "C:/Users/mmarikar/ONEDRI~1/Profile/Desktop/myR/DUSTIN~1/KAPRAU~1/KAPRAU~1/mod.exe v6.1.0"
    -- a model preprocessor by Don Maszle
   ___________________________________________________

   Copyright (c) 1993-2019 Free Software Foundation, Inc.

   Model calculations for compartmental model:

   8 States:
     A_mf = 0.0,
     A_i = 0.0,
     AUC_m = 0.0,
     AUC_i = 0.0,
     d_m = 0.0,
     d_i = 0.0,
     T_in = 0.0,
     T_out = 0.0,

   9 Outputs:
    "C_m",
    "C_i",
    "M_mf",
    "M_i",
    "D_m",
    "D_i",
    "R_milk",
    "r_m_mf",
    "A_bal",

   4 Inputs:
     M_mf_in (forcing function)
     M_i_in (forcing function)
     R_milk_in (forcing function)
     r_m_mf_in (forcing function)

   8 Parameters:
     half_life = 70.0,
     F_m = 0.094,
     F_milk = 0.154,
     r_f_m = 0.35,
     F_abs = 0.9,
     n_i = 10,
     food_dose = 0,
     k = 0.0,
*/

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/* Model variables: States */
#define ID_A_mf 0x00000
#define ID_A_i 0x00001
#define ID_AUC_m 0x00002
#define ID_AUC_i 0x00003
#define ID_d_m 0x00004
#define ID_d_i 0x00005
#define ID_T_in 0x00006
#define ID_T_out 0x00007

/* Model variables: Outputs */
#define ID_C_m 0x00000
#define ID_C_i 0x00001
#define ID_M_mf 0x00002
#define ID_M_i 0x00003
#define ID_D_m 0x00004
#define ID_D_i 0x00005
#define ID_R_milk 0x00006
#define ID_r_m_mf 0x00007
#define ID_A_bal 0x00008

/* Parameters */
static double parms[8];

#define half_life parms[0]
#define F_m parms[1]
#define F_milk parms[2]
#define r_f_m parms[3]
#define F_abs parms[4]
#define n_i parms[5]
#define food_dose parms[6]
#define k parms[7]

/* Forcing (Input) functions */
static double forc[4];

#define M_mf_in forc[0]
#define M_i_in forc[1]
#define R_milk_in forc[2]
#define r_m_mf_in forc[3]

/* Function definitions for delay differential equations */

int Nout=1;
int nr[1]={0};
double ytau[1] = {0.0};

static double yini[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; /*Array of initial state variables*/

void lagvalue(double T, int *nr, int N, double *ytau) {
  static void(*fun)(double, int*, int, double*) = NULL;
  if (fun == NULL)
    fun = (void(*)(double, int*, int, double*))R_GetCCallable("deSolve", "lagvalue");
  return fun(T, nr, N, ytau);
}

double CalcDelay(int hvar, double dTime, double delay) {
  double T = dTime-delay;
  if (dTime > delay){
    nr[0] = hvar;
    lagvalue( T, nr, Nout, ytau );
}
  else{
    ytau[0] = yini[hvar];
}
  return(ytau[0]);
}

/*----- Initializers */
void initmod (void (* odeparms)(int *, double *))
{
  int N=8;
  odeparms(&N, parms);
}

void initforc (void (* odeforcs)(int *, double *))
{
  int N=4;
  odeforcs(&N, forc);
}


/* Calling R code will ensure that input y has same
   dimension as yini */
void initState (double *y)
{
  int i;

  for (i = 0; i < sizeof(yini) / sizeof(yini[0]); i++)
  {
    yini[i] = y[i];
  }
}

void getParms (double *inParms, double *out, int *nout) {
/*----- Model scaling */

  int i;

  for (i = 0; i < *nout; i++) {
    parms[i] = inParms[i];
  }


  k = log ( 2 ) / half_life ;

  for (i = 0; i < *nout; i++) {
    out[i] = parms[i];
  }
  }
/*----- Dynamics section */

void derivs (int *neq, double *pdTime, double *y, double *ydot, double *yout, int *ip)
{

  yout[ID_M_mf] = M_mf_in ;

  yout[ID_M_i] = M_i_in ;

  yout[ID_R_milk] = R_milk_in ;

  yout[ID_r_m_mf] = r_m_mf_in ;

  yout[ID_C_m] = y[ID_A_mf] / yout[ID_M_mf] * ( yout[ID_M_i] > 0 ? 1.0 : yout[ID_r_m_mf] ) ;

  yout[ID_C_i] = ( yout[ID_M_i] > 0 ? y[ID_A_i] / yout[ID_M_i] : r_f_m * yout[ID_C_m] ) ;

  yout[ID_C_i] = ( (*pdTime) < 0 ? 0 : yout[ID_C_i] ) ;

  yout[ID_D_m] = ( food_dose ? 0.065 * pow ( yout[ID_M_mf] , 0.7919 ) / yout[ID_M_mf] : 1.0 ) * y[ID_d_m] ;

  yout[ID_D_i] = ( food_dose ? 0.065 * pow ( yout[ID_M_i] , 0.7919 ) / ( yout[ID_M_i] > 0 ? yout[ID_M_i] : 1.0e-8 ) : 1.0 ) * y[ID_d_i] ;

  yout[ID_A_bal] = y[ID_T_in] - y[ID_A_mf] - n_i * y[ID_A_i] - y[ID_T_out] ;

  ydot[ID_A_mf] = F_abs * yout[ID_D_m] * yout[ID_M_mf] - k * y[ID_A_mf] - n_i * yout[ID_R_milk] * F_milk * y[ID_A_mf] / ( F_m * yout[ID_M_mf] ) ;

  ydot[ID_A_i] = F_abs * yout[ID_D_i] * yout[ID_M_i] - k * y[ID_A_i] + yout[ID_R_milk] * F_milk * y[ID_A_mf] / ( F_m * yout[ID_M_mf] ) ;

  ydot[ID_AUC_m] = yout[ID_C_m] ;

  ydot[ID_AUC_i] = yout[ID_C_i] ;

  ydot[ID_d_m] = 0.0 ;

  ydot[ID_d_i] = 0.0 ;

  ydot[ID_T_in] = F_abs * yout[ID_D_m] * yout[ID_M_mf] + n_i * F_abs * yout[ID_D_i] * yout[ID_M_i] ;

  ydot[ID_T_out] = k * ( y[ID_A_mf] + n_i * y[ID_A_i] ) ;

} /* derivs */


/*----- Jacobian calculations: */
void jac (int *neq, double *t, double *y, int *ml, int *mu, double *pd, int *nrowpd, double *yout, int *ip)
{

} /* jac */


/*----- Events calculations: */
void event (int *n, double *t, double *y)
{

} /* event */

/*----- Roots calculations: */
void root (int *neq, double *t, double *y, int *ng, double *gout, double *out, int *ip)
{

} /* root */

