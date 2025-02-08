//--------------------------------------------------------------------------------------------------
//  For standalone testing (uncomment main block first)
//  gcc -o gamma_integral gamma_integral.c -lm
//
//  For use with IDL's call external routine
//  gcc -Wall -fPIC -c -I../../include gamma_integral.c
//  gcc -shared -Wl,-soname,gamma_integral.so.1 -o gamma_integral.so.1.0   *.o
//
//  Calling from within IDL: 
//  result = call_external('gamma_integral.so.1.0','gamma_integral',p1,p2,chi0,eta,mu,eps,i,err,nfunc)
//--------------------------------------------------------------------------------------------------
#include <math.h>
#include <stdio.h>

#include "gamma_integral.h"
#include "intde1.h"

//#include "/usr/local/itt/idl71/external/include/idl_export.h" /* IDL external definitions */
//#include "/usr/local/exelis/idl/external/include/idl_export.h" /* IDL external definitions */

static long nfunc;

static double p;
static double eta;
static double mu;

double fgam(double);
double fgam_prime(double);
double fgam_asymptotic(double);

double fgam1(double);
double fgam2(double);
double fgam3(double);

double fgam(double x)
{
    extern long nfunc;
    double trig;
    
    nfunc++;
   
    trig = sin(x*eta/2.0)*sin(x*mu/2.0);
    return pow(x,-p)*trig*trig;
}

double fgam_prime(double x)
{
    extern long nfunc;
    double c,r,trig;
    
    nfunc++;
	
	if (eta >= mu){
      c = pow(2.0/mu,1.0-p);
      r = eta/mu;
	} else {
      c = pow(2.0/eta,1.0-p);
      r = mu/eta;
	}
   
    trig = sin(x)*sin(x*r);
    return c*pow(x,-p)*trig*trig;
}

double fgam_asymptotic(double x)
{
    extern long nfunc;
    double c,r, trig, term;
    
    nfunc++;
	
	if (eta >= mu){
      c = pow(2.0/mu,1.0-p);
      r = eta/mu;
	} else {
      c = pow(2.0/eta,1.0-p);
      r = mu/eta;
	}
   
    trig = sin(x);
	term = pow(r*x,2)/(1.0+2.0*pow(r*x,2));
    return c*pow(x,-p)*trig*trig*term;
}

double fgam_fluctuations(double x)
{
    extern long nfunc;
    
    nfunc++;

    return fgam_prime(x) - fgam_asymptotic(x);
}


void get_eta_mu(double eta, double mu, double *etap, double *mup)
{
   if (eta < mu){
       // mu is the larger
       *etap = eta;
       *mup = mu;
   }else{
      // exchange eta and mu so that mu is the larger
       *mup = eta;
      *etap = mu;
   }
   return;
}

// ----------------------------------------------------------------
// we split the integrand fgam so that fgam = fgam1 + fgam2 + fgam3
// ----------------------------------------------------------------

// oscillates primarily with wavenumber mu
double fgam1(double x)
{
    extern long nfunc;
    double etap, mup;
    double argm, arge;
    double result;
     
    nfunc++;
   
    get_eta_mu(eta, mu, &etap, &mup);

    argm = x*mup /2.0;
    arge = x*etap/2.0;
                 
    result =  pow(x,-p) * pow(sin(arge),2.0) * 
      ( pow(sin(argm),2.0) - pow(argm,2.0) / (1.0 + 2.0*pow(argm,2.0) ) );
    //result = 2.D0 * 8.D0  / (2.D0*!dpi) * result
    return result;
}

// oscillates primarily with wavenumber eta
double fgam2(double x)
{
    extern long nfunc;
    double etap, mup;
    double argm, arge;
    double result;
    
    nfunc++;
   
    get_eta_mu(eta, mu, &etap, &mup);

    argm = x*mup /2.0;
    arge = x*etap/2.0;
    
    result =  pow(x,-p)  * pow(argm,2.0) / (1.0 + 2.0*pow(argm,2.0) ) * 
     (  pow(sin(arge),2.0) - pow(arge,2.0) / (1.0 + 2.0*pow(arge,2.0) ) );
     
    //result = 2.D0 * 8.D0  / (2.D0*!dpi) * result
    return result;
}

// does not oscillate
double fgam3(double x)
{
    extern long nfunc;
    double etap, mup;
    double argm, arge;
    double result;
    
    nfunc++;
   
    get_eta_mu(eta, mu, &etap, &mup);

    argm = x*mup /2.0;
    arge = x*etap/2.0;
    
    result = pow(x,-p) * pow(arge,2.0) / (1.0 + 2.0*pow(arge,2.0))  * pow(argm,2.0) / (1.0 + 2.0*pow(argm,2.0));
     
    //result = 2.D0 * 8.D0  / (2.D0*!dpi) * result
    return result;
}

int gamma_integral_natural(double *p1, double *p2, double *chi0, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out){

    double chi0_,eps_;
    double i1, i2, err1, err2;
    long nfunc1, nfunc2;

    extern long nfunc;

    eta   = *eta_in;
    mu    = *mu_in;
    eps_  = *eps;
    chi0_ = *chi0;

    //intdei(fgam, 0.0, eps, &i, &err);

    if (*p1 == *p2){
       p = *p1;
       nfunc = 0;
       intdei(fgam, 0.0, eps_, i, err);
       //printf("    I_gam= %lg, err= %lg, N= %ld\n", *i, *err, nfunc);
       *nfunc_out = nfunc;

    }else{

       // If this branch proves to be significantly faster, we could use for p1=p2 case also

       p = *p1;
       nfunc = 0;
       intde(fgam, 0.0, chi0_, eps_, &i1, &err1);
       nfunc1 = nfunc;
       //printf("    I_gam= %lg, err= %lg, N= %ld\n", i1, err1, nfunc1);

       #if 1 
       // we should decide which branch to take based on value of chi0?

       p = *p2;
       nfunc = 0;
       intdei(fgam, chi0_, eps_, &i2, &err2);
       nfunc2 = nfunc;
       //printf("    I_gam= %lg, err= %lg, N= %ld\n", i2, err2, nfunc2);
       #else
       {
       double ia, ib, erra, errb;
       long nfunca, nfuncb;
       p = *p2;
       nfunc = 0;
       intdei(fgam, 0.0, eps_, &ia, &erra);
       nfunca = nfunc;
       //printf("    I_gam= %lg, err= %lg, N= %ld\n", ia, erra, nfunca);
       nfunc = 0;
       intde(fgam, 0.0, chi0_, eps_, &ib, &errb);
       nfuncb = nfunc;
       //printf("    I_gam= %lg, err= %lg, N= %ld\n", ib, errb, nfuncb);
       i2 = ia - ib;
       err2 = fabs(erra) + fabs(errb);
       if ((erra < 0.0) || (errb < 0.0)){ *err = -(*err); } // make err<0 if error condition
       }
       #endif

       *i = i1 + i2;
       *err = fabs(err1) + fabs(err2);
       if ((err1 < 0.0) || (err2 < 0.0)){ *err = -(*err); } // make err<0 if error condition
       *nfunc_out = nfunc1 + nfunc2;
    }

    return 1;
}

int gamma_integral(int argc, void *argv[]){
   if (argc!=9) {
      printf("Error: calling gamma_prime with incorrect number of arguments\n");
      return 0;
   }

   return gamma_integral_natural((double *) argv[0], (double *) argv[1], (double *) argv[2], (double *) argv[3], 
        (double *) argv[4], (double *) argv[5], (double *) argv[6], (double *) argv[7], (IDL_LONG *) argv[8]);

}

int gamma_infinite_natural(double *p1, double *a, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out){

    double a_,eps_;

    extern long nfunc;

    eta   = *eta_in;
    mu    = *mu_in;
    eps_  = *eps;
    a_    = *a;

       p = *p1;
       nfunc = 0;
       intdei(fgam, a_, eps_, i, err);
       //printf("    I_gam= %lg, err= %lg, N= %ld\n", *i, *err, nfunc);
       *nfunc_out = nfunc;

    return 1;
}

int gamma_infinite(int argc, void *argv[]){
   if (argc!=8) {
      printf("Error: calling gamma_infinite with incorrect number of arguments\n");
      return 0;
   }

   return gamma_infinite_natural((double *) argv[0], (double *) argv[1], (double *) argv[2], 
        (double *) argv[3], (double *) argv[4], (double *) argv[5], (double *) argv[6], (IDL_LONG *) argv[7]);

}


int gamma_prime_infinite_natural(double *p1, double *a, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out){

    double a_,eps_;
    //double i1, i2, err1, err2;
    //long nfunc1, nfunc2;

    extern long nfunc;

    eta   = *eta_in;
    mu    = *mu_in;
    eps_  = *eps;
    a_    = *a;

    //intdei(fgam, 0.0, eps, &i, &err);

       p = *p1;
       nfunc = 0;
       intdei(fgam_prime, a_, eps_, i, err);
       //printf("    I_gam= %lg, err= %lg, N= %ld\n", *i, *err, nfunc);
       *nfunc_out = nfunc;

    return 1;
}

int gamma_prime_infinite(int argc, void *argv[]){
   if (argc!=8) {
      printf("Error: calling gamma_prime with incorrect number of arguments\n");
      return 0;
   }

   return gamma_prime_infinite_natural((double *) argv[0], (double *) argv[1], (double *) argv[2], 
        (double *) argv[3], (double *) argv[4], (double *) argv[5], (double *) argv[6], (IDL_LONG *) argv[7]);

}

int gamma_finite_natural(double *p1, double *a, double *b, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out){

    double a_,b_,eps_;
    //double i1, i2, err1, err2;
    //long nfunc1, nfunc2;

    extern long nfunc;

    eta   = *eta_in;
    mu    = *mu_in;
    eps_  = *eps;
    a_    = *a;
    b_    = *b;

    //intdei(fgam, 0.0, eps, &i, &err);

       p = *p1;
       nfunc = 0;
       intde(fgam, a_, b_, eps_, i, err);

       //printf("    I_gam= %lg, err= %lg, N= %ld\n", *i, *err, nfunc);
       *nfunc_out = nfunc;

    return 1;
}

int gamma_finite(int argc, void *argv[]){
   if (argc!=9) {
      printf("Error: calling gamma_finite with incorrect number of arguments\n");
      return 0;
   }

   return gamma_finite_natural((double *) argv[0], (double *) argv[1], (double *) argv[2], 
        (double *) argv[3], (double *) argv[4], (double *) argv[5], (double *) argv[6],  (double *) argv[7], (IDL_LONG *) argv[8]);

}

int gamma_prime_finite_natural(double *p1, double *a, double *b, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out){

    double a_,b_,eps_;
    //double i1, i2, err1, err2;
    //long nfunc1, nfunc2;

    extern long nfunc;

    eta   = *eta_in;
    mu    = *mu_in;
    eps_  = *eps;
    a_    = *a;
    b_    = *b;

    //intdei(fgam, 0.0, eps, &i, &err);

       p = *p1;
       nfunc = 0;
	   intde(fgam_prime, a_, b_, eps_, i, err);

       //printf("    I_gam= %lg, err= %lg, N= %ld\n", *i, *err, nfunc);
       *nfunc_out = nfunc;

    return 1;
}

int gamma_prime_finite(int argc, void *argv[]){
   if (argc!=9) {
      printf("Error: calling gamma_prime with incorrect number of arguments\n");
      return 0;
   }

   return gamma_prime_finite_natural((double *) argv[0], (double *) argv[1], (double *) argv[2], 
        (double *) argv[3], (double *) argv[4], (double *) argv[5], (double *) argv[6],  (double *) argv[7], (IDL_LONG *) argv[8]);

}

int gamma_asymptotic_infinite_natural(double *p1, double *a, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out){

    double a_,eps_;
    //double i1, i2, err1, err2;
    //long nfunc1, nfunc2;

    extern long nfunc;

    eta   = *eta_in;
    mu    = *mu_in;
    eps_  = *eps;
    a_    = *a;

    //intdei(fgam, 0.0, eps, &i, &err);

       p = *p1;
       nfunc = 0;
       intdei(fgam_asymptotic, a_, eps_, i, err);
       //printf("    I_gam= %lg, err= %lg, N= %ld\n", *i, *err, nfunc);
       *nfunc_out = nfunc;

    return 1;
}

int gamma_asymptotic_infinite(int argc, void *argv[]){
   if (argc!=8) {
      printf("Error: calling gamma_asymptotic with incorrect number of arguments\n");
      return 0;
   }
   
   return gamma_asymptotic_infinite_natural((double *) argv[0], (double *) argv[1], (double *) argv[2], 
        (double *) argv[3], (double *) argv[4], (double *) argv[5], (double *) argv[6], (IDL_LONG *) argv[7]);
}


int gamma_asymptotic_finite_natural(double *p1, double *a, double *b, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out){

    double a_,b_,eps_;
    //double i1, i2, err1, err2;
    //long nfunc1, nfunc2;

    extern long nfunc;

    eta   = *eta_in;
    mu    = *mu_in;
    eps_  = *eps;
    a_    = *a;
    b_    = *b;

    //intdei(fgam, 0.0, eps, &i, &err);

       p = *p1;
       nfunc = 0;
	   intde(fgam_asymptotic, a_, b_, eps_, i, err);

       //printf("    I_gam= %lg, err= %lg, N= %ld\n", *i, *err, nfunc);
       *nfunc_out = nfunc;

    return 1;
}

int gamma_asymptotic_finite(int argc, void *argv[]){
   if (argc!=9) {
      printf("Error: calling gamma_prime with incorrect number of arguments\n");
      return 0;
   }

   return gamma_asymptotic_finite_natural((double *) argv[0], (double *) argv[1], (double *) argv[2], 
        (double *) argv[3], (double *) argv[4], (double *) argv[5], (double *) argv[6], (double *) argv[7], (IDL_LONG *) argv[8]);

}

int gamma_fluctuations_infinite_natural(double *p1, double *a, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out){

    double a_,eps_;
    //double i1, i2, err1, err2;
    //long nfunc1, nfunc2;

    extern long nfunc;

    eta   = *eta_in;
    mu    = *mu_in;
    eps_  = *eps;
    a_    = *a;

    //intdei(fgam, 0.0, eps, &i, &err);

       p = *p1;
       nfunc = 0;
       intdei(fgam_fluctuations, a_, eps_, i, err);
       //printf("    I_gam= %lg, err= %lg, N= %ld\n", *i, *err, nfunc);
       *nfunc_out = nfunc;

    return 1;
}

int gamma_fluctuations_infinite(int argc, void *argv[]){
   if (argc!=8) {
      printf("Error: calling gamma_asymptotic with incorrect number of arguments\n");
      return 0;
   }
   
   return gamma_fluctuations_infinite_natural((double *) argv[0], (double *) argv[1], (double *) argv[2], 
        (double *) argv[3], (double *) argv[4], (double *) argv[5], (double *) argv[6], (IDL_LONG *) argv[7]);
}


int gamma_fluctuations_finite_natural(double *p1, double *a, double *b, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out){

    double a_,b_,eps_;
    //double i1, i2, err1, err2;
    //long nfunc1, nfunc2;

    extern long nfunc;

    eta   = *eta_in;
    mu    = *mu_in;
    eps_  = *eps;
    a_    = *a;
    b_    = *b;

    //intdei(fgam, 0.0, eps, &i, &err);

       p = *p1;
       nfunc = 0;
	   intde(fgam_fluctuations, a_, b_, eps_, i, err);

       //printf("    I_gam= %lg, err= %lg, N= %ld\n", *i, *err, nfunc);
       *nfunc_out = nfunc;

    return 1;
}

int gamma_fluctuations_finite(int argc, void *argv[]){
   if (argc!=9) {
      printf("Error: calling gamma_prime with incorrect number of arguments\n");
      return 0;
   }

   return gamma_fluctuations_finite_natural((double *) argv[0], (double *) argv[1], (double *) argv[2], 
        (double *) argv[3], (double *) argv[4], (double *) argv[5], (double *) argv[6], (double *) argv[7], (IDL_LONG *) argv[8]);

}

// -----------------

// This one does not appear to be useful
int gamma_finite_left_natural(double *p1, double *a, double *b, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out){

    double a_,b_,eps_;
    double ia1,ia2,ia3;
    double erra1,erra2,erra3;
    
    extern long nfunc;

    eta   = *eta_in;
    mu    = *mu_in;
    eps_  = *eps;
    a_    = *a;
    b_    = *b;

    p = *p1;
    nfunc = 0;
              
    intde(fgam1, a_, b_, eps_, &ia1, &erra1);
    intde(fgam2, a_, b_, eps_, &ia2, &erra2);
    intde(fgam3, a_, b_, eps_, &ia3, &erra3);
       
    *i = ia1 + ia2 + ia3;
       
    *err = fabs(erra1) + fabs(erra2) + fabs(erra3) ;

    if ((erra1 < 0.0) || (erra2 < 0.0) || (erra3 < 0.0)){ *err = -(*err); } // make err<0 if error condition

    //printf("    I_gam= %lg, err= %lg, N= %ld\n", *i, *err, nfunc);
    *nfunc_out = nfunc;

    return 1;
}

int gamma_finite_left(int argc, void *argv[]){
   if (argc!=9) {
      printf("Error: calling gamma_finite_left with incorrect number of arguments\n");
      return 0;
   }

   return gamma_finite_left_natural((double *) argv[0], (double *) argv[1], (double *) argv[2], 
        (double *) argv[3], (double *) argv[4], (double *) argv[5], (double *) argv[6],  (double *) argv[7], (IDL_LONG *) argv[8]);
}

// requires p1>1 otherwise integral diverges
int gamma_finite_right_natural(double *p1, double *a, double *b, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out){

    double a_,b_,eps_;
    double ia1,ia2,ia3, ib1,ib2,ib3;
    double erra1,erra2,erra3, errb1,errb2,errb3;
    double smaller, larger;

    extern long nfunc;

    eta   = *eta_in;
    mu    = *mu_in;
    eps_  = *eps;
    a_    = *a;
    b_    = *b;

    if (eta < mu){
       // mu is the larger
       smaller = eta;
       larger = mu;
   }else{
      // exchange eta and mu so that mu is the larger
       larger = eta;
       smaller = mu;
   }

       p = *p1;
       nfunc = 0;

       intdeo(fgam1, a_, larger, eps_, &ia1, &erra1);  //printf("nfunc1a: %ld",nfunc);
       intdeo(fgam2, a_, smaller, eps_, &ia2, &erra2); //printf("nfunc1a: %ld",nfunc);
       intdei(fgam3, a_, eps_, &ia3, &erra3);          //printf("nfunc1a: %ld",nfunc);
       
       intdeo(fgam1, b_, larger, eps_, &ib1, &errb1);  //printf("nfunc1b: %ld",nfunc);
       intdeo(fgam2, b_, smaller, eps_, &ib2, &errb2); //printf("nfunc1b: %ld",nfunc);
       intdei(fgam3, b_, eps_, &ib3, &errb3);          //printf("nfunc1b: %ld",nfunc);
       
       *i = (ia1 + ia2 + ia3) - (ib1 + ib2 + ib3);
       
       *err = fabs(erra1) + fabs(erra2) + fabs(erra3) + fabs(errb1) + fabs(errb2) + fabs(errb3) ;

       if ((erra1 < 0.0) || (erra2 < 0.0) || (erra3 < 0.0) || (errb1 < 0.0) || (errb2 < 0.0) || (errb3 < 0.0)){ *err = -(*err); } // make err<0 if error condition

       //printf("    I_gam= %lg, err= %lg, N= %ld\n", *i, *err, nfunc);
       *nfunc_out = nfunc;

    return 1;
}

int gamma_finite_right(int argc, void *argv[]){
   if (argc!=9) {
      printf("Error: calling gamma_finite_right with incorrect number of arguments\n");
      return 0;
   }

   return gamma_finite_right_natural((double *) argv[0], (double *) argv[1], (double *) argv[2], 
        (double *) argv[3], (double *) argv[4], (double *) argv[5], (double *) argv[6],  (double *) argv[7], (IDL_LONG *) argv[8]);
}

int gamma_infinite_right_natural(double *p1, double *a, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out){

    double a_,eps_;
    double ia1,ia2,ia3;
   // double ia31,ia32;
    double erra1,erra2,erra3;
   // double erra31,erra32;
    double smaller, larger;
   // double breakpoint;
    extern long nfunc;

    eta   = *eta_in;
    mu    = *mu_in;
    eps_  = *eps;
    a_    = *a;

    if (eta < mu){
       // mu is the larger
       smaller = eta;
       larger = mu;
   }else{
      // exchange eta and mu so that mu is the larger
       larger = eta;
       smaller = mu;
   }

       p = *p1;
       nfunc = 0;

       intdeo(fgam1, a_, larger, eps_, &ia1, &erra1);  //printf("nfunc1: %ld",nfunc);
       intdeo(fgam2, a_, smaller, eps_, &ia2, &erra2); //printf("nfunc2: %ld",nfunc);
       
     //  breakpoint = M_PI/(2.0*larger);
     //  if (a_ < breakpoint){ // add a breakpoint at pi/(2*larger)
     //     intde(fgam3, a_, breakpoint,eps_, &ia31, &erra31);
     //     intdei(fgam3, breakpoint, eps_, &ia32, &erra32);
     //     ia3 = ia31 + ia32;
     //     erra3 = fabs(erra31) +fabs(erra32);
     //  } else {
          intdei(fgam3, a_, eps_, &ia3, &erra3);          //printf("nfunc3: %ld",nfunc);
     //  }
       
       *i = ia1 + ia2 + ia3;
       
       *err = fabs(erra1) + fabs(erra2) + fabs(erra3) ;

       if ((erra1 < 0.0) || (erra2 < 0.0) || (erra3 < 0.0)){ *err = -(*err); } // make err<0 if error condition

       //printf("    I_gam= %lg, err= %lg, N= %ld\n", *i, *err, nfunc);
       *nfunc_out = nfunc;

    return 1;
}

int gamma_infinite_right(int argc, void *argv[]){
   if (argc!=8) {
      printf("Error: calling gamma_infinite_right with incorrect number of arguments\n");
      return 0;
   }

   return gamma_infinite_right_natural((double *) argv[0], (double *) argv[1], (double *) argv[2], 
        (double *) argv[3], (double *) argv[4], (double *) argv[5], (double *) argv[6], (IDL_LONG *) argv[7]);
}

int expm_natural(double *x, double *y){
    *y = expm1(*x);
    return 1;
}

double expm1(double);

int expm(int argc, void *argv[]){  
   return expm_natural((double *) argv[0], (double *) argv[1]);
}

 
#ifdef MAIN
main()
{
    extern long nfunc;
    double i, err;
    long nfunc_out;
    int status;

    double p1, p2, chi0, eta, mu, eps;

    p1 = 2.5;
    p2 = 2.5;
    chi0 = 1.0;
    eta = 0.01;
    mu = 1000.;
    eps = 1.e-6;

    printf("One-part integration:\n");
    status= gamma_integral_natural(&p1, &p2, &chi0, &eta, &mu, &eps, &i, &err, &nfunc_out);
    printf(" I_gam= %lg, err= %lg, N= %ld\n", i, err, nfunc_out);

    // The two-part integral seems less accurate--requires a smaller eps (and therefore takes longer)

    // One remedy might be to do the (0,inf) integration and subtract from it (0,chi) to get second part
    // This would mean three integrations instead of two, however, but it is indeed faster

    printf("\n\nTwo-part integration:\n");
    p1 = 2.5001;

    status =  gamma_integral_natural(&p1, &p2, &chi0, &eta, &mu, &eps, &i, &err, &nfunc_out);
    printf(" I_gam= %lg, err= %lg, N= %ld\n", i, err, nfunc_out);
}
#endif
