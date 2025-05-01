
//============================================================================================
// Computes the normalized intensity SDF for a plane wave that traverses a phase screen
// Written by Charles Carrano, Boston College (charles.carrano@bc.edu)
//
// Execute with no command line options for a description of the program and its usage. 
// For additional details, see
//
// Carrano, C. S., and C. L. Rino (2016), A theory of scintillation for two-component\n");
//    power law irregularity spectra: Overview and numerical results, Radio Sci., \n");
//    51, 789-813, doi:10.1002/2015RS005903.\n\n");
//    
// To compile this program under Linux:
//   gcc -c ispectrum.c -lm
//   gcc -o ispectrum ispectrum.o gam_quadrature.o gamma_integral.o gam_analytic.o intde1.o -lm
//
// To do:
// [ ] Prevent spurious I(mu) for very large mu with inner scale is used
// [ ] Separate gamma_integral.c routines that are specific to IDL (rename to avoid awkward "_natural")
// [ ] Place under version control CVS repository
// [ ] Try applying log change of variable when integrating I(mu) to get S4 (just when U large?)
//============================================================================================

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "gam_quadrature.h"
#include "gam_analytic.h"
#include "gamma_integral.h"
#include "intde1.h"

#define DEBUG 1

#ifndef M_PI
   #define M_PI 3.1415926535897932
#endif

static double U1;
static double U2;
static double p1;
static double p2;
static double mu0;
static double mu_outer;
static double mu_inner;

static double mu;

static long nfunc=0;

static int debug=0;

static unsigned int numv=0;

//#define NMAX 65536
#define NMAX 65535

static double muv[NMAX], intv[NMAX];
static long nfuncv[NMAX];

/* holds address of array of which the sorted index order needs to be found */
static double *base_arr;
 
// Compare function for sorting arrays. Value passed to this
// function by qsort are the idx array elements.
static int compar (const void *a, const void *b){
  int aa = *((int *) a), bb = *((int *) b);
  if (base_arr[aa] < base_arr[bb])   return -1;
  if (base_arr[aa] == base_arr[bb])  return 0;
  if (base_arr[aa] > base_arr[bb])   return 1;
  
  return 0;
}
       
// x = eta
double gam(double x){ // move this to gamma_integral.c ?
   double gam;
   double err_gam;
   long nfunc_gam;

   gam = gam_new(x, mu, 
               U1, U2, p1, p2, mu0, mu_outer, mu_inner, debug,
               &err_gam, &nfunc_gam);
   //printf("mu: %f, eta: %f, mu-eta: %f, gam: %f\n",mu,x,mu-x,gam);
   return gam;
}

// x = eta
double spectrum_integrand(double x){
   nfunc++;
   return expm1(-gam(x))*cos(x*mu);
}

// transformation is eta-> x/mu
double spectrum_integrand_quadrature_transformed(double x){
   nfunc++;
   return expm1(-gam(x/mu))* cos(x); 
}

// x = mu
double spectrum_old(double x){
   double eps = 1.e-8;
   double Int;
   double err;
   mu = x;
   nfunc = 0;
   intdeo(spectrum_integrand,0.0,mu,eps, &Int,&err);
   //printf("nfunc: %ld\n",nfunc);
   return 2.0*Int;
}

// x = mu
double spectrum(double x){
   double eps = 1.e-8;
   double Int;
   double err;    
   double Int1, Int2;
   double err1, err2;
   double hack_thresh;
   
   mu = x; 
   
   nfunc = 0;
   
   hack_thresh = mu_outer/5.0;
   //printf("mu_outer: %f, hack_thresh: %f\n",mu_outer,hack_thresh);

   if (mu < 1.0){
       
       if (mu < hack_thresh){
         // a hack to prevent challenging evaluation of I(mu) for mu << mu_outer
         Int1 = U1 * pow(mu_outer,-p1) * pow(mu,4.0); // asymptote for outer scale region
         err = 0.0;
         err1 = 0.0;
         Int2 = 0.0;
         err2 = 0.0;         
         //printf("Invoking quick approx for mu: %e, I: %e\n",mu,Int1);
         
       } else {
                
                // use eta*mu->y, upper limit eta=mu0 -> y(eta|eta=mu0) = mu0*m0 = mu0^2 

                intde(spectrum_integrand_quadrature_transformed,0.0,mu0*mu0,eps, &Int1,&err1);
                Int1 = 2.0*Int1/mu;
                //Int1 = 2.D0*qintde('spectrum_integrand_quadrature_transformed',0.D0,mu0^2.D0,epsrel,err1,nfunc=nfunc1)/mu
                
                intdeo(spectrum_integrand_quadrature_transformed,mu0*mu0,1.0,eps, &Int2,&err2);
                Int2 = 2.0*Int2/mu;
                //Int2 = 2.D0*qintdeo('spectrum_integrand_quadrature_transformed',mu0^2.D0,1.D0, epsrel,err2,nfunc=nfunc2)/mu
               
       }

       Int = Int1 + Int2;
       err = sqrt(err1*err1 + err2*err2);
       //nfunc = nfunc1 + nfunc2;
         
   } else {
 
      intdeo(spectrum_integrand,0.0,mu,eps, &Int,&err);
      Int = 2.0*Int;
      //Int = 2.D0*qintdeo('spectrum_integrand_quadrature',0.D0,mu,epsrel,err,nfunc=nfunc);
   
   }

   //printf("nfunc: %ld\n",nfunc);
   if (Int < 0.0){ Int = 0.0; }
   
   if (numv < NMAX){
      muv[numv] = mu;
      intv[numv] = Int;
      nfuncv[numv] = nfunc;
      #ifdef DEBUG
         printf("%f %20.15e %ld\n",mu,Int,nfunc);
      #endif
      numv++;
   }
   return Int;
}

// transformation x->sqrt(2t) makes integrand oscillate with wavenumber 2
// --this helps to accurately resolve Fresnel nulls in weak scatter
double spectrum_transformed_weak(double t){
    double phi;
    phi = sqrt(2.0*t);
    return spectrum(phi)/phi;
}

// transformation x->(1-t)/t maps [0,inf] to [0,1]
double spectrum_transformed(double t){
    return spectrum((1.0-t)/t)/pow(t,2.0);
}

// spectral density function of phase
double pspectrum(double x){
   if (mu_outer>0 && x<mu_outer){ return U1*pow(mu_outer,-p1);}
   if (x<mu0){ return U1*pow(x,-p1); }
   if (x<mu_inner){ return U1*pow(mu0,p2-p1)*pow(x,-p2); } 
   if (mu_inner>0){ return 0.0; } else { return U1*pow(mu0,p2-p1)*pow(x,-p2); }
}

double spectrum_weak_smoothed(double x){
   double x4, filterfunction;
   x4 = pow(x,4.0);
   filterfunction = x4/(4.0+2.0*x4);
   return (4.0*filterfunction)*pspectrum(x);
}

// this function oscillates about zero with wavenumber 2
double integrand_transformed_weak(double t){
    double x;
    x = sqrt(2.0*t);
    return (spectrum(x) - spectrum_weak_smoothed(x))/x;
}

double spectrum_transformed_weak_smoothed(double t){
    double x;
    x = sqrt(2.0*t);
    return spectrum_weak_smoothed(x)/x;
}


// Compute the rms phase fluctuation
//    sigma_phi = int_{-inf}^{inf} I_phi(k) dk/(2pi)
//              = int_{-inf}^{inf} P(mu) dmu/(2pi)
//              = 2 int_{0}^{inf} P(mu) dmu/(2pi)
//    (We return -1.0 if sigma phi does not exist)
double calc_sigphi(void){

  double term1,term2;
  double ans; 
  double pi = M_PI;
  double outer_term, inner_term;

  if ((mu_outer == 0) && (p1>=1)) { return -1.0; }
  if ((mu_inner == 0) && (p2<=1)) { return -1.0; }

  if (mu_outer > 0) { outer_term = pow(mu_outer,1.0-p1)*p1; } else { outer_term = 0.0; }
  if (mu_inner > 0) { inner_term = pow(mu_inner,1.0-p2);    } else { inner_term = 0.0; }
  
  // requires p1<1 if have_outer = 0
  term1 = U1 * 2.0/(p1-1.0) * (outer_term - pow(mu0,1.0-p1));
  
  // requires p2>1 if have_inner = 0
  term2 = U2 * 2.0/(p2-1.0) * (pow(mu0,1.0-p2) - inner_term );

  printf("outer_term: %e, inner_term: %e, term1: %e, term2: %e\n",outer_term,inner_term,term1,term2);
  ans = term1 + term2;
    
  return sqrt(ans/(2.0*pi));
}

// Compute the rms density fluctuation. The density variance is given by 
//   var_N = int{-inf}^{inf} int{-inf}^{inf} I_N(kx,kz) dkx dkz/(2*pi)^2 
//         = int{-inf}^{inf} int{-inf}^{inf} P(mux,muz) dmux dmuz/(2*pi)^2 / (re^2 lambda^2 * L  * rhof)
//         = int{0}^{int} P(mu) mu dmu/(2*pi) / (re^2 lambda^2 * L  * rhof)
// The RMS density fluctuation is sigma_N = sqrt(var_N)   
// This subroutine actually computes the quantity sigma_N * re*lambda*L^(1/2)rhof^(1/2),
// In other words, you must divide the output of this subroutine by 
// re*lambda*L^(1/2)rhof^(1/2) to obtain the value of sigma_N.
// (We return -1.0 if the RMS density fluctuation does not exist)
double calc_sigNfac(void){

  double term1,term2;
  double ans;
  double pi = M_PI;
  double outer_term, inner_term;
  
  if ((mu_outer == 0) && (p1>=2)) { return -1.0; }
  if ((mu_inner == 0) && (p2<=2)) { return -1.0; }
  
  if (mu_outer > 0) { outer_term = pow(mu_outer,2.0-p1)*p1/2.0; } else { outer_term = 0.0; }
  if (mu_inner > 0) { inner_term = pow(mu_inner,2.0-p2);        } else { inner_term = 0.0; }
  
   if (p1==2.0){
      term1 = U1 * ( 0.5 + log(mu0/mu_outer) );
      //printf("Nfac term1a: %f\n",term1);
   }else{
      // requires p1<2 if have_outer = 0
      term1 = U1 / (p1-2.0) * ( outer_term - pow(mu0,2.0-p1) );
      //printf("Nfac term1b: %f\n",term1);
   }
  
  if (p2==2.0){
     term2 = U2 * log(mu_inner/mu0);
     //printf("Nfac term2a: %f\n",term2);
  }else{
     // requires p2>2 if have_inner = 0
     term2 = U2 / (p2-2.0) * ( pow(mu0,2.0-p2) - inner_term );
     //printf("Nfac term2b: %f\n",term2);
  }

  //printf("Nfac term1: %f, term2: %f\n",term1,term2);
  ans = term1 + term2;
    
  return sqrt(ans/(2.0*pi));
}

// ============
// Main program
// ============

int main(int argc, char *argv[]){
    
    double result;
    double Ustar;
    double mui;

    FILE *fp,*fl;
    
    if (argc < 7 || argc> 10){
    printf("Program Description: \n");
    printf("\nComputes the normalized intensity SDF for a plane wave that traverses a\n");
    printf("phase screen\n");
    printf("\nI(mu) = int_{-inf}^{inf} exp[-gamma(eta,mu)] exp(-i eta mu) d eta\n\n");
    printf("The normalized screen SDF is specified as a piecewise power law\n\n"); 
 
    printf("           {   muo^(-p1),          if  0<= mu <= muo,\n");
    printf("P(mu) = U1 {    mu^(-p1),          if muo < mu <= mub,\n");
    printf("           { mub^(p2-p1) mu^(-p2), if mub < mu <= mui,\n");

    printf("\nwhere U1 = Cp rhof^(p1-1) and \n\n");
    printf("    Cp   = phase spectral strength\n");
    printf("    rhof = sqrt(z/wavk) is the Fresnel scale\n");
    printf("    z    = distance past the screen\n");
    printf("    wavk = signal wavenumber\n");
    printf("    k    = transverse wavenumber\n");
    printf("    mu   = rhof*k  normalized transverse wavenumber\n");
    printf("    muo  = rhof*ko normalized outer scale wavenumber\n");
    printf("    mub  = rhof*kb normalized break scale wavenumber\n");
    printf("    mui  = rhof*ki normalized inner scale wavenumber\n");

    printf("\nWe assume that 0 < muo < mub < mui and also that muo << muf << mui, where\n");
    printf("muf = 2*pi is the normalized wavenumber corresponding to the Fresnel scale.\n");
    
    printf("\nThe structure interaction function gamma(eta,mu) is given by\n");

    printf("\ngamma(eta,mu) = \n");
    printf("   16 U1 int_{  0,muo} muo^(-p1)   sin^2(chi eta/2) sin^2(chi mu/2) d chi/(2 pi)\n");
    printf(" + 16 U1 int_{muo,mub}             sin^2(chi eta/2) sin^2(chi mu/2) d chi/(2 pi)\n");
    printf(" + 16 U1 int_{muo,mui} mub^(p2-p1) sin^2(chi eta/2) sin^2(chi mu/2) d chi/(2 pi)\n\n");
    printf("Universal scattering strength U equals U1 if mub>=1 and U1 mub^(p2-p1) otherwise\n");
    
    printf("\nThis program fully supports the limiting cases muo->0 and/or mui->infinity\n");
    
    printf("\n--------------------------------------------------------------------------------\n");
    
    printf("\nProgram Usage: ispectrum U p1 p2 mub muo mui [mu] | ([mu_min] [mu_max] [mu_num])\n");
    printf("\nCalling options: \n");
    printf(" 1: ispectrum U p1 p2 mub muo mui\n");
    printf(" 2: ispectrum U p1 p2 mub muo mui mu\n");
    printf(" 3: ispectrum U p1 p2 mub muo mui mu_min mu_max mu_num\n");
    printf("\nNotes:\n");
    printf(" * Setting muo and/or mui to zero omits them from the model\n");
    printf(" * Option1 uses adaptive quadrature to compute S4\n");
    printf(" * Option2 computes I(mu) at the specified value of mu\n");
    printf(" * Option3 computes I(mu) at log spaced mu values between mu_min and mu_max\n");
    printf(" * The intensity SDF I(mu) is written to ispectrum.dat\n");
    printf(" * Parameters and moments are written to ispectrum.log\n");
    printf("\nWritten by Charles Carrano, Boston College (charles.carrano@bc.edu)\n");
    printf("\nFor additional details, see\n");
    printf("\n   Carrano, C. S., and C. L. Rino (2016), A theory of scintillation for\n");
    printf("        two-component power law irregularity spectra: Overview and numerical\n");
    printf("        results, Radio Sci., 51, 789-813, doi:10.1002/2015RS005903.\n\n");
    
    exit(0);
    }
   
    fp=fopen("ispectrum.dat","w");
    fl=fopen("ispectrum.log","w");

    Ustar    = atof(argv[1]);
    p1       = atof(argv[2]);
    p2       = atof(argv[3]);
    mu0      = atof(argv[4]);
    mu_outer = atof(argv[5]);
    mu_inner = atof(argv[6]);
    
    
    if (mu0 == 0) mu0 = 1.0; // so we can compute U2 = U1 mu0^(p2-p1), even if mu0 was set to zero

    if (mu0 >= 1){
       U1 = Ustar;
       U2 = U1 * pow(mu0,p2-p1);
    } else {
       U2 = Ustar;
       U1 = U2 * pow(mu0,p1-p2);
    }

    if (argc == 8){
       // compute I(mu) for a single mu
       mui  = atof(argv[7]);
       result = spectrum(mui);
       #ifdef DEBUG
          printf("%f %20.15e %ld\n",mui,result,nfunc);
	   #endif
       fprintf(fp,"%f %20.15e %ld\n",mui,result,nfunc);
       return 0;
    }
    
    if (argc == 10){
           // compute I(mu) for a range of mu values
           double mu_min;// = 1.e-6;
           double mu_max;// = 1.e6;
           double ratio;
           int num;// = 50;
           int i;
           mu_min  = atof(argv[7]);
           mu_max  = atof(argv[8]);
           num     = atoi(argv[9]);

           mui = mu_min;
           ratio = pow(10.0,log10(mu_max/mu_min)/((double)num-1));
           //printf("ratio: %f\n",ratio); exit(1);
           result = spectrum(mui);
		   // lines below for experimentation only
           // result = pspectrum(mui);
           // result = spectrum_weak_smoothed(mui);
           // result = integrand_transformed_weak(mui); // commented out 8 Mar 2017

		   #ifdef DEBUG
              printf("%f %20.15e %ld\n",mui,result,nfunc);
		   #endif
           fprintf(fp,"%f %20.15e %ld\n",mui,result,nfunc);
           numv = 0;
           for (i=1;i<num;i++){
               mui = mui * ratio;
               // lines below for experimentation only   
               // result = pspectrum(mui);
               // result = spectrum_weak_smoothed(mui);
               // result = integrand_transformed_weak(mui); // commented out 8 Mar 2017
              result = spectrum(mui);
			   #ifdef DEBUG
               //printf("%f %20.15e %ld\n",mui,result,nfunc);
               fprintf(fp,"%e %20.15e %ld\n",mui,result,nfunc);
			   #endif
           }
           fclose(fp);
		   #ifdef DEBUG
              printf("number of wavenumbers evaluated: %d\n",numv+1);
		   #endif
           return 0;
    }
    
    
    if (argc == 7){
       // compute S4 using adaptive quadrature
       double ans,err;
       double S4;
       double sigphi,sigNfac;
       double eps;
       double pi = M_PI;
       int i;
       int idx[NMAX];

       double thresh_weak; 

	   // use a threshold for weak scatter that depends on p2 in an attempt
	   // to maintain a smooth result as p2 changes (crosses the threshold)
	   // this works well for one-component spectra, but probably should 
	   // depend also on p1 and mub for two component spectra

       thresh_weak = -0.25*p2 + 0.75; // gives 0.5 for p2=2, 0.25 for p2=4

       #if 1 
       //if (Ustar<0.5){ // works great for p=2
       //if (Ustar<0.8){ -- tried for p=4, gives discontinuity at U=0.8
       //if (Ustar<0.35){ 
       if (Ustar<thresh_weak){ 

          double ans1, ans2;
          eps  = 0.0001;

          // option below is very efficient when scatter is weak 
          intdeo(integrand_transformed_weak, 0.0, 2.0, eps, &ans1, &err);
          intdei(spectrum_transformed_weak_smoothed, 0.0, eps, &ans2, &err);
          //printf("ans1: %f, ans2: %f\n",ans1,ans2);

          ans = ans1 + ans2;
          
       } else {

	      // eps = 0.01;
	      if (p1 == p2){ 
				eps  =  0.0001; 
			} else { 
				// two-component spectra are CPU intensive, so compute with less precision
				eps = 0.001; 
			}
          intde(spectrum_transformed, 0.0, 1.0, eps, &ans, &err);

       }
       #endif
       
       S4 = sqrt(2.0*ans/(2.0*pi));

       #if 0
          // four ways to compute S4 -- for experimentation only, do not use!
       
          int numvw,numv0,numv1,numv2,numv3;
          double ans1,ans2;
		  double S4w, S40,S41,S42,S43;
       
          eps = 0.01;
          //eps = 0.0001;
   
          // option below is very efficient when scatter is weak (gives nan when p=3 and Ustar~1)
          numv=0;
          intdeo(integrand_transformed_weak, 0.0, 2.0, eps, &ans1, &err);
          intdei(spectrum_transformed_weak_smoothed, 0.0, eps, &ans2, &err);
          ans = ans1 + ans2;
          S4w = sqrt(2.0*ans/(2.0*pi));
          numvw=numv;
       
          // option below gives pretty spectrum but doesn't give accurate S4 when spectrum is shallow (p<3)
          numv=0;
          intdeo(spectrum_transformed_weak, 0.0, 2.0, eps, &ans, &err);
          S40 = sqrt(2.0*ans/(2.0*pi));
          numv0=numv;
       
          // this option saves some function evaluations in weak scatter
          numv=0;
          intdei(spectrum_transformed_weak, 0.0, eps, &ans, &err);
          S41 = sqrt(2.0*ans/(2.0*pi));
          numv1=numv;

          numv=0;
          intde(spectrum_transformed, 0.0, 1.0, eps, &ans, &err);
          S42 = sqrt(2.0*ans/(2.0*pi));       
          numv2=numv;
       
          numv=0;
          intdei(spectrum, 0.0, eps, &ans, &err);
          S43 = sqrt(2.0*ans/(2.0*pi));
          numv3=numv;
       
          printf("ans1 %f, ans2: %f\n",ans1,ans2); 
          printf("intdew S4: %f, numv: %d\n",S4w,numvw); 
          printf("intdeo S4: %f, numv: %d\n",S40,numv0); 
          printf("intdei S4: %f, numv: %d\n",S41,numv1); 
          printf("intde  S4: %f, numv: %d\n",S42,numv2); 
          printf("intdei S4: %f, numv: %d\n",S43,numv3); 
       
       #endif
       
	   // Sort the intensity results by in increasing wavenumber

       /* initialize initial index permutation of unmodified `arr'
       */
       for (i = 0; i < numv; i++)
       {
           idx[i] = i;
        }
    
       /* Assign the address of out original array to the static global
        * pointer, this will be used by the compare function to index 
        * into the original array using `idx' values
       */
       base_arr = muv;
       
       // sort by wavenumber
       qsort (idx, numv, sizeof (int), compar);

       for (i=1;i<numv;i++){
          //printf("%f %20.15e %ld\n",muv[i],intv[i],nfunc);
          //we need to sort these first
          fprintf(fp,"%e %20.15e %ld\n",muv[idx[i]],intv[idx[i]],nfuncv[idx[i]]);
        }
        fclose(fp);

        //sigphi  =  calc_sigphi(Ustar,p1,p2,mu0,mu_outer,mu_inner);
        //sigNfac = calc_sigNfac(Ustar,p1,p2,mu0,mu_outer,mu_inner);

        sigphi  = calc_sigphi();
        sigNfac = calc_sigNfac();

		#ifdef DEBUG
           printf("number of wavenumbers evaluated: %d\n",numv+1);
           //printf("Sigma-phi: %f %f\n",sigphi,calc_sigphi_old(Ustar,p1,p2,mu0,mu_outer,mu_inner));
           printf("Sigma-phi: %f\n",sigphi);
           printf("SigmaNfac: %f\n",sigNfac);
           printf("S4: %f\n",S4);
        #endif

        fprintf(fl,"#   Ustar         U1           U2         p1       p2      mu_break     mu_outer     mu_inner      S4        sigphi   sigNfac num\n");
        fprintf(fl,"%e %e %e %f %f %e %e %e %f %e %e %d\n",Ustar,U1,U2,p1,p2,mu0,mu_outer,mu_inner,S4,sigphi,sigNfac,numv+1);
        fclose(fl);
        //printf("U1: %e, U2: %e\n",U1,U2);
        return 0;
    }
    return 0;
}
