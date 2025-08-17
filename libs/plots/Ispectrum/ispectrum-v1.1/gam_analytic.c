
//-----------------------------------------------------------------------------
// Compute 8/pi * Integral_0^Infinity sin(eta*x0/2)^2 * sin(mu*x0/2)^2 / x^p dx
// valid for 1 < p < 5, the case p=3 requiring special treatment
//
// Call with double precision arguments for a double precision result!
//
// Don't forget to multiply result by U1 to obtain gamexact with proper scaling!
//-----------------------------------------------------------------------------

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "gam_analytic.h"

#ifndef M_PI
   #define M_PI 3.1415926535897932
#endif

// GNU implementation of gamma function
double tgamma(double x);

double gam_analytic(double p, double eta, double mu){

  double pm, r, c, prev_result, result;
  double m, e;
  double atanhr;
  double c0, Cp1;
  double coef,term;
  int i,k,fact;
  double pi = M_PI;
  double eps = 1.e-15;
  int max_iter = 100;
  double r_thresh = 1.E-6; // Threshold for applying Talyor series expansion, determined via experiment

  if ((p <= 1.0) || (p >= 5.0)){
	    printf("gam_analytic is not defined for this value of p: %f\n",p);
        exit(1);
  }

  pm = p - 1.0;
  
  if (p == 3.0){
	
	if (eta == mu){

		c = pow(fabs(2.0/eta),1.0-p);
		result = 0.25 * pow(eta,2.0) * log(2.0);
        
        // we need to return here!
        
	} else {
		if (eta > mu){
			e = eta;
			m = mu;
			c = pow(fabs(2.0/eta),1.0-p);
		} else {
			m = eta;
			e = mu;
			c = pow(fabs(2.0/mu),1.0-p);
		}
		r = m/e;
		atanhr = 0.5*log((1.0+r)/(1.0-r)); // formula valid for r<1

		Cp1 = 1.0/8.0;
		result = -pow(fabs(eta),2.0) * log(fabs(eta)) - pow(fabs(mu),2.0) * log(fabs(mu)) 
			+ 0.50*pow(fabs(eta-mu),2.0) * log(fabs(eta-mu)) + 0.50*pow(fabs(eta+mu),2.0) * log(fabs(eta+mu));
  		result = result * Cp1;

		// ------------------------------------------------------------------------------------
		// The Taylor expansion below mitigates roundoff errors when eta/mu or mu/eta are large
		// ------------------------------------------------------------------------------------
		if (pow(r,pm) < r_thresh){

			// Taylor series in r about r=0 to order 2
             result = c * (3.0/4.0 - log(r)/2.0 ) * pow(r,2.0);

			// For purposes of error control, we add subsequent higher order terms until convergence

			for (i=2;i<=max_iter;i++){
				k = 2*i;
				coef = 4.0*i*(i-1)*(2*i-1);
				term = c * pow(r,k) / coef;
				//printf("i %d term %e\n",i,term);
				result = result - term;
				if (fabs(term) < eps*fabs(result)){ break; }
			}
			if (fabs(term) > eps*fabs(result)){
				printf("Warning, maximum number of iterations exceeded in gamma_analytic\n");
			}
			//return result*8/pi;
		}


        // previous line commented out 9 Oct 2016
		//return result*8/pi;

	} 
    
     // line added out 9 Oct 2016
    return result*8/pi;
    
 }

  // ----------------
  // p not equal to 3
  // ----------------
  
  c0 = 1.0/(sqrt(pi)*tgamma(p/2.0)*(p-1.0)*pow(2.0,p-1.0));
  
  Cp1 = c0 *  tgamma((3.0-p)/2.0); 
  
  result = 2.0*pow(fabs(eta),pm) + 2.0*pow(fabs(mu),pm) - pow(fabs(eta+mu),pm) - pow(fabs(eta-mu),pm);

  result = result * Cp1;
  //printf("c0 %f, cp1 %f, pm %f result %f\n",c0,Cp1,pm,result);

   // ------------------------------------------------------------------------------------
   // The Taylor expansion below mitigates roundoff errors when eta/mu or mu/eta are large
   // ------------------------------------------------------------------------------------
   
	if (mu > eta){
		r = eta/mu;
		c = pow(fabs(mu),pm);
	} else {
		r = mu/eta;
		c = pow(fabs(eta),pm);
	}

	if (pow(r,pm) < r_thresh){
		// gamma(r>0; mu) = abs(mu)^pm * [2 r^pm - (p-2)(p-1)*r^2*2/2!
		//                             - (p-4)(p-3)(p-2)(p-1)*r^4*2/4!
		//                   - (p-6)(p-5)(p-4)(p-3)(p-2)(p-1)*r^6*2/6!
		//         - (p-8)(p-7)(p-6)(p-5)(p-4)(p-3)(p-2)(p-1)*r^8*2/8! - ...]

		// Taylor series in r about r=0, add terms until convergence reached
		prev_result = result;
		result = 2.0*pow(fabs(r),pm);
		coef = 1.0;
		fact = 1;
		for (i=1;i<=max_iter;i++){
			k = 2*i;
			fact = fact*(k-1)*k;
			coef = coef * (pm-k+1.0)*(pm-k+2.0);
			term = coef * pow(r,k) * 2.0/fact;
			//printf("i %d term %e\n",i,term);
			result = result - term;
			if (fabs(term) < eps*fabs(result)){ 
				result = result * Cp1 * c;
				return result;
			}
		}
		printf("Warning, maximum number of iterations exceeded in gam_analytic\n");
		result = result * Cp1 * c;		
	}
   
 return result;
}

#ifdef MAIN
int main(int argc, char *argv[]){
double p;
double eta;
double mu;
double gam;
if (argc != 4){
   printf("Usage: gam_analytic p eta mu\n");
   exit(1);
}
p   = atof(argv[1]);
eta = atof(argv[2]);
mu  = atof(argv[3]);
//printf("%s %s %s\n",argv[1],argv[2],argv[3]);
//printf("p: %f, eta: %f, mu: %f\n",p,eta,mu);
gam = gam_analytic(p,eta,mu);
printf("%20.16f\n",gam);
return 0;
}
#endif
