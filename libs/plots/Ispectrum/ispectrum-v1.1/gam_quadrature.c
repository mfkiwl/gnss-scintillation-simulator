// Compile like this:
//gcc -c gam_quadrature.c -lm
//gcc -o gam_quadrature gam_quadrature.o gamma_integral.o gam_analytic.o -lm

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "gamma_integral.h"
#include "gam_analytic.h"

#ifndef M_PI
   #define M_PI 3.1415926535897932
#endif

double gam_outer_r(double r, double y0){
   double y2;
   double r2,r4;
   double thresh1,thresh2;
   double r2m1,rp1,rm1;  
   double numer,denom;
   double result;
   
   y2 = 2.0*y0;
     
   if (r == 1.0){
       result = (6.0 * y2 - 8.0 * sin(y2) + sin(2.0*y2)) / 32.0;
       return result;
   }
   
   thresh1 = 5.e-1;
   thresh2 = 1.0/thresh1;
   
   r2 = pow(r,2.0);

   if (y0 < thresh1){ 
       // 10th order Taylor series in y0 about zero
       r4 = pow(r2,2.0);
       result = (r2 * pow(y0,5.0))/5.0 - r2 * (1.0 + r2) * pow(y0,7.0) / 21.0 + 
                 r2 * (2.0 + 5.0 * r2 + 2.0 * r4) * pow(y0,9.0) / 405.0;
       return result;
   }
  
    // intermediate values of r
    r2m1 = r2 - 1.0;
    rp1  = r + 1.0;
    rm1  = r - 1.0;

    numer = -2.0 * r2m1 * (-y2 * r + r * sin(y2) + sin(y2 * r)) + 
             rm1 * r * sin(y2 * rp1) + r * rp1 * sin(y2 * rm1);
    
    denom =  16.0 * r * r2m1;
    
    result = numer / denom; 
    return result;
}

double gam_outer(double eta, double mu, double mu_outer){
   double r,c,result;
   double y0;
   
   //eta = eta_in;
   //mu = mu_in;
     
   if (eta <= mu){
      r = eta/mu; // < 1
      c = mu/2.0;
   } else {
      r = mu/eta; // < 1
      c = eta/2.0;
   }
     
   y0 = c*mu_outer ; // transformed outer scale

   result = gam_outer_r(r,y0)/c;
   
   return result;
}
    
    
double gam_new(double eta, double mu, 
               double U1, double U2, double p1, double p2, double mu0, double mu_outer, double mu_inner,
               int debug, double *err_, long *nfunc_){

   double eps;
   double pi = M_PI;
   double thresh1, thresh2, thresh3;
   long iter1,iter2,iter3;
   long nfunc3;
   double ans0a,ans0b,ans3q, result;
   double ans0q,ans1q,ans2q;
   double a;
   double err;
   //long nfunc;

   *err_ = 0.0;
   *nfunc_ = 0;
                
   eps =  1.e-10; //(set too low and the outer integration converges more slowly)

   thresh1 = 1.e-1; 
   thresh2 = 1.e1; 
  
   thresh3 = 1.e5;
    
   if (p1 == p2){ // *** 1-component spectrum ***
         
         iter1=0;
         iter2=0;
         iter3=0;
         if ((eta / mu < thresh1 || eta / mu > thresh2 || mu > thresh3) && p1 > 1.0){
                if (mu_outer == 0.0){
                    ans0a = 0.0;
                    ans0b = 0.0;
                } else {
                    a = 0.0;
                    (void) gamma_finite_right_natural(&p1, &a, &mu_outer, &eta, &mu, &eps, &ans0a, &err, &iter1);
                    ans0a = ans0a* 8.0/pi ;
                    //printf("ans0a: %e, iter1: %ld\n",ans0a,iter1);
                    //ans0a = gamma_finite_right(p1, 0.0, mu_outer, eta, mu, eps, nfunc=iter1)* 8.0/pi ;   
                   
                    ans0b = gam_outer(eta, mu, mu_outer)* 8.0/pi; 
                    //printf("ans0b: %e\n",ans0b);

                    // note: gam_outer(eta,mu,mu_outer)*8.0/pi*U1*mu_outer^(-p1) = gamma_series_outer(eta,mu)
                }
                
                if (mu_inner == 0.0){
                    ans3q = 0.0;
                    nfunc3 = 0;
                } else {
                    (void) gamma_infinite_right_natural(&p2, &mu_inner, &eta, &mu, &eps, &ans3q, &err, &iter3);
                    ans3q = ans3q * 8.0/pi;
                    //printf("ans3q: %e\n",ans3q);
                    //ans3q = gamma_infinite_right(p2, mu_inner, eta, mu, eps, nfunc=iter3)* 8.0/pi;
                }
                   
                   
                if (mu_outer > 0.0){                   
                    result = U1*(gam_analytic(p1,eta,mu) - ans0a) + U1*pow(mu_outer,-p1)*ans0b - U2*ans3q;                    
                } else {              
                    result = U1*gam_analytic(p1,eta,mu) - U2*ans3q;                    
                }
                
                if (debug > 0){ printf("a %e %e %ld %ld %ld %e\n",eta,mu,iter1,iter2,iter3,result); }

        } else { // ; eta ~ mu
 
                if (mu_outer == 0){
                    ans0a = 0.0;
                    ans0b = 0.0;
                } else {
                    a = 0.0;
                    (void) gamma_finite_natural(&p1, &a, &mu_outer, &eta, &mu, &eps, &ans0a, &err, &iter1);
                    ans0a = ans0a * 8.0/pi;
                    //ans0a = gamma_finite(p1, 0.0, mu_outer, eta, mu, eps, nfunc=iter1)* 8.0/pi;                   
                    ans0b = gam_outer(eta, mu, mu_outer)* 8.0/pi;
                }    
                    
                if (mu_outer > 0.0){
                    if (mu_inner > 0.0){
                       (void) gamma_finite_natural(&p2, &mu_outer, &mu_inner, &eta, &mu, &eps, &result, &err, &iter1);
                       result = result*U1* 8.0/pi + U1*pow(mu_outer,-p1)*ans0b;
                       //result = U1*gamma_finite(p2, mu_outer, mu_inner, eta, mu, eps, nfunc=iter1)* 8.0/pi + U1*pow(mu_outer,-p1)*ans0b;
                       
                    } else {
                       (void) gamma_infinite_natural(&p2, &mu_outer, &eta, &mu, &eps, &result, &err, &iter2);
                       result = result* U1* 8.0/pi + U1*pow(mu_outer,-p1)*ans0b;
                       //result = U1*gamma_infinite(p2, mu_outer, eta, mu, eps, nfunc=iter2)* 8.0/pi + U1*pow(mu_outer,-p1)*ans0b;
                    }                

                } else {
                
                    if (mu_inner > 0.0){
                        a = 0.0;
                       (void) gamma_finite_natural(&p2, &a, &mu_inner, &eta, &mu, &eps, &result, &err, &iter2); 
                       result = result*U1* 8.0/pi;
                       //result = U1*gamma_finite(p2, 0.0, mu_inner, eta, mu, eps, nfunc=iter2)* 8.0/pi;
                    } else {
                       result = U1*gam_analytic(p1,eta,mu);
                       //if (isnan(result)){
                       //    printf("result is nan: p1 :%f, eta: %f, mu: %f\n",p1,eta,mu);   
                       //}
                    }
                    
                }

                if (debug > 0){ printf("b %e %e %ld %ld %ld %e\n",eta,mu,iter1,iter2,iter3,result); }
                //if result[i] ne result[i] then stop

        }
      
      *err_ = 0.0;
      *nfunc_ = iter1 + iter2 + iter3;
                
      //}
      if (result < 0.0) { result = 0.0; }
      return result;
      
   } else { // ***2-component spectrum ***
             
            iter1=0;
            iter2=0;
            iter3=0;
            
            if (mu_outer == 0.0){
                ans0q = 0.0;
                a = 0.0;
            } else {
                a = 0.0;
                (void) gamma_finite_natural(&a, &a, &mu_outer, &eta, &mu, &eps, &ans0q, &err, &iter1); 
                ans0q = ans0q* 8.0/pi;
                //ans0q = gamma_finite(0.0, 0.0, mu_outer, eta, mu, eps, nfunc=iter1)* 8.0/pi;
                a = mu_outer;   
            }
              
            if ((eta / mu < thresh1 || eta / mu > thresh2 || mu > thresh3) && p1 > 1.0){
                                
                if (mu_inner == 0.0){
                   ans3q = 0.0;
                   iter3 = 0;
                } else {
                   (void) gamma_infinite_right_natural(&p2, &mu_inner, &eta, &mu, &eps, &ans3q, &err, &iter3);
                   ans3q = ans3q*8.0/pi;
                   //ans3q = gamma_infinite_right(p2, mu_inner, eta, mu, eps, nfunc=iter3)* 8.0/pi;
                }
                
                (void) gamma_finite_right_natural(&p1, &a, &mu0, &eta, &mu, &eps, &ans1q, &err, &iter1);
                ans1q = ans1q* 8.0/pi;
                //ans1q = gamma_finite_right(p1, a, mu0, eta, mu, eps, nfunc=iter1)* 8.0/pi; 
               
                (void) gamma_infinite_right_natural(&p2, &mu0, &eta, &mu, &eps, &ans2q, &err, &iter2);
                ans2q = ans2q* 8.0/pi ;
                //ans2q = gamma_infinite_right(p2, mu0, eta, mu, eps, nfunc=iter2)* 8.0/pi ;
                
                if (mu_outer > 0.0){
                   result = U1*pow(mu_outer,-p1)*ans0q + U1*ans1q + U2*ans2q - U2*ans3q;
                } else {
                   result =  U1*ans1q + U2*ans2q - U2*ans3q;
                }
                *err_ = 0.0;
                *nfunc_ = iter1 + iter2 + iter3;
                
                if (debug > 0){ printf("a %e %e %ld %ld %ld %e\n",eta,mu,iter1,iter2,iter3,result); }
                
                if (isnan(result)){
                   // print diagnostics
                   printf("pow(mu_outer,-p1): %f, ans0q: %f, ans1q: %f, ans2q: %f, ans3q: %f\n",
                           pow(mu_outer,-p1), ans0q, ans1q, ans2q, ans3q);
                }

            } else { // ; eta ~ mu
            
                if (mu_inner == 0.0){
                    
                    if (mu_inner == 0.0){
                       ans3q = 0.0;
                       iter3 = 0;
                    } else {
                       (void) gamma_infinite_right_natural(&p2, &mu_inner, &eta, &mu, &eps, &ans3q, &err, &iter3);
                       ans3q = ans3q * 8.0/pi;
                       //ans3q = gamma_infinite_right(p2, mu_inner, eta, mu, eps, nfunc=iter3)* 8.0/pi;
                    }
                
                    (void) gamma_finite_natural(&p1, &a, &mu0, &eta, &mu, &eps, &ans1q, &err, &iter1); 
                    ans1q = ans1q* 8.0/pi;
                    //ans1q = gamma_finite(p1, a, mu0, eta, mu, eps, nfunc=iter1)* 8.0/pi;
                    
                    (void) gamma_finite_natural(&p2, &a, &mu0, &eta, &mu, &eps, &ans2q, &err, &iter2); 
                    ans2q = ans2q* 8.0/pi;
                    //ans2q = gamma_finite(p2, a, mu0, eta, mu, eps, nfunc=iter2)* 8.0/pi;
                  
                    if (mu_outer > 0.0){
                       result = U1*pow(mu_outer,-p1)*ans0q + U1*ans1q + 
                                U2*( gam_analytic(p2,eta,mu) - (pow(mu_outer,-p2)*ans0q + ans2q) ) - U2*ans3q;
                    } else {
                        result = U1*ans1q + U2*( gam_analytic(p2,eta,mu) - ans2q ) - U2*ans3q;
                    }
                
                }
                
                if (mu_inner > 0.0){
                    
                   if (mu_outer > 0.0){
                      (void) gamma_finite_natural(&p1, &a, &mu0, &eta, &mu, &eps, &ans1q, &err, &iter1);
                      ans1q = ans1q * 8.0/pi;
                      //ans1q = gamma_finite(p1, a, mu0, eta, mu, eps, nfunc=iter1)* 8.0/pi;
                      
                      (void) gamma_finite_natural(&p2, &mu0, &mu_inner, &eta, &mu, &eps, &ans2q, &err, &iter2);
                      ans2q = ans2q* 8.0/pi;
                      //ans2q = gamma_finite(p2, mu0, mu_inner, eta, mu, eps, nfunc=iter2)* 8.0/pi;
                      iter3 = 0;
                      result = U1*pow(mu_outer,-p1)*ans0q + U1*ans1q + U2*ans2q;
                   } else {
                      (void) gamma_finite_natural(&p1, &a, &mu0, &eta, &mu, &eps, &ans1q, &err, &iter1);
                      ans1q = ans1q* 8.0/pi;
                      //ans1q = gamma_finite(p1, a, mu0, eta, mu, eps, nfunc=iter1)* 8.0/pi;
                      
                      (void) gamma_finite_natural(&p2, &mu0, &mu_inner, &eta, &mu, &eps, &ans2q, &err, &iter2);
                      ans2q = ans2q* 8.0/pi;
                      //ans2q = gamma_finite(p2, mu0, mu_inner, eta, mu, eps, nfunc=iter2)* 8.0/pi;
                      iter3 = 0;
                      result = U1*ans1q + U2*ans2q;
                   }
                }
                
                *err_ = 0.0;
                *nfunc_ = iter1 + iter2 + iter3;

                if (debug > 0){ printf("b %e %e %ld %ld %ld %e\n",eta,mu,iter1,iter2,iter3,result); }

            }
            //;print,'eta,mu,result,iter1,iter2: ',eta,mu,result[i],iter1,iter2
      
   }

   if (result < 0.0){ result = 0.0; }
   return result;
}

// To do:
// [ ] add description of gamma in terms of phi(chi) and give the latter (all options)
// [ ] check that mu_outer < mu0 < mu_inner
// [ ] if mu0 = 0 and p1=p2 then set mu0 = 1 (so that code U2 = U1 mu0^(p2-p1) gives U2=U1)

#ifdef MAIN
int main(int argc, char *argv[]){
    
    double result;
    double err;
    long nfunc;
    double Ustar;
    
    double eta;
    double mu;
    double U1;
    double U2;
    double p1;
    double p2;
    double mu0;
    double mu_outer;
    double mu_inner;
    
    int debug = 1;
 
    if (argc != 9){
    printf("Usage: gam_quadrature Ustar p1 p2 mu0 mu_outer mu_inner eta mu\n");
    exit(1);
    }

    Ustar  = atof(argv[1]);
    p1  = atof(argv[2]);
    p2  = atof(argv[3]);
    mu0 = atof(argv[4]);
    mu_outer = atof(argv[5]);
    mu_inner = atof(argv[6]);
    eta = atof(argv[7]);
    mu  = atof(argv[8]);

    if (mu0 >= 1){
       U1 = Ustar;
       U2 = U1 * pow(mu0,p2-p1);
    } else {
       U2 = Ustar;
       U1 = U2 * pow(mu0,p1-p2);
    }
    
    result = gam_new(eta, mu, 
               U1, U2, p1, p2, mu0, mu_outer, mu_inner, debug,
               &err, &nfunc);
    printf("gam_new: %22.15e, nfunc: %ld\n",result,nfunc);
    return 0;
}
#endif
