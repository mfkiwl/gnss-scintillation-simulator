
//#include "/usr/local/itt/idl71/external/include/idl_export.h" /* IDL external definitions */
//#include "/usr/local/exelis/idl/external/include/idl_export.h" /* IDL external definitions */

#ifndef IDL_LONG
   typedef long int IDL_LONG; 
#endif

double fgam(double x);
double fgam_prime(double x);
double fgam_asymptotic(double x);
double fgam_fluctuations(double x);
void get_eta_mu(double eta, double mu, double *etap, double *mup);

double fgam1(double x);
double fgam2(double x);
double fgam3(double x);

int gamma_integral_natural(double *p1, double *p2, double *chi0, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out);
int gamma_integral(int argc, void *argv[]);
int gamma_infinite_natural(double *p1, double *a, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out);
int gamma_infinite(int argc, void *argv[]);
int gamma_prime_infinite_natural(double *p1, double *a, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out);
int gamma_prime_infinite(int argc, void *argv[]);
int gamma_finite_natural(double *p1, double *a, double *b, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out);
int gamma_finite(int argc, void *argv[]);
int gamma_prime_finite_natural(double *p1, double *a, double *b, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out);
int gamma_prime_finite(int argc, void *argv[]);
int gamma_asymptotic_infinite_natural(double *p1, double *a, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out);
int gamma_asymptotic_infinite(int argc, void *argv[]);
int gamma_asymptotic_finite_natural(double *p1, double *a, double *b, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out);
int gamma_asymptotic_finite(int argc, void *argv[]);
int gamma_fluctuations_infinite_natural(double *p1, double *a, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out);
int gamma_fluctuations_infinite(int argc, void *argv[]);
int gamma_fluctuations_finite_natural(double *p1, double *a, double *b, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out);
int gamma_fluctuations_finite(int argc, void *argv[]);

int gamma_finite_left_natural(double *p1, double *a, double *b, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out);
int gamma_finite_left(int argc, void *argv[]);
int gamma_finite_right_natural(double *p1, double *a, double *b, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out);
int gamma_finite_right(int argc, void *argv[]);
int gamma_infinite_right_natural(double *p1, double *a, double *eta_in, double *mu_in, double *eps, 
    double *i, double *err, IDL_LONG *nfunc_out);
int gamma_infinite_right(int argc, void *argv[]);
int expm_natural(double *x, double *y);
double expm1(double);
int expm(int argc, void *argv[]);

    



