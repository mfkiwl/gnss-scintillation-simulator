/* test of intde1.c */

#include <math.h>
#include <stdio.h>

#include "intde1.c"

int nfunc;

main()
{
    double f1(double), f2(double), f3(double), f4(double), 
        f5(double), f6(double), f10(double), fgam(double);
    void intde(double (*f)(double), double a, double b, double eps, 
        double *i, double *err);
    void intdei(double (*f)(double), double a, double eps, 
        double *i, double *err);
    void intdeo(double (*f)(double), double a, double omega, double eps, 
        double *i, double *err);
    extern int nfunc;
    double i, err;

#if 1    
    //nfunc = 0;
    //intde(f1, 0.0, 1.0, 1.0e-15, &i, &err);
    //printf("I_1=int_0^1 1/sqrt(x) dx\n");
    //printf(" I_1= %lg\t, err= %lg\t, N= %d\n", i, err, nfunc);
    //nfunc = 0;
    //intde(f2, 0.0, 2.0, 1.0e-15, &i, &err);
    //printf("I_2=int_0^2 sqrt(4-x*x) dx\n");
    //printf(" I_2= %lg\t, err= %lg\t, N= %d\n", i, err, nfunc);
    
    nfunc = 0;
    intdei(f3, 0.0, 1.0e-15, &i, &err);
    printf("I_3=int_0^infty 1/(1+x*x) dx\n");
    printf(" I_3= %lg\t, err= %lg\t, N= %d\n", i, err, nfunc);
	printf(" I_3= %22.16e\t, err= %22.16e\t, N= %d\n", i, err, nfunc);

    nfunc = 0;
    intdei(f4, 0.0, 1.0e-15, &i, &err);
    printf("I_4=int_0^infty exp(-x)/sqrt(x) dx\n");
    printf(" I_4= %lg\t, err= %lg\t, N= %d\n", i, err, nfunc);
    printf(" I_4= %22.16e\t, err= %22.16e\t, N= %d\n", i, err, nfunc);

	//nfunc = 0;
    //intdeo(f5, 0.0, 1.0, 1.0e-15, &i, &err);
    //printf("I_5=int_0^infty sin(x)/x dx\n");
    //printf(" I_5= %lg\t, err= %lg\t, N= %d\n", i, err, nfunc);
	//printf(" I_5= %22.15e\t, err= %22.15e\t, N= %d\n", i, err, nfunc);


    //nfunc = 0;
    //intdeo(f6, 0.0, 1.0, 1.0e-15, &i, &err);
    //printf("I_6=int_0^infty cos(x)/sqrt(x) dx\n");
    //printf(" I_6= %lg\t, err= %lg\t, N= %d\n", i, err, nfunc);
	//printf(" I_6= %22.15e\t, err= %22.15e\t, N= %d\n", i, err, nfunc);
	
	//nfunc = 0;
    //intdeo(f10, 0.0, 100.0, 1.0e-15, &i, &err);
    //printf("I_10=int_0^infty exp(- ( (2*abs(x) + 2*abs(mu) -abs(x-mu) -abs(x+mu) )) * cos(x*mu) dx, mu=100\n");
    //printf(" I_10= %lg\t, err= %lg\t, N= %d\n", i, err, nfunc);
	//printf(" I_10= %22.15e\t, err= %22.15e\t, N= %d\n", i, err, nfunc);

#endif
#if 0 
    printf("\n");
    nfunc = 0;
    intdei(fgam, 0.0, 1.0e-6, &i, &err);
    printf("I_gam=int_0^infty x^(-p) * sin(x*eta/2)^2 * sin(x*mu/2)^2 dx\n");
    printf(" I_gam= %lg\t, err= %lg\t, N= %d\n", i, err, nfunc);
	printf(" I_gam= %22.15e\t, err= %22.15e\t, N= %d\n", i, err, nfunc);

#endif
}


double f1(double x)
{
    extern int nfunc;
    
    nfunc++;
    return 1 / sqrt(x);
}


double f2(double x)
{
    extern int nfunc;
    
    nfunc++;
    return sqrt(4 - x * x);
}


double f3(double x)
{
    extern int nfunc;
    
    nfunc++;
    return 1 / (1 + x * x);
}


double f4(double x)
{
    extern int nfunc;
    
    nfunc++;
    return exp(-x) / sqrt(x);
}


double f5(double x)
{
    extern int nfunc;
    
    nfunc++;
    return sin(x) / x;
}


double f6(double x)
{
    extern int nfunc;
    
    nfunc++;
	printf("nfunc: %d\n",nfunc);
    return cos(x) / sqrt(x);
}

double f10(double x)
{
    extern int nfunc;
    double mu = 100.;
    
    nfunc++;
    return 2.0*exp( -( (2.0*fabs(x) + 2.0*fabs(mu) - fabs(x-mu) - fabs(x+mu)) ) ) * cos(x*mu);
}

double fgam(double x)
{
    extern int nfunc;
    double eta = 1000.;
    double mu = 0.01;
    double p = 2.5;
    double trig;
    
    nfunc++;
   
    trig = sin(x*eta/2.0)*sin(x*mu/2.0);
    return pow(x,-p)*trig*trig;
}

