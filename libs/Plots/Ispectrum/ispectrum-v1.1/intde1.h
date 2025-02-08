
//#include "intde1.c"

void intde(double (*f)(double), double a, double b, double eps, 
        double *i, double *err);
void intdei(double (*f)(double), double a, double eps, 
        double *i, double *err);
void intdeo(double (*f)(double), double a, double omega, double eps, 
        double *i, double *err);
