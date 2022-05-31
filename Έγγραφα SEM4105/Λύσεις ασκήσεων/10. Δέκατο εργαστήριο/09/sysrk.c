#include <stdio.h>
#include <stdlib.h>


double force(double x)
{
    return x * (1.0 - 0.01*x*x);
}



/*
  m x'' = F(x) =>    x' = v     =>    y0 = x, y1 = v
                     v' = F/m
*/
void f(double t, double const y[], double dy[])
{
    double const mass = 2.0;
    
    dy[0] = y[1];
    dy[1] = force(y[0]) / mass;
}



void ralston(size_t n, double t0, double const y0[], double t1, double y1[], 
	     void (*f)(double t, double const y[], double dy[]))
{
  double h = t1-t0;

  double *k1 = malloc(n * sizeof(double));
  double *k2 = malloc(n * sizeof(double));  

    f(t0, y0, k1);
    for (size_t i = 0; i < n; ++i) {
      k1[i] *= h;
    }


    double *p = malloc(n * sizeof(double));
    
    for (size_t i = 0; i < n; ++i) {
      p[i] = y0[i] + 2.0/3.0 * k1[i];
    }

    f(t0+2.0/3.0*h, p, k2);
    for (size_t i = 0; i < n; ++i) {
      k2[i] *= h;
    }

    for (size_t i = 0; i < n; ++i) {
	y1[i] = y0[i] + (k1[i] + 3.0 * k2[i])/4.0;
    }

    free(p);
    free(k1);
    free(k2);
}



int main()
{
    double t0 = 0.0;
    double x0 = 2.5e-2;
    double v0 = 0.0;
    
    size_t const nsteps = 100000;
    double const h = 0.001;

    double y0[] = {x0,v0};
    double y1[2];

    FILE * out =fopen("sysrk.txt", "w");

    for (size_t i = 0; i <= nsteps; ++i) {
      fprintf(out, "%.5lf %.15lf %.15lf\n", t0, y0[0], y0[1]);
	
      ralston(2, t0, y0, t0+h, y1, f);
	
	t0+=h;
	for (int j=0; j < 2; ++j) {
	  y0[j] = y1[j];
	}
    }

    fclose(out);
}
