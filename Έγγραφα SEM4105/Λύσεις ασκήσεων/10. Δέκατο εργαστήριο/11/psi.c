#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
  y' = z
  z' = (x*x-5) * y   
*/

void f(double x, double const y[], double dy[])
{
  dy[0] = y[1];
  dy[1] = (x*x-5.0) * y[0];
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
  double const pi = 3.14159265358979323846;

#define nsteps  101

  double const a = -2.0, b = 2.0, h = (b-a)/(nsteps-1);
  int const offset = nsteps/2;

  double x0, y0[2], y1[2], x[nsteps], y[nsteps], v[nsteps];

  // nsteps/2 , 0.0, nsteps/2
  
  y0[0] = -1.0/sqrt(2.0*sqrt(pi));
  y0[1] =  0.0;

  y[offset] = y0[0];
  v[offset] = y0[1];

  for (int j = 0; j < nsteps; ++j) {
    x[j] = a + j * h;
  }

  
  for (int j = 0; j < nsteps/2; ++j) {
    x0 = x[offset+j];
    ralston(2, x0, y0, x0+h, y1, f);

    y0[0] = y1[0];
    y0[1] = y1[1];
    
    y[offset+(j+1)] = y1[0];
    v[offset+(j+1)] = y1[1];    
  }

  y0[0] = -1.0/sqrt(2.0*sqrt(pi));
  y0[1] =  0.0;

  for (int j = 0; j < nsteps/2; ++j) {
    x0 = x[offset-j]; 
    ralston(2, x0, y0, x0-h, y1, f);

    y0[0] = y1[0];
    y0[1] = y1[1];
    
    y[offset-(j+1)] = y1[0];
    v[offset-(j+1)] = y1[1];    
  }


  
  FILE * out = fopen("psi.txt", "w");
  for (int j = 0; j < nsteps; ++j) {
    fprintf(out, "%lf %.15lf %.15lf\n", x[j], y[j], v[j]);
  }
    
  fclose(out);

  return 0;
}
