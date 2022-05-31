#include <math.h>
#include <assert.h>
#include <complex.h>
#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>




void fillv(double complex v[], size_t n)
{
  for (size_t k=0; k < n/2; ++k) {
    v[k] = 2.0 / (1-4.0*k*k) -1.0 / (n*n-1.0 + n%2);
  }

  v[n/2] = (n-3.0)/( n - n%2 -1) - 1
    + 1.0 / (n*n-1.0 + n%2) * ((2 - n%2)*n-1);
    
  for (size_t k=1; k <= (n-1)/2; ++k) {
    v[n-k] = v[k];
  }
}



double const pi = 3.14159265358979323846;

bool isPowerOf2(size_t n)
{
  return n && !(n & (n - 1));
}



// with auxiliary vectors
void fft(double complex const f[], size_t N, double complex C[])
{
  if (N == 1) {
    C[0] = f[0];
    
    return;
  }
  
  assert(isPowerOf2(N));
 
  double complex * fhalf = malloc(sizeof(double complex) * N/2);
  double complex * ce = malloc(sizeof(double complex) * N/2);
  double complex * co = malloc(sizeof(double complex) * N/2);

  for (size_t i = 0; i < N/2; ++i) {
    fhalf[i] = f[2*i];
  }
  fft(fhalf, N/2, ce);

  for (size_t i = 0; i < N/2; ++i) {
    fhalf[i] = f[2*i+1];
  }
  fft(fhalf, N/2, co);

  free(fhalf);
    
  double complex pol = cexp(-2.0*I*pi/N);
  double complex coef=1.0;
  for (size_t i = 0; i < N/2; ++i) {
    co[i] *= coef;	
    C[i]     = 0.5*(ce[i] + co[i]);
    C[i+N/2] = C[i] - co[i];	
    coef *= pol;
  }
    
  free(co);
  free(ce);
}



// no auxiliary vectors
void fft_nv(double complex const *f, size_t N, double complex *C, size_t inc) 
{
  size_t const n = N / inc;
  if (n == 1) {
    *C = *f;
    return;
  }
  assert(isPowerOf2(n));

  double complex const * Feven = f;
  double complex const * Fodd  = f + inc;
  
  double complex * Ceven = C;
  double complex * Codd  = C + n/2;

  fft_nv(Feven, N, Ceven, 2*inc);
  fft_nv(Fodd,  N, Codd,  2*inc);
  
  double complex const pol = cexp(-2.0*I*pi/n);
  double complex coef=1.0;
  for (size_t i = 0; i < n/2; ++i) {
    *Codd *= coef;
    *Ceven = 0.5*(*Ceven + (*Codd));
    
    *Codd = *Ceven - (*Codd);

    ++Ceven;
    ++Codd;
    coef *= pol;
  }
}




void coefs(double complex w[], size_t n)
{
  --n;
  double complex * f = malloc(sizeof(double complex) * n);
  fillv(f,n);
    
  //fft(f, n, w);
  fft_nv(f, n, w, 1);

  w[n] = w[0];

  free(f);
}





/*
  \int_{-2}^{2} \frac{1}{1+x^2} \D x = \int_{-1}^{1} \frac{2}{1+4x^2} \D x
*/
double f(double x)
{
  return 2.0/ (1.0+4.0*x*x);
}



double calc(int n)
{
  double complex * w = malloc(sizeof(double complex) * (n+1));
  coefs(w,n+1);
    
  double complex s = 0.0;
  for (size_t i =0; i < n+1; ++i) {
    double x = cos(i*pi/n);
    s += w[i] * f(x);
  }

  free(w);
    
  return creal(s);  
}
    

int main()
{
  double const correct = 2.0*atan(2.0);
  size_t n = 2;
  while (true) {
    double s = calc(n);
    double diff = s-correct;
    printf("%.15g %.15g %.15g\n", s, correct, diff);
	
    if (fabs(diff) < 1e-12) {
      break;	    
    }
    n*=2;
  }

  printf("%d\n", n);
}

