#include <assert.h>
#include <complex.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>



double const pi = 3.14159265358979323846;


bool isPowerOf2(size_t n)
{
    return n && !(n & (n - 1));
}



// with auxiliary vectors
void fft(size_t N, double complex const F[], double complex C[])
{
  assert(isPowerOf2(N));

  if (N == 1) {
    C[0] = F[0];
    
    return;
  }

  double complex * fhalf = malloc(N/2 * sizeof(double complex));
  double complex * ce    = malloc(N/2 * sizeof(double complex));
  double complex * co    = malloc(N/2 * sizeof(double complex));

  
  for (size_t i = 0; i < N/2; ++i) {
    fhalf[i] = F[2*i];
  }
  fft(N/2,fhalf, ce);
  
  for (size_t i=0; i < N/2; ++i) {
    fhalf[i] = F[2*i+1];
  }
  fft(N/2,fhalf, co);
  

  double complex const pol = cexp(-I*2.0*pi/N);
  double complex coef = 1.0;
  
  for (size_t i = 0; i < N/2; ++i) {

    co[i] *= coef;
    
    C[i]     = 0.5*(ce[i] + co[i]);
    
    C[i+N/2] = C[i] - co[i];
    
    coef *= pol;
  }
  
  free(co);
  free(ce);
  free(fhalf);
}



// with no auxiliary vectors
void fft_novec(size_t N, double complex const F[], double complex C[],
	       size_t inc)
{
  size_t n = N / inc;
  assert(isPowerOf2(n));

  if (n == 1) {
    C[0] = F[0];
    return;
  }

  double complex const * Feven = F;
  double complex const * Fodd = F+inc;

  double complex * Ceven = C;
  double complex * Codd = C+n/2;

  fft_novec(N, Feven, Ceven, 2*inc);
  fft_novec(N, Fodd,  Codd,  2*inc);

  double complex const pol = cexp(-I*2.0*pi/n);
  double complex coef = 1.0;

  for (size_t i = 0; i < n/2; ++i) {

    *Codd *= coef;

    *Ceven = 0.5*(*Ceven + (*Codd));

    *Codd = *Ceven -  (*Codd);
    
    ++Ceven;
    ++Codd;
    coef *= pol;
  }
}


int main()
{
#define N 1024

  double complex f[N], C[N];

  for (size_t i = 0 ; i < N; ++i) {
    f[i] = -0.5 + ((double) i)/N;
  }

   fft(N,f,C);
  //fft_novec(N,f,C,1);

    
  FILE * out = fopen("fourier.dat", "w");
  
  for (size_t i=0 ; i < N; ++i) {
    double complex correct = I /(2.0*i*pi);
    if (i==0) correct = 0.0;
    fprintf(out, "%lf %lf\t %lf %lf\n", creal(C[i]), cimag(C[i]),
	    creal(correct), cimag(correct));
  }

  fclose(out);
}
