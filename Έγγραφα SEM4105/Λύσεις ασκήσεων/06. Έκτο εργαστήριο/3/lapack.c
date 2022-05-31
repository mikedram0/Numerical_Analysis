#include <complex.h>
#include <stdio.h>
#include <stdlib.h>



void 
dgeev_(const char * jobvl, const char * jobvr, const int * n, 
       double a[], const int * lda, double wr[], double wi[], double vl[], 
       const int * ldvl, double vr[], const int * ldvr, double work[], 
       const int * lwork, int * info);


    
void eival(double a[], int n, double complex w[])
{
  double * wr = malloc(n*sizeof(double));
  double * wi = malloc(n*sizeof(double));
  
  int const ldvl = 1;
  double * vl = malloc(ldvl*n*sizeof(double));

  int const ldvr = 1;
  double * vr = malloc(ldvr*n*sizeof(double));

  double temp;
  int lwork = -1;
  int info = 0;
  char jobvl = 'n';
  char jobvr = 'n';
    
  dgeev_(&jobvl, &jobvr, &n, a, &n, wr, wi,
	 vl, &ldvl, vr, &ldvr, &temp, &lwork, &info);
    
  lwork = (int) temp;
  double * work = malloc(lwork*sizeof(double));

  info = 0;
    
  dgeev_(&jobvl, &jobvr, &n, a, &n, wr, wi,
	 vl, &ldvl, vr, &ldvr, work, &lwork, &info);

  if (info == 0) {
    for(int i=0; i < n; ++i) {
      w[i] = CMPLX(wr[i], wi[i]);
    }
  }

  free(work);
  free(vr);
  free(vl);
  free(wr);
  free(wi);
    
}


int main(void)
{
  size_t const n = 3;
  double * a = malloc(n*n*sizeof(double));

  a[0+3*0] = 6.3;
  a[0+3*1] = 2.1;
  a[0+3*2] = 4.15;

  a[1+3*0] = 3.1;
  a[1+3*1] = 5.14;
  a[1+3*2] = 1.03;

  a[2+3*0] = -11.0;
  a[2+3*1] = 12.3;
  a[2+3*2] = -8.8;


  double complex * w = malloc(n*sizeof(double complex));
    
  eival(a,n,w);

  for (size_t i = 0 ; i < n; ++i) {
    printf("%g %g\n", creal(w[i]), cimag(w[i]));
  }

  free(w);
  free(a);
    
  return 0;
}
