#include <stdio.h>
#include <math.h>
#include <stdlib.h>


// max abs value of a(row,:)
double maxval(int n, double const a[], int row)
{
  int mj = 0;
    
  for (int j = 1; j < n; ++j) {
    if (fabs(a[row+j*n]) > fabs(a[row+mj*n])) {
      mj=j;
    }
  }

  return fabs(a[row+mj*n]);
}

void swap(double * a, double *b)
{
  double t = *a;
  *a = *b;
  *b = t;
}

void pivot(int k, int n, double a[],int * change)
{
  int mi = k;
  double mvmi = maxval(n,a,mi);
  for (int i = k+1; i < n; ++i) {
    double  mvi = maxval(n,a,i);
    
    if (fabs(a[i+k*n] / mvi) > fabs(a[mi+k*n] / mvmi) ) {
      mi = i;
      mvmi = mvi;
    }
  }

  if (mi!=k) {
    for (int j=k; j < n; ++j) {
      swap(&a[mi+j*n], &a[k+j*n]);
    }
    ++*change;
  }
  
}



void triang(int n, double a[], int * change)
{
  for (int k = 0; k < n-1; ++k) {

    pivot(k,n,a,change);
	
    for (int i = k+1; i < n; ++i) {
      double const g = -a[i+n*k] / a[k+n*k];
      for (int j = k; j < n; ++j) {
	a[i+j*n] += g * a[k +j*n];
      }
    }
  }
}






double det(int n, double a[])
{
  int change = 0;
  
  triang(n,a,&change);

  double d = 1.0;
  for (int i=0; i < n; ++i) {
    d *= a[i+i*n];
  }

  return (change%2 == 0? d : -d);
}


double f(int n, double x, double a[])
{
  double * aa = malloc(n*n*sizeof(double));

  for(int i= 0; i < n*n; ++i) {
    aa[i] = a[i];
  }
  
  for(int i= 0; i < n; ++i) {
	aa[i+n*i] -= x;
    }

    return det(n,aa);
}


double
secant(double x1, double x2, double toler,
       double (*func)(int n, double x, double a[]),
       int n, double a[])
{
    double f1 = func(n,x1,a);
    double f2 = func(n,x2,a);
    double x;
    
    do {
	x = x2 - f2 * (x2 - x1) / (f2 - f1);

	x1 = x2;
	f1 = f2;

	x2 = x;
	f2 = func(n,x,a);

    } while (fabs(f2) > toler);

    return x;
}


int
main()
{
    int const n = 4;
    double a[n*n];

    a[0+n*0] = 2.1;
    a[0+n*1] = 3.9;
    a[0+n*2] = 0.3;
    a[0+n*3] = -4.1;
    a[1+n*0] = 4.3;
    a[1+n*1] = -1.3;
    a[1+n*2] = 0.8;
    a[1+n*3] = 1.5;
    a[2+n*0] = 1.0;
    a[2+n*1] = -2.8;
    a[2+n*2] = 4.3;
    a[2+n*3] = -8.1;
    a[3+n*0] = 2.4;
    a[3+n*1] = 6.1;
    a[3+n*2] = -1.1;
    a[3+n*3] = 12.5;

    double x1 = -1.0;
    double x2 = 1.0;

    double x = secant(x1, x2, 1e-8, f, n, a);

    printf("An eigenvalue is %g with determinant value %e\n", x, f(n,x,a));
}

