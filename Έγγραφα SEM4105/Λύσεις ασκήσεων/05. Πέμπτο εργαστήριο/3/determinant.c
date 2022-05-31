#include <stdio.h>
#include <math.h>



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


int
main()
{
  int const n=4;
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

  printf("%g\n", det(n,a));
}

