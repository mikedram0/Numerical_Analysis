#include <stdio.h>
#include <math.h>
#include <stdlib.h>



/*     PRINT MATRIX */
void print(size_t n, double const a[])
{
  for (size_t i= 0 ; i < n; ++i) {
    for (size_t j= 0; j < n; ++j) {
      printf("%g ", a[i+n*j]);
    }
    printf("\n");
  }
  printf("\n");
}



void normalize(size_t n, double *a, double *b, size_t k)
{
    for (size_t j = k; j < n; ++j) {
      double maxv = fabs(b[j]);
	for (size_t p = k; p < n; ++p) {
	    double const v = fabs(a[j+n*p]);
	    if (v > maxv) {
		maxv = v;
	    }
	}
	
	for (size_t p = k; p < n; ++p) {
	    a[j+n*p] /= maxv;
	}
	for (size_t p = 0; p < n; ++p) {
	  b[j+n*p] /= maxv;
	}
    }
}


void swap(double *a, double *b)
{
    double tmp = *a;
    *a = *b;
    *b = tmp;
}


void swapa(size_t n, double *a, double *b, size_t k, size_t p)
{
    for (size_t j = k; j < n; ++j) {
      swap(&a[p + n*j], &a[k + n *j]);
    }
    for (size_t j = 0; j < n; ++j) {
      swap(&b[p+n*j], &b[k+n*j]);
    }
}


void pivot(size_t n, double *a, double *b, size_t k)
{

  normalize(n,a,b,k);

    size_t maxp = k;
    for (size_t p = k+1; p < n; ++p) {
	if (fabs(a[p+n*k]) > fabs(a[maxp+n*k])) {
	    maxp = p;
	}
    }

    if (maxp != k) {
      swapa(n,a,b,k,maxp);
    }
}

/*
  A : Ο πίνακας (δεδομένο)
  B : Ο αντίστροφός του (το "αποτέλεσμα" της υπορουτίνας).
*/  
void invert(size_t n, double a[], double b[])
{
  for (size_t i = 0; i < n*n; ++i) {
    b[i] = 0.0;	
  }
    
  for (size_t i = 0; i < n; ++i) {
    b[i+n*i] = 1.0;
  }


  for (size_t k = 0; k < n-1; ++k) {
    pivot(n, a,b,k);

    for (size_t i = k+1; i<n; ++i) {
      double g = -a[i+n*k] / a[k+n*k];

      for (size_t j = k; j < n; ++j) {
	a[i+j*n] += a[k+j*n] * g;
      }

      for (size_t j = 0; j < n; ++j) {
	b[i+j*n] += b[k+j*n] * g;
      }
    }
  }
    
  double * x = malloc(n*n*sizeof(double));
    
  for (size_t i = n; i-- > 0; ) {
    for (size_t j= 0; j < n; ++j) {

      x[i+j*n] = b[i+j*n];
      for (size_t k = i+1; k < n; ++k) {
	    
	x[i+j*n] -= a[i+k*n] * x[k+j*n];
      }
      x[i+j*n] /= a[i+i*n];
    }
  }
    
  for (size_t i = 0; i < n*n; ++i) {
    b[i] = x[i];
  }
    

  free(x);
}


int main()
{
  size_t const n = 4;
  double *a = malloc(n*n*sizeof(double));
  double *ainv = malloc(n*n*sizeof(double));

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

  invert(n,a, ainv);
  print(n,ainv);

  free(a);
  free(ainv);

  return 0;
}
