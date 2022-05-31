#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double const pi = 3.14159265358979323846;

void rk4(double x, double const y[], double h,
	 double (*f)(int i, double x, double const y[]), 
	 double * xnew, double ynew[], int n)
{

  double *k1 = malloc(n*sizeof(double));
  double *k2 = malloc(n*sizeof(double));
  double *k3 = malloc(n*sizeof(double));
  double *k4 = malloc(n*sizeof(double));

  double *p  = malloc(n*sizeof(double));

  
  for (int i=0; i < n; ++i) {
    k1[i] = h * f(i,x,y);
  }

  for (int i=0; i < n; ++i) {
    p[i] = y[i] + k1[i]/2.0;
  }

  for (int i=0; i < n; ++i) {	
    k2[i] = h * f(i,x+h/2.0,p);

  }

  for (int i=0; i < n; ++i) {
    p[i] = y[i] + k2[i]/2.0;
  }

  for (int i=0; i < n; ++i) {	
    k3[i] = h * f(i,x+h/2.0,p);
  }

  for (int i=0; i < n; ++i) {
    p[i] = y[i] + k3[i];
  }
    
  for (int i=0; i < n; ++i) {	
    k4[i] = h * f(i,x+h,p);
  }

  for (int i=0; i < n; ++i) {	
    ynew[i] = y[i] + (k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]) / 6.0;
  }
    
  *xnew = x + h;


  free(p);
  free(k4);
  free(k3);
  free(k2);
  free(k1);
}



/*
  The system is 
  D theta / D t  = z
  D z / D t      = -sin(theta)

  theta -> y[0], z-> y[1].  D theta / D t -> f[0], D z / D t -> f[1]. t -> x
*/

double f(int i, double x, double const y[])
{
  switch (i) {
  case 0:
    return y[1];
  case 1:
    return -sin(y[0]);
  default:
    fprintf(stderr, "error in f\n");
    exit(-1);
  }
}


/*
  solution of theta'' = - theta:   theta = A cos(t) + B sin(t)
  theta(0) = 45 deg. , theta'(0) = 0 => 
*/
void solution(double t, double y[])
{
  double const a = 45.0 * pi/180.0; 

  y[0] =  a * cos(t);
  y[1] = -a * sin(t);
}


int main()
{
#define N 2
  
  double const tinit = 0.0;
  double const tfin = 10.0;
  double const h = 0.1;

  double const yinit[N] = {45.0 / 180.0 * pi, 0.0};

  double told = tinit;
  double yold[N];
  double ynew[N];
  double yapprox[N];
  double tnew;
    
  for (int i=0; i < N; ++i) {
    yold[i] = yinit[i];
  }

  int const niter = (tfin-tinit)/h;
    
  printf("# t theta omega theta_approx omega_approx\n");
  
  for (int k=0; k < niter; ++k) {

    rk4(told,yold,h,f,&tnew,ynew, N);
    told = tnew;
    for (int i=0; i < N; ++i) {
      yold[i] = ynew[i];
    }

    solution(tnew, yapprox);
    printf("%g ", tnew);

    for (int i=0; i < N; ++i) {
      printf("%g ", ynew[i]);
    }
	
    for (int i=0; i < N; ++i) {
      printf("%g ", yapprox[i]);
    }
    printf("\n");
  }
}
