#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mylib.h"
double polynomial_interpolation( double *x, double *y, int N, double xi, char method){
    N = N - 1;
    double sum = 0;
    if(method == 'l'){
        double product = 1;
        int i,j = 0;
        for( i = 0; i<= N ; ++i){
            product = 1;
            for(j = 0; j <= N; ++j){
                if( i == j ){ continue; }
                product *= (xi - x[j])/(x[i]-x[j]);
            }
            sum += y[i]*product;
    
        }
    }

    if(method == 'n'){
        double sum2 = 0;
        double q = 1;
        double product = 1;
        int i = 0;
        int j = 0;
        int k = 0;
        double a[N+1];
        a[0] = y[0];
        sum += a[0];

        for(i = 1; i <= N; ++i){
            product *= (xi - x[i-1]);
            sum2 = 0;
            q = 1;
            sum2 += y[0]; 
            for( j = 1 ; j <= i-1 ; ++j){
                q = q*(x[i] - x[j-1] );
                sum2 += a[j]*q;
            }

            a[i] = 1/(q*(x[i] - x[i-1])) * ( y[i] - sum2  );
            sum += a[i]*product;
        }
    
    }
}
void polynomialgauss (double* x, double* y, int N, double* xi,int Ni, double* result){
 
        N = N-1;
        double A[(N+1)*(N+1)];
        double b[N+1];
        copyto(1,N+1,y,b);
        int i = 0;
        int j = 0;
        for( i = 0; i<= N ; ++i){
            for( j = 0; j <= N; ++j){
                A[i + (N+1)*j] = pow(x[i],j);
            }
        }
    
        gauss(N+1,A,1,b,1e-14);
        double sum = 0;
        for(j = 0; j < Ni ; ++j ){
            sum = 0;
            for( i = 0; i<= N; ++i){
                sum += b[i]*pow(xi[j],i);
            }
            result[j] = sum;
        }
    }

void poldiv(double *x, double *y, int M,int N,int Ni, double *xi, double *result){
    //N = N-1;
    //M = M-1;
    int i = 0;
    int j = 0;
    double b[N+M+1];
    copyto(1,N+M+1,y,b);
    double A[ (N+M+1)*(N+M+1) ];
    for (i = 0; i <= N+M ; ++i){
        for(j = 0; j <= M; ++j){
            A[i + (N+M+1)*j] = pow(x[i],j);
    
        }
        for(j = M + 1; j <= M + N; j++ ){
            A[i + (N+M+1)*j] = -y[i]*pow(x[i], j - M );
           
        }
    }
    //print(N+M+1,N+M+1,A);
    gauss(N+M+1,A,1,b,1e-14);
    //print(N+M+1,1,b);
    if(Ni = 0){
	copyto(1,N+M+1,b,result);	
	}else{
	
    double sum1 = 0;
    double sum2 = 1;
    for(i = 0; i < Ni ; ++i ){
        sum1 = 0;
        sum2 = 1;
        for(j = 0; j <= M ; ++j ){
            sum1 += b[j]*pow(xi[i],j);
        }
    
        for(j = 1; j <= N ; ++j ){
            sum2 += b[j+M]*pow(xi[i],j);
        }
        result[i] = sum1/sum2;
    }
    }
}
double inter2d( double *x, int Nx, double *y, int Ny, double* f, double xi, double yi ){

	int i,j;
	double tempf[Ny];
	double row[Nx];
	double result[1];
	double xxi[] = {xi};
	double yyi[] = {yi};
	for(i = 0; i < Ny; ++i){
		for(j = 0; j < Nx; ++j){
			row[j] = f[i + Ny*j];
		}
		polynomialgauss (x, row, Nx, xxi, 1 , result);
		tempf[i] = result[0];

	}
	polynomialgauss (y, tempf, Ny, yyi, 1 , result);
	return result[0];

}
double mean(double* x, int n){
    int i = 0;
    double sum = 0;
    for(i = 0; i < n ; ++i){
        sum += x[i];
    }
    sum /= n;
    return sum;
}

void lsq(double* x, double* y, int n, double *a, double *b){
    int i = 0;
    double xy = 0;
    double x2 = 0;
    double y2 = 0;
    double xm = mean(x,n);
    double ym = mean(y,n);


    for(i = 0; i < n; ++i){
        xy += x[i]*y[i];
        x2 += x[i]*x[i];
        y2 += y[i]*y[i];
    }
    xy /= n;
    x2 /= n;
    y2 /= n;

    *a = (xy - xm*ym)/(x2 - xm*xm);
    *b = (ym-(*a)*xm);

    double r2 = (*a)*(*a)*(x2 - xm*xm)/(y2 - ym*ym);
    printf("r^2 = %g\n", r2);
   
}

void fit( double *x, double *y, int n, double *a, double *b, char type){
    int i;
    double yy[n];
    double xx[n];
    switch(type){
        case 'l':
            lsq(x,y,n,a,b);
            break;

        case 'e':
            for(i = 0; i < n; ++i){
                xx[i] = exp(x[i]);
            }
            lsq(xx,y,n,a,b);
            break;

        case 'p':
            for(i = 0; i<n ; ++i){
                yy[i] = log(y[i]);
                xx[i] = log(x[i]);
            }
            lsq(xx,yy,n,a,b);
            *b = exp(*b);
            break;

        case 'g':
            for(i = 0; i<n; ++i){
                xx[i] = log(x[i]);
            }
            lsq(xx,y,n,a,b);
            break;
    
    }

}


void polynomial( double *x, double *y, double *res,int M, int N ){
    
    double *A = malloc((M+1)*(M+1)*sizeof(double));
    //double *B = malloc((M+1)*sizeof(double));

    int i,j,k;
    double sum1 = 0;
    double sum2 = 0;

    for(i = 0; i <= M ; ++i ){
        sum1 = 0;
        for(k = 0; k <= N-1 ; ++k){
            sum1 += pow(x[k],i)*y[k];
        }
        res[i] = sum1;
        for(j = 0; j <= M; ++j){
            sum2 = 0;
            for( k = 0; k <= N-1 ; ++k){
                sum2 += pow(x[k],i+j);
            }
            A[i + (M+1)*j] = sum2; 
        }
    }
    print((M+1),(M+1),A);
    print((M+1),1,res);
    gauss((M+1),A,1,res,1e-14);
    

}
void linspace(double a, double b, int N, double *x){

    int i = 0;
    double h = (b-a)/(N-1);
    for( i = 0 ; i < N ; ++i ){
    
        x[i] = a + i*h;
    }
}
void range(double a, double b, double h, double *x){

    int i = 0;
    int N = (int) (b-a)/h + 1;
    for( i = 0 ; i < N ; ++i ){
    
        x[i] = a + i*h;
    }
}
