#include <stdio.h>
#include <stdlib.h>
#include "../mylib.h"
#include <math.h>

void linspace(double a, double b, int N, double *x){

    int i = 0;
    double h = (b-a)/(N-1);
  //printf("%g\n", h);
    for( i = 0 ; i < N ; ++i ){
    
        x[i] = a + i*h;
    }
}

void applyfunc( double (*func)(double), double *x, double *y, int N  ){
    
  int i =0;
  for(i =0 ; i < N ; ++i){
    y[i] = func(x[i]);
  }
}

double func(double x){
    return sin(x);
}

double polynomial( double *x, double *y, int N, double xi, char method){
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
void copyto( int M, int N, double begin[], double end[]  ){
    int i = 0;
    for( i =0 ; i<= M*N - 1; ++i){
        end[i] = begin[i];
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
    
        gauss(N+1,A,1,b);
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
    gauss(N+M+1,A,1,b);
    //print(N+M+1,1,b);
    if(Ni = 0){
	copyto(b,result)	
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

int main(){

   // int N = 15;
   // double x[N]; 
   // double y[N];
   // linspace(2,4,N,x);
   // applyfunc(func, x,y,N);
   // print(1,N,x);
   // print(1,N,y);
    double x[] = {0.9,1.1,1.5,2.0,2.9,3.5};
    double y[] = {5.607,4.576,3.726,3.354,3.140,3.087};
    int i = 0;
    double x2[100];
    //linspace(2,4,100,x2);
    linspace(0.9,3.5,100,x2);
    double result[100];

    FILE* fptr = fopen("data.txt","w");
    if(fptr == NULL){
        printf("Could not open file");
        return -1;
    }
  //  polynomialgauss(x,y,N,x2,100,result);
    poldiv(x,y,2,3,100,x2,result);
    for( i = 0; i < 100; ++i){
        fprintf(fptr, "%g %g\n", x2[i],result[i] );
       // fprintf(fptr ,"%g %g\n",x2[i],polynomial(x,y,N,x2[i],'g'));
    }
    fclose(fptr);
    return 0;
}
