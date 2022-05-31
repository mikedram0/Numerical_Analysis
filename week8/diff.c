#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../mylib.h"
double fact( int x){
    double res = 1;
    int i;
    for (i = 1; i <= x; ++i){
        res *= i;
    }

    return res;
}
int main(){




    FILE *fpt = fopen("points.txt","r");
    if(fpt == NULL){
       printf("Could not open file, exiting");
       return -1; 
    }
    int N;
    fscanf(fpt, "%d", &N);
    double *x = malloc(N*sizeof(double));
    double *y = malloc(N*sizeof(double));
    int i;
    for( i = 0 ; i < N ; ++i){
        fscanf(fpt,"%lf", &x[i]); 
        fscanf(fpt,"%lf", &y[i]); 
        }

    fclose(fpt);

   // print(1,N,x);
   // print(1,N,y);
    double A[N*N];
    int j = 0;
    for( i = 0; i < N ; ++i){
        for(j = 0; j < N ; ++j){
            A[i + N*j] = pow(x[j],i); 
        }
    }

    int m = 1;
    double xi = 2;
    double b[N];
    for(i = 0; i<m ; ++i){
        b[i] = 0;
    }
    b[m] = fact(m);
    for( i = m + 1; i < N ; ++i){
        b[i] = b[i-1]*(i)*xi/(i - m); 
    }

    // print(N,N,A);
    // print(N,1,b);

    gauss(N,A,1,b,1e-14);
    double diff = 0;
    for(i = 0; i<N; ++i){
        diff += b[i]*y[i];
    }

    printf("%g\n", diff);

    return 0;
}
