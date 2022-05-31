#include <stdio.h>
#include "../mylib.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>

typedef double complex COMPL;

int main(){
/*
    int N1 = 5;
    COMPL A[] = {0,1,0,0,1,1,0,5,0,0,5,1,0,2,0,0,5,1,0,1,1,0,0,1,0};
    //COMPLprint(N,N,A);

    COMPL roots[] = {5.021785902779,0,0,0,0};
    COMPL result[N1];

    eigenvec(N1,A,roots,result, 1e-5);
    COMPLprint(N1,1,result);
*/

    FILE *fptr = fopen("data.txt","r");
    if(fptr == NULL){
        printf("Could not open file");
        return -1;
    }
    int N2;
    fscanf(fptr,"%d", &N2);
    //printf("%d\n",N2);
    double *x = malloc(N2*sizeof(double));
    double *y = malloc(N2*sizeof(double));

    int i;
    for(i = 0; i < N2 ; ++i){
        fscanf(fptr,"%lf %lf", &x[i], &y[i]);
        //fscanf(fptr,"%g", &y[i]);
    }
    
    fclose(fptr);
 
    //print(N2,1,x);
    //print(N2,1,y);

    int M = 3;
    int N = 1;

    i = 0;
    int j = 0;
    double b[N+M+1];
    copyto(1,N+M+1,y,b);
    double A2[ (N+M+1)*(N+M+1) ];
    for (i = 0; i <= N+M ; ++i){
        for(j = 0; j <= M; ++j){
            A2[i + (N+M+1)*j] = pow(x[i],j);
    
        }
        for(j = M + 1; j <= M + N; j++ ){
            A2[i + (N+M+1)*j] = -y[i]*pow(x[i], j - M +1);
           
        }
    }
    print(N+M+1,N+M+1,A2);
    print(N+M+1,1,b);
    gauss(N+M+1,A2,1,b,1e-10);
    print(M+N+1,1,b);

    COMPL f(COMPL x){
        int i;
        COMPL sum1 = 0;
        COMPL sum2 = 0;
        for(i =0 ; i <= M; ++i ){
            sum1 += b[i]*cpow(x,i);
           // printc(sum1);
        }
        sum2 = 1 + b[4]*cpow(x,2);
        return sum1/sum2;
    }
    
//    COMPL f2(COMPL x){
//       return cpow(x,2) + 1;
//    }

   //printf("%g,%g,%g,%g,%g",b[0],b[1], b[2],b[3],b[4]);
   printc(temn(3,4,f,1e-10));
   printc(f(0.144639));
    return 0;
}
