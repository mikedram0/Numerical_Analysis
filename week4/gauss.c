#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void matmul(int N, double A[], double B[], double C[]){
    int i,j = 0;
    double sum = 0;
    for(i = 0; i <= N-1; ++i){
        sum = 0;
        for( j =0; j<= N-1 ; ++i){
            sum += A[i+N*j]*B[j];
        }
        C[i] = sum;
    }

}

void triangle_pivot( int N, double A[], int Nb, double B[]){
    int i,j,k, mi = 0;
    double l,max ,maxl, temp;
    for(k = 0 ; k <= N-2 ; ++k){

        maxl = 0;
        for(i = k; i <= N-1 ; ++i){
            max = 0;
            for( j = k ; j <= N-1 ; j++ ){
                if( fabs(A[i+N*j]) > max)
                    max = A[i+N*j];
            }
            if(fabs(A[i+N*k]/max) > maxl ){
                maxl = A[i+N*k]/max;
                mi = i ;
            }
        }

        temp = 0;
        for(j = 0; j<= N-1; ++j){
            temp = A[k+N*j];
            A[k+N*j] = A[mi+N*j];
            A[mi+N*j] = temp;
        }
        for(j = 0; j <= Nb - 1; j++){ 
            temp = B[k + N*j];
            B[k+N*j] = B[mi+ N*j];
            B[mi+N*j] = temp;
        }

        for(i = k + 1; i<= N-1 ; ++i){
            l = -1.0*A[i+N*k]/A[k+N*k];
            for(j = 0; j <= Nb-1 ; ++j){
                B[i+N*j] += l*B[k+N*j]; 
            }
            A[i+N*k] = 0;
            for(j = k + 1 ; j <= N-1 ; ++j ){
                A[i+N*j] += l*A[k+N*j];
            }    
        }   
    }



}
/*
void triangle( int N, double A[], double B[]){
    int i,j,k = 0;
    double l;
    for(k = 0 ; k <= N-2 ; ++k){
        for(i = k + 1; i<= N-1 ; ++i){
            l = -1.0*A[i+N*k]/A[k+N*k];
            B[i] += l*B[k];
            A[i+N*k] = 0;
            for(j = k + 1 ; j <= N-1 ; ++j ){
                A[i+N*j] += l*A[k+N*j];
            }    
        }   
    }



}
*/

int backtrack( int N, double A[], int Nb, double B[], double epsilon){

    int i,j,k = 0;
    double sum,rhs = 0;
    for(k = 0; k <= Nb - 1; ++k){
    for(i = N-1 ; i >= 0 ; --i ){
        sum = 0;
        for(j = i+1 ; j <= N-1 ; ++j){
            sum += A[i+N*j]*B[j+N*k];
        } 
        rhs = (B[i+N*k] - sum);
        if(fabs(A[i+N*i]) < epsilon){
            if(fabs(rhs) < epsilon){
                B[i+N*k] = 1;
                printf("Encountered Infinite solutions for solution number %d, choosing x[%d,%d] = 1\n", i , i, k);
            }else{
                printf("No solution");
                return -1;
            }
        }else{
            B[i+N*k] = (1/(A[i+N*i]))*rhs;
        }
    }
    }
    return 1;
}

int gauss( int N, double A[],int Nb, double B[]){
    triangle_pivot(N,A,Nb,B);
    return backtrack(N,A,Nb,B,1e-14);
}


void print( int M,int N, double A[]){
    int i,j = 0;
    for(i = 0 ; i<= M-1 ; ++i ){
        for(j = 0; j <= N-1 ; ++j){
            printf("%g ", A[i+M*j]);
        }
        putchar('\n');
    }
    putchar('\n');
}


int main(){
    FILE *fpt = fopen("gauss2.txt","r");
    if(fpt == NULL){
       printf("Could not open file, exiting");
       return -1; 
    }
    int N;
    int Nb = 2;
    fscanf(fpt, "%d", &N);
    double *A = malloc(N*N*sizeof(double));
    double *B = malloc(Nb*N*sizeof(double));
    int i,j = 0 ;
    for( i = 0 ; i<= N-1 ; ++i){
        for ( j = 0 ; j <= N-1 ; ++j){
            fscanf(fpt,"%lf", &A[i+N*j]); 
        }
    }
    for(i = 0; i<= N-1 ; ++i){
        fscanf(fpt, "%lf", &B[i]);
    }

    fclose(fpt);
    B[N] = 1;
    B[N+1] = 0;
    B[N+2] = 0;

    printf("File read succesfully\n");
    
    print(N,N,A);
    print(N,Nb,B);
    putchar('\n');
    
    if(gauss(N,A,Nb,B) > 0 ){
        print(N,N,A);
        print(N,Nb,B);
    }else{
        printf("No solution");
    }
    
    free(A);
    free(B);
    return 0;
}
