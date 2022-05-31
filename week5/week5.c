#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "../mylib.h"
#define ITERMAX 10000
typedef double complex COMPL;
void printc( COMPL x){
    printf("%g%+gi\n", creal(x),cimag(x));
}
void copyto( int M, int N, double begin[], double end[]  ){
    int i = 0;
    for( i =0 ; i<= M*N - 1; ++i){
        end[i] = begin[i];
    }
}

void copyto_complex( int M, int N, COMPL begin[], COMPL end[]  ){
    int i = 0;
    for( i =0 ; i<= M*N - 1; ++i){
        end[i] = begin[i];
    }
}

double absmax( int M, int N, double A[]){
    int i = 0;
    double max = 0;
    double temp;
    for(i = 0; i<= M*N - 1; ++i){
        if( (temp = fabs(A[i]))  > max)
            max = temp;
    }
    return max;
}


int matdiff( int M, int N, double A[], double B[]){
    int i = 0;
    for(i == 0; i<= N*M - 1 ; ++i){
        A[i] -= B[i];
    }

}

double matsum_abs( int M, int N, double A[]){
    double sum = 0;
    int i =0;
    for(i = 0; i <= M*N -1 ; ++i){
        sum += fabs(A[i]);
    }
    return sum;
}


int matmul(int MA, int NA, double A[], int MB, int NB, double  B[], double result[]){
    if( NA != MB){
        printf("Invalid matrix size, inner dimensions must be equal");
        return(-1);
    }
    int i,j,k = 0;
    double sum = 0;

    for(i = 0 ; i <= MA -1 ; ++i){
        for( j =0 ; j <= NB - 1; ++j){
            sum = 0;
            for( k = 0; k <= NA - 1; ++k){
                sum += A[i + MA*k]*B[k+MB*j];
            }
            result[i + MA*j] = sum;
        }
    }
    return(1);
}

void GS(int N,double A[], double B[], double x[], double epsilon){

    int i,j,k = 0;
    double sum = 1;
    double result[N]; 
    matmul(N,N,A,N,1,x,result);
    matdiff(N,1,result,B);
    while(absmax(N,1,result) > epsilon){
        for(i = 0; i<= N-1 ; ++i){
            sum = 0;
            for(j = 0; j<= N-1 ; ++j){
                sum += A[i+N*j]*x[j];
            }
            x[i] = x[i] + (1/A[i+N*i])*(B[i] - sum) ;
        }
        matmul(N,N,A,N,1,x,result);
        matdiff(N,1,result,B);
        k++;
    }
    printf("Guass Seidel iterations for %e accuracy: %d\n",epsilon, k);
}


void jacobi(int N,double A[], double B[], double x[], double epsilon){

    double xnew[N];

    int i,j,k = 0;
    double sum = 1;
    double result[N]; 
    matmul(N,N,A,N,1,x,result);
    matdiff(N,1,result,B);
    while(absmax(N,1,result) > epsilon){
        for(i = 0; i<= N-1 ; ++i){
            sum = 0;
            for(j = 0; j<= N-1 ; ++j){
                sum += A[i+N*j]*x[j];
            }
            xnew[i] = x[i] + (1/A[i+N*i])*(B[i] - sum) ;
        }
        copyto(N,1,xnew,x);
        matmul(N,N,A,N,1,x,result);
        matdiff(N,1,result,B);
        k++;
    }
    printf("Jacobi iterations for %e accuracy: %d\n",epsilon, k);
}
void setmat(int M, int N, double A[],double val){
    int i = 0;
    for(i =0 ; i<= N*M -1 ; ++i){
        A[i] = val;
    }
}
void setmatc(int M, int N, COMPL A[],COMPL val){
    int i = 0;
    for(i =0 ; i<= N*M -1 ; ++i){
        A[i] = val;
    }
}
void changediag(int M, COMPL B[], COMPL num ){
    int i = 0;
    for(i = 0 ; i <= M-1 ; ++i){
        B[i+M*i] += num;
    }
    
}
int maketriangle_complex( int N, COMPL A[]){
    int i,j,k, mi = 0;
    int rowchanges = 1;
    COMPL l, temp;
    double max, maxl;
    for(k = 0 ; k <= N-2 ; ++k){

        maxl = 0;
        for(i = k; i <= N-1 ; ++i){
            max = 0;
            for( j = k ; j <= N-1 ; j++ ){
                if( cabs(A[i+N*j]) > max)
                    max = A[i+N*j];
            }
            if(cabs(A[i+N*k]/max) > maxl ){
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
	rowchanges *= -1;

        for(i = k + 1; i<= N-1 ; ++i){
            l = -1.0*A[i+N*k]/A[k+N*k];
            A[i+N*k] = 0;
            for(j = k + 1 ; j <= N-1 ; ++j ){
                A[i+N*j] += l*A[k+N*j];
            }    
        }   
    }

    return rowchanges;

}
COMPL det(int N, COMPL A[], COMPL lambda){
    COMPL B[N*N];
    copyto_complex(N,N,A,B);
    changediag(N, B, -1*lambda);
    COMPL det = 1;
    det = maketriangle_complex(N, B);
    int i = 0;
    for( i = 0; i<= N-1 ; ++i){
        det *= B[i+N*i]; 
    }
    return det;

}

COMPL mod_temn_matrix( COMPL a, COMPL b, COMPL (*func)(COMPL,COMPL*,int,int,COMPL*), double epsilon, COMPL* roots, int n, int N, COMPL A[]){
	int count = 0;
	COMPL f;
	COMPL funca = func(a,roots,n,N,A);
	COMPL funcb = func(b,roots,n,N,A);
	COMPL c = (b*funca - a*funcb)/(funca - funcb);
	while( cabs(f = func(c,roots,n,N,A)) > epsilon){
        if( count >= ITERMAX){ 
            printf("Could not converge to within ITERMAX, f(c) = %g%+gi\n", creal(f),cimag(f));  
            c = 0;
            break;
        }
		a = b;
		b = c;
		funca = func(a,roots,n,N,A);
		funcb = f;
		c = (b*funca - a*funcb)/(funca - funcb);
		count++;
	}
	printf("count: %d\n",count);
	return c;
}
COMPL char_pol (COMPL x, COMPL roots[], int n, int N, COMPL A[]){
	int i = 0;
    COMPL result = det(N, A, x);
	for(i = 0 ; i < n ; ++i)
		result /= (x-roots[i]);

	return result;
}
void eigenval(int N, COMPL A[], COMPL roots[] , COMPL a, COMPL b ){

    int i = 0;
    for(i = 0; i <= N-1  ; ++i){
 	    roots[i] = mod_temn_matrix(a,b,char_pol,1e-6,roots,i, N,A);        
	    printc(roots[i]);
    }

}
void eigenvec(int N, COMPL A[], COMPL roots[], COMPL result[]){
    int i = 0;
    int j = 0;
    COMPL Anew[N*N];
    COMPL B[N];
    for(i = 0; i < N ; ++i){
        setmatc(N,1,B,0);
        copyto_complex(N,N, A, Anew);
        changediag(N,Anew,-1*roots[i]);
        COMPLgauss(N, Anew, 1, B, 1e-6);
        COMPL sum = 0;
        for(j = 0 ; j < N ; ++j){
            sum += pow(cabs(B[j]),2);
        }
        sum = 1/(sqrt(sum));
        for(j = 0; j < N; ++j){
            result[j + i*N] = B[j]*sum;
        }
        
    
    }


}
int main(){
    
    int N1 = 4;
    COMPL A1[] = {2.1,4.3,1.0,2.4,3.9,-1.3,-2.8,6.1,0.3,0.8,4.3,-1.1,-4.1,1.5,-8.1,12.5};
    COMPL A2[] = {12.1,4.3,1.0,2.4,2.9,-11.3,-2.8,6.1,0.3,0.8,14.3,-1.1,-4.1,1.5,-8.1,12.5};
    double B2[] = {1.2,2.3,3.4,4.5};
    //double result1[] = {1,2,3,4};
    //double result2[] = {1,2,3,4};
    COMPL B[] = {1,2,3,4};
    
 

    //setdiag(N1,B,1);
    
    //gauss(N1,A,N1,B,1e-14);

    //matmul(N1,N1,A,N1,1,B,result);
    /*
    printf("%g\n", det(N1,A2));
    print(N1,N1,A2);
    print(N1,1,B2);

    jacobi(N1, A2, B2, result1, 1e-7);

    print(N1,1,result1);

    GS(N1,A2,B2,result2, 1e-7);
    print(N1,1,result2);
    */

    // COMPL A3[] = {0,-1,1,0};
    // COMPL result[2] ;
    COMPL result[N1*N1];
    eigenval(N1,A1,B,1+I,2+0.5*I);
    printf("\n");
    eigenvec(N1,A1,B,result);
    COMPLprint(N1,N1,result);
    

    return 0;
}
