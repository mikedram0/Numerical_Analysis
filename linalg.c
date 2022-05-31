#include <stdio.h>
#include <math.h>
#include <complex.h>
#define ITERMAX 100000
typedef double complex COMPL;


int COMPLtriangle_pivot( int N, COMPL A[], int Nb, COMPL B[]){
    int i,j,k,mi = 0;
    COMPL l, temp, maxc;
    double max,maxl;
    for(k = 0 ; k <= N-2 ; ++k){

        maxl = 0;
        for(i = k; i <= N-1 ; ++i){
            max = 0;
            for( j = k ; j <= N-1 ; j++ ){
                if( cabs(A[i+N*j]) > max){
                    max = cabs(A[i+N*j]);
		    maxc = A[i+N*j];		
		}
            }
            if(cabs(A[i+N*k]/maxc) > maxl ){
                maxl = cabs(A[i+N*k]/maxc);
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

    return l;

}

int triangle_pivot( int N, double A[], int Nb, double B[]){
    int i,j,k,mi = 0;
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

    return l;

}
int maketriangle( int N, double A[]){
    int i,j,k, mi = 0;
    int rowchanges = 1;
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
int COMPLbacktrack( int N, COMPL A[], int Nb, COMPL B[], double epsilon){

    int i,j,k = 0;
    COMPL sum,rhs = 0;
    for(k = 0; k <= Nb - 1; ++k){
    for(i = N-1 ; i >= 0 ; --i ){
        sum = 0;
        for(j = i+1 ; j <= N-1 ; ++j){
            sum += A[i+N*j]*B[j+N*k];
        } 
        rhs = (B[i+N*k] - sum);
        if(cabs(A[i+N*i]) < epsilon){
            if(cabs(rhs) < epsilon){
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

int gauss( int N, double A[],int Nb, double B[],double epsilon){
    triangle_pivot(N,A,Nb,B);
    return backtrack(N,A,Nb,B,epsilon);
}

int COMPLgauss( int N, COMPL A[],int Nb, COMPL B[], double epsilon){
    COMPLtriangle_pivot(N,A,Nb,B);
    return COMPLbacktrack(N,A,Nb,B,epsilon);
}
void printc( COMPL x){
    printf("%g%+gi\n", creal(x),cimag(x));
}

void COMPLprint( int M,int N, COMPL A[]){
    int i,j = 0;
    for(i = 0 ; i<= M-1 ; ++i ){
        for(j = 0; j <= N-1 ; ++j){
            printc(A[i+M*j]);
        }
        putchar('\n');
    }
    putchar('\n');
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
void eigenvec(int N, COMPL A[], COMPL roots[], COMPL result[],double epsilon){
    int i = 0;
    int j = 0;
    COMPL Anew[N*N];
    COMPL B[N];
    for(i = 0; i < N ; ++i){

        setmatc(N,1,B,0);
        copyto_complex(N,N, A, Anew);
        changediag(N,Anew,-1*roots[i]);
        COMPLgauss(N, Anew, 1, B, epsilon);
        COMPL sum = 0;
        for(j = 0 ; j < N ; ++j){
            sum += pow(cabs(B[j]),2);
        }
        sum = 1/(sqrt(sum));
        for(j = 0; j < N; ++j){
            result[j + i*N] = B[j]*sum;
	    printc(result[j+i*N]);
        }
        printf("DONE\n");
    
    }


}

int nonlinear(double *xi,  double (*F[])(double*, int), int N, double epsilon){
    double A[N*N];
    double b[N];
    int i,j;
    int k = 0;

    for( i = 0; i< N; ++i){
        b[i] = F[i](xi,0);
    }
    while( absmax(1,N,b) > epsilon){
        for(i = 0; i < N; ++i){
            for(j = 0; j < N; ++j){
                A[i+N*j] = F[i](xi,j+1);
            }
        }
        gauss(N,A,1,b,epsilon);
        matdiff(1,N,xi,b);
        for( i = 0; i< N; ++i){
            b[i] = F[i](xi,0);
        }
        k++;
    }
    printf("%d\n",k);
}
