#include <stdio.h>
#include "../mylib.h"

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



int main(){
	int N = 3;
	double xi = 0.5;
	double yi = 0.5;
	double x[] = {0.1,0.3,0.5};
	double y[] = {0.2,0.4,0.6};
	double A[] = {0.2955,0.4794,0.6442,0.4794,0.6442,0.7833,0.6442,0.7833,0.8912};


	printf("%g\n",inter2d(x,N,y,N,A,xi,yi));

	

	return 0;
}