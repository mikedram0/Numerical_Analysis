#include <stdio.h>
#include <math.h>

void deriv(double x, double y, double z, double dy[], double dz[])
{
  dy[0] = y + z*z - x*x*x;
  dz[0] = y*y*y + z + cos(x);

  dy[1] =-3.0 * x*x + dy[0] + 2.0 * z * dz[0];
			 
  dz[1] = 3.0 * y*y * dy[0] + dz[0] - sin(x);
 
  dy[2] = -6.0 * x + dy[1] + 2.0 * (dz[0]*dz[0] + z * dz[1]);
  dz[2] = 3.0 * y * (2.0 * dy[0]*dy[0] + y * dy[1]) + dz[1] - cos(x);

  dy[3] = -6.0 + dy[2] + 2.0 * (3.0 * dz[0] * dz[1] + z * dz[2]);
  
  dz[3] = dz[2] + sin(x) + 6.0 * dy[0] * dy[0] * dy[0]
    + 3.0 * y * (6.0 * dy[0] * dy[1] + y * dy[2]);

}




int main()
{
  double const a = 0.0;
  double const b = 1.0;
  double h = 0.1;

  double const yinit = 0.3;
  double const zinit = 0.1;

  double xnew, ynew, znew;
  double dy[4], dz[4];

  double xold = a;
  double yold = yinit;
  double zold = zinit;

  int nsteps = (b-a)/h;


  for (int k=0; k < nsteps; ++k) {
    if (xold+h > b) {
      h = b-xold;
    }
    xnew = xold + h;
        
    deriv(xold, yold, zold, dy, dz);
        
    ynew = yold + (dy[0] + ( dy[1] / 2 + (dy[2] / 6
					    + dy[3] / 24 * h) * h) * h) * h;
        
    znew = zold + (dz[0] + ( dz[1] / 2 + (dz[2] / 6
					    + dz[3] / 24 * h) * h) * h) * h;
        
    printf ("%g %g %g\n", xnew, ynew, znew);
    
    xold = xnew;
    yold = ynew;
    zold = znew;
  }

  return 0;
}
