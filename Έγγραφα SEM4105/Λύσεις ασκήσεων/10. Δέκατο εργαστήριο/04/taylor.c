#include <math.h>
#include <stdio.h>


void difeq(double x, double y, double dy[])
{    
  dy[0] = cos(x)-sin(y) + x*x;
  
  dy[1] = 2.0 * x - sin(x) - cos(y) * dy[0];
  
  dy[2] = 2.0 - cos(x) - cos(y) * dy[1] + sin(y) * dy[0]*dy[0];
  
  dy[3] = sin(x) + 3.0 * sin(y) * dy[0] * dy[1] +
    cos(y) * (pow(dy[0],3) - dy[2]);
  
  
  dy[4] = cos(x) + cos(y) * (6.0 * dy[0] * dy[0] * dy[1] - dy[3]) 
    + (3.0 * dy[1] * dy[1] + dy[0] * (4.0 * dy[2] - pow(dy[0],3)))
    * sin(y);
}

double taylorstep(double x0, double y0, double x1)
{
  double dy[5];
    
  difeq(x0, y0, dy);

  double y1 = y0;
  double h = x1-x0;
  double term = h;
    
  for (size_t i = 0; i < 5; ++i) {
    y1 += dy[i] * term;
    term *= h / (i+2);
  }

  return y1;
    
}


void taylor(double a, double b, double h, double ya, double * yb)
{
  double xold = a;
  double yold = ya;

  while (xold < b) {
    double xnew = xold + h;
    if (xnew > b) {
      xnew = b;
    }
	
    double ynew = taylorstep(xold, yold, xnew);

    xold = xnew;
    yold = ynew;
  }

  *yb = yold;
}


int main()
{
  double xa = -1.0;    
  double xb = 1.0;
  double h = 0.01;
  double ya = 3.0;
  double yb;
    
  taylor(xa, xb, h, ya, &yb);

  printf("%.5lf %.15lf\n", xb, yb);

  return 0;
}
