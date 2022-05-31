#include <stdio.h>
#include <math.h>

// expseries 
int main(void)
{
  printf("Dose x: ");
  double x;
  scanf("%lf", &x);

  double sum = 0.0;
  double a = 1.0; //  oros gia n=0. 
  
  for (int n = 1; sum +a != sum; ++n) {  
    sum += a;
    
    // o epomenos oros einai 
    a *= x / n;
  }

  printf("To athroisma sto  %f einai %.12f\n", x, sum);
  printf("H swsti timi einai %.12f\n", exp(x));

  return 0;
}
