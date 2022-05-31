#include <stdio.h>
#include <float.h>


int main(void)
{
  float epsre = 1.0f;
  
  while (1.0f + epsre != 1.0f) {
    epsre *= 0.5f;
  }
  
  /*     here epsre is half the epsilon. */
  printf("%s %e\n", "float epsilon is ", 2.0f * epsre);
  printf("%s %e\n", "float epsilon from C is ", FLT_EPSILON);

  
  double epsdo = 1.0;
  
  while (1.0 + epsdo != 1.0) {
    epsdo *= 0.5;
  }
  
  /*     here epsdo is half the epsilon. */
  printf("%s %e\n", "double epsilon is ", 2.0 * epsdo);
  printf("%s %e\n", "double epsilon from C is ", DBL_EPSILON);

  return 0;
}
