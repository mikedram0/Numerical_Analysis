#include <stdio.h>


int main(void)
{
	/* Lathos tropos */
  double s = 0.0;
  for (double x=0.0; x <=1.9; x+=0.1) {
    s += x;
  }

  /* Swstos tropos */
  double p = 0.0;
  for (int i=0; i <=19; ++i) {
    p += 0.1 * i;
  }

  printf("%s %f %f\n", "To athroisma einai: ", s, p);

  return 0;
}
