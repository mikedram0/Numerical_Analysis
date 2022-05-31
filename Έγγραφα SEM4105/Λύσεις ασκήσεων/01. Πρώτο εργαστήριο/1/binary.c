#include <stdio.h>


int main(void)
{
  printf("Dose mh arnhtiko akeraio: ");

  int k;
  scanf("%d", &k);
  
  int b[32];
  
  int i = 0;

  do {
    b[i] = k%2;
    k /= 2;
    ++i;	
  } while (k != 0);

  do {
    --i;  
    printf("%d", b[i]);
  } while (i > 0);
  
  printf("\n");

  return 0;
}
