#include <stdio.h>
#include <string.h>

int len( char* str){

int i = 0;

	while(str[++i] != '\0'){
		if (isdigit(str[i]) == 0)  return -1;

}

	return i; 
}


int main(int argc, char* argv[]){

	if (argc != 2 ||  len(argv[1]) < 0){
		printf("Invalid parameters, Usage: dtb <num>\n");
		return 1;
}

	int num = atoi(argv[1]);
	int length =0;
	int tempnum = num;
	while(tempnum > 0){
		tempnum /= 2;
		length++;
}

	int digits[length];
	int i = 0;
	for(i = 0 ; i < length ; ++i ){
		int temp = num;
		digits[i] = temp % 2;
		num = num /2;
}
	for (i = length - 1 ; i >= 0 ; --i) 
		printf("%d", digits[i]);

	putchar('\n');
	return 0;
}

