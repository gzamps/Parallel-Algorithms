

 
#include <stdio.h>
#include <stdlib.h>
int main(){
	

int *array=malloc(1000000*sizeof(int));
unsigned long int i;
int a=20000;
float c=0;
for (i = 0; i < 1000000; i++)
{
    if (rand() % a == 0)
    {
        //array[i] = rand() % 100;
        array[i] = 1;
        c++;
    }
    else
    {
        array[i] = 0;
    }
}
for (i = 0; i < 1000000; i++){

	printf("%d" , array[i] , ",");
	//printf("\n");
}
printf("\n");
printf("sparse factor = %f", c/10 );
printf("\n");
return 0;
}