#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <math.h>

void error_message(){

char *help = "Error using PageRank: Three arguments required\n"
  "First: size of Graph\n"
  "Second: number of maximum connections\n"
  "Third: numder of clusters\n";

  printf("%s", help);

}




typedef struct {
  double *graph;
  float *P;
  float *E;
  int leading_dim; // max sundeseis
  int secondary_dim; // shmeia grafou
} data_struct;









int main(int argc, char **argv){

 



  if(argc<4){
    error_message();
    return 0;
    //printf("Error using kmeans: Three arguments required\n");
  }
  //printf("se kathe process!\n");


  int points = atoi(argv[1]);
  int max_connections = atoi(argv[2]);
  int i = 0;
  int j = 0;



printf("LOLOL\n");
data_struct data_in;


 /*=======Memory Allocation=========*/  //  
  data_in.leading_dim = max_connections;
  data_in.secondary_dim = points;
  data_in.graph = (double*)malloc(points*max_connections*sizeof(double));
  data_in.P = (float*)malloc(points*sizeof(float));
  data_in.E = (float*)malloc(points*sizeof(float));





int N = 1000;
int d = 15;
int pick = 0;

printf("here\n");
  // OPEN FILE GRAPH
  FILE *fp;
  double* data = (double*)malloc(N*d*sizeof(double));

  fp = fopen("G1000.bin", "rb");
  if(fp==NULL){printf("Error opening the file\n");}

  size_t w = fread(data, sizeof(double), N*d, fp);
  if(w!=N*d){printf("Error reading the data from G\n");}

  fclose(fp);


  for(i=0; i<data_in.secondary_dim; i++){
    for(j=0; j<data_in.leading_dim; j++){
      data_in.graph[i*max_connections + j] = data[i*max_connections + j]; 
    }
  }

//OPEN FILE P
  
  double* data2 = (double*)malloc(N*sizeof(double));

  fp = fopen("P1000.bin", "rb");
  if(fp==NULL){printf("Error opening the file\n");}

   w = fread(data2, sizeof(double), N, fp);
  if(w!=N){printf("Error reading the data from P\n");}

  fclose(fp);


  for(i=0; i<data_in.secondary_dim; i++){
    
      data_in.P[i] = data2[i] ;
    
  }

//OPEN FILE E
  
  double* data3 = (double*)malloc(N*sizeof(double));

  fp = fopen("E1000.bin", "rb");
  if(fp==NULL){printf("Error opening the file\n");}

   w = fread(data3, sizeof(double), N, fp);
  if(w!=N){printf("Error reading the data from E\n");}

  fclose(fp);


  for(i=0; i<data_in.secondary_dim; i++){
    
      data_in.E[i] = data3[i] ;
    
  }



return 1;
}
