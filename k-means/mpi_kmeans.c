

//=============================================================================//
//                        k-means clustering with MPI                          //
//        benchmark file: kdd-cup-1999.bin (812900 points x34 dimensions)      //
//                            11 - 1 - 2014                                    //
//                           Zampokas Giorgos                                  //
//=============================================================================//



#include <stdio.h>
#include <stdlib.h>
#include <time.h>
//#include "kmeans.h"
//#include "cluster.h"
#include  <time.h>
#include <sys/time.h>
#include "mpi.h"
//#define max_iterations 50
#include <float.h>
#include <math.h>
#define threshold 0.001
//double threshold = 0.01;

typedef struct {
  double *dataset;
  unsigned int *members;
  int leading_dim;
  int secondary_dim; 
} data_struct;

//extern void cluster(data_struct *data_in, data_struct *clusters, int max_iterations);





void error_message(){

char *help = "Error using kmeans: Three arguments required\n"
  "First: number of elements\n"
  "Second: number of attributes (dimensions)\n"
  "Third: numder of clusters\n"
  "Fourth: number of max iterations\n";

  printf("%s", help);

}

void print(data_struct* data2print){

  int i, j = 0;
  int n = data2print->leading_dim;
  int m = data2print->secondary_dim;
  double *tmp_dataset = data2print->dataset;

  
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      printf("%f ", tmp_dataset[i*n + j]);
    }
    printf("\n");
  }
  
}


void save(data_struct* data2save, char *filename1, char *filename2){

  int i, j = 0;
  FILE *outfile;
  int n = data2save->leading_dim;
  int m = data2save->secondary_dim;
  double *tmp_dataset = data2save->dataset;
  unsigned int *tmp_members = data2save->members;

  printf("Saving data to files: "); printf( "%s",filename1); printf(" and "); printf("%s",filename2); printf("\n");

  /*===========Save to file 1===========*/
  if((outfile=fopen(filename1, "wb")) == NULL){
    printf("Can't open output file\n");
  }

  fwrite(tmp_dataset, sizeof(double), m*n, outfile);

  fclose(outfile);

  /*===========Save to file 2========*/

  if((outfile=fopen(filename2, "wb")) == NULL){
    printf("Can't open output file\n");
  }

  fwrite(tmp_members, sizeof(unsigned int), m, outfile);

  fclose(outfile);

}

void clean(data_struct* data1){

  free(data1->dataset);
  free(data1->members);
}



double euclidean_distance(double *v1, double *v2, int length){

  int i = 0;
  double dist = 0;

  for(i=0; i<length; i++){
    dist += (v1[i] - v2[i])*(v1[i] - v2[i]); 
  }

  return(dist);
}





int main(int argc, char **argv){

 //struct timeval first, second, lapsed;
  //struct timezone tzp;
clock_t begin, end;
double time_spent;



  if(argc<5){
    error_message();
    return 0;
    //printf("Error using kmeans: Three arguments required\n");
  }
  


  int numObjects = atoi(argv[1]);
  if (numObjects<=0 && numObjects>=1048576){
    printf("invalind number of objects\n");
    return 1;
  }
  int numAttributes = atoi(argv[2]);
  if (numAttributes<1 && numAttributes>128){
    printf("invalid number of attributes\n");
    return 2;
  }
  int numClusters = atoi(argv[3]);
  if (numClusters<2 && numClusters>32){
    printf("invalid number of clusters\n");
    return 3;
  }
  int max_iterations = atoi(argv[4]);
  if (max_iterations<1 && max_iterations>50){
    printf("invalid number of max_iterations\n");
  }
  
  int i = 0;
  int j = 0;

int SelfTID, NumTasks, t, data,wsize;
  MPI_Status mpistat;
  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &NumTasks );
  MPI_Comm_rank( MPI_COMM_WORLD, &SelfTID );




  char *file1_0 = "centroids.bin";
  char *file1_1 = "ClusterSize.bin";
  char *file2_0 = "dataset.bin";
  char *file2_1 = "Index.bin"; 



  data_struct data_in;
  data_struct clusters;

  

  //printf("Hello World from %i of %i\n", SelfTID , NumTasks);
  if (SelfTID==0 ){
    printf("There are %d points of %d dimensions.They will be seperated to  %d clusters using %d MPI processes\n  ", 
      numObjects, numAttributes ,numClusters , NumTasks);
  }

  
MPI_Barrier(MPI_COMM_WORLD);

  int part = numObjects/NumTasks; //number of points for each process
  int elems_remained = numObjects - NumTasks* part ; //number of points that remain after division.Those points will be 
  // taken by the last process
  

  //printf("elems_remained: %d\n", elems_remained);


if (SelfTID==0){
  /*=======Memory Allocation=========*/  // root 
  data_in.leading_dim = numAttributes;
  data_in.secondary_dim = numObjects;
  data_in.dataset = (double*)malloc(numObjects*numAttributes*sizeof(double));
  data_in.members = (unsigned int*)malloc(numObjects*sizeof(unsigned int));


  clusters.leading_dim = numAttributes;
  clusters.secondary_dim = numClusters;
  clusters.dataset = (double*)malloc(numClusters*numAttributes*sizeof(double));
  clusters.members = (unsigned int*)malloc(numClusters*sizeof(unsigned int)); 
 
}else if (SelfTID==  NumTasks -1){

  /*=======Memory Allocation=========*/   // last 
  data_in.leading_dim = numAttributes;
  data_in.secondary_dim = numObjects;
  data_in.dataset = (double*)malloc((numObjects*numAttributes/NumTasks + elems_remained*numAttributes )*sizeof(double));
  data_in.members = (unsigned int*)malloc((numObjects)*sizeof(unsigned int));

  clusters.leading_dim = numAttributes;
  clusters.secondary_dim = numClusters;
  clusters.dataset = (double*)malloc(numClusters*numAttributes*sizeof(double));
  clusters.members = (unsigned int*)malloc(numClusters*sizeof(unsigned int)); 

}else {

   /*=======Memory Allocation=========*/   // others
  data_in.leading_dim = numAttributes;
  data_in.secondary_dim = numObjects;
  data_in.dataset = (double*)malloc((numObjects*numAttributes/NumTasks  )*sizeof(double));
  data_in.members = (unsigned int*)malloc((numObjects)*sizeof(unsigned int));

  clusters.leading_dim = numAttributes;
  clusters.secondary_dim = numClusters;
  clusters.dataset = (double*)malloc(numClusters*numAttributes*sizeof(double));
  clusters.members = (unsigned int*)malloc(numClusters*sizeof(unsigned int)); 
}

   

if (SelfTID==NumTasks -1 ){
    printf("--> Diergasia %d: %d shmeia.  \n  ",SelfTID , numObjects/NumTasks + elems_remained );
  }else {
  printf("--> Diergasia %d: %d shmeia.  \n  ",SelfTID , numObjects/NumTasks );
}
 
   int iter;
  double SumOfDist = 0, new_SumOfDist = 0;
  double* newCentroids;
 
  
newCentroids = (double*)malloc(clusters.leading_dim*clusters.secondary_dim*sizeof(double));

 MPI_Barrier(MPI_COMM_WORLD);


  /*=============initialize ==========*/
  if (SelfTID==0){


  




  FILE *fp;



/*
Example of importing the kdd-cup-1999 data into a C program.
The data are a [819200 x 34] matrix stored column-wise.
author: Nikos Sismanis 
date: Jan 2014
*/


int N = numObjects;
int d = numAttributes;
int pick = 0;





  double* data = (double*)malloc(N*d*sizeof(double));

  fp = fopen("kdd-cup-1999.bin", "rb");
  if(fp==NULL){printf("Error opening the file\n");}

  size_t w = fread(data, sizeof(double), N*d, fp);
  if(w!=N*d){printf("Error reading the data\n");}

  


  for(i=0; i<data_in.secondary_dim; i++){
    data_in.members[i] = 0;
    for(j=0; j<data_in.leading_dim; j++){
      data_in.dataset[i*data_in.leading_dim + j] = data[i*data_in.leading_dim + j];
    }
    
  }



  int step = numObjects / numClusters;     
  
  for(i=0; i<clusters.secondary_dim; i++){
    for(j=0; j<clusters.leading_dim; j++){
      clusters.dataset[i*clusters.leading_dim + j] = data_in.dataset[pick * clusters.leading_dim + j];
    }
    pick += step; 
  }
fclose(fp);
  free(data);




  //random_initialization(&data_in);

  //initialize_clusters(&data_in, &clusters);

  begin = clock();  // Start of execution timer

  /*=================================*/
  MPI_Bcast( clusters.dataset, numClusters*numAttributes, MPI_DOUBLE, 0,
          MPI_COMM_WORLD);
  
} 


  for (j=0; j<numObjects ; j++){              //initialization of point memberships
   data_in.members[i]=0;
}



    if (SelfTID==0){

      for (i=1;i<=NumTasks-2;i++){
      MPI_Send(data_in.dataset+ i*numAttributes  ,numAttributes*part, MPI_DOUBLE, i,
         0, MPI_COMM_WORLD);
      }
         MPI_Send (data_in.dataset+ i*numAttributes   , numAttributes*part + numAttributes*elems_remained , MPI_DOUBLE, NumTasks-1,
          0, MPI_COMM_WORLD);
    }
   

    if (SelfTID != 0  && SelfTID < NumTasks -1 ){
        MPI_Recv( data_in.dataset, numAttributes*part, MPI_DOUBLE, 0,
         MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (SelfTID==NumTasks -1 ){
        MPI_Recv( data_in.dataset, numAttributes*part + numAttributes*elems_remained, MPI_DOUBLE, 0,
         MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }




int iters;
if(SelfTID==NumTasks-1){
  iters=part+elems_remained;
}else {
  iters=part;
}
//int totaliters=0;
//MPI_Reduce(  &iters , &totaliters, 1 ,
     //  MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);




  
   int k;
  double tmp_dist = 0;
  int tmp_index = 0;
  double min_dist = 0;
 
  
  
  
double part_distance;
   

unsigned int *dp=(unsigned int*)malloc(numObjects*sizeof(unsigned int));

unsigned int *cp=(unsigned int*)malloc(clusters.secondary_dim*sizeof(unsigned int)); 

double * newC = (double*)malloc(clusters.leading_dim*clusters.secondary_dim*sizeof(double));



  for(iter=0; iter<max_iterations; iter++){

    new_SumOfDist = 0;
    part_distance = 0;

    for(i=0; i<clusters.secondary_dim; i++){
      for(j=0; j<clusters.leading_dim; j++){
  newCentroids[i * clusters.leading_dim + j] = 0;
  newC[i * clusters.leading_dim + j] = 0;
      }
    }

    



for(i=0; i<clusters.secondary_dim; i++){
    cp[i] = 0;                                                 
    
}
for (i=0; i<numObjects;i++){
  dp[i]=0;
}


  for(i=0; i<iters; i++){
    tmp_dist = 0;
    tmp_index = 0;
    min_dist = FLT_MAX;
    /*find nearest center*/
    for(k=0; k<clusters.secondary_dim; k++){
      tmp_dist = euclidean_distance(data_in.dataset+  i*data_in.leading_dim, clusters.dataset +k*clusters.leading_dim  , data_in.leading_dim);
      if(tmp_dist<min_dist){
  min_dist = tmp_dist;
  tmp_index = k;     //which cluster is closer
      }
    }
   
    dp[SelfTID*part + i] = tmp_index;  //which cluster the point belongs to       
    part_distance += min_dist;  //total distance in this iteration
    cp[tmp_index]++;  //increment of cluster membership     
    for(j=0; j<data_in.leading_dim; j++){
      newC[tmp_index * clusters.leading_dim + j] += data_in.dataset[ i * data_in.leading_dim + j];  
     
    }
    
  }

 // Bcast + Reduce//

MPI_Reduce(  &part_distance , &new_SumOfDist, 1 ,
       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


if (SelfTID==0){
  printf("distance: %f\n", new_SumOfDist);
  for (i=1 ; i<=NumTasks-1 ; i++){
  MPI_Send(&new_SumOfDist ,1 , MPI_DOUBLE, i,
         0, MPI_COMM_WORLD);
  }

}else {
  MPI_Recv( &new_SumOfDist, 1, MPI_DOUBLE, 0,
         MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}


MPI_Reduce(dp,  data_in.members , data_in.secondary_dim,
           MPI_UNSIGNED, MPI_SUM, 0,
           MPI_COMM_WORLD);

MPI_Reduce(cp, clusters.members, clusters.secondary_dim,
           MPI_UNSIGNED, MPI_SUM, 0,
           MPI_COMM_WORLD);

MPI_Reduce(newC, newCentroids, clusters.leading_dim*clusters.secondary_dim,
           MPI_DOUBLE, MPI_SUM, 0,
           MPI_COMM_WORLD);




if (SelfTID==0){ //new cluster coordinates

  for(k=0; k<clusters.secondary_dim; k++){
    for(j=0; j<data_in.leading_dim; j++){
      clusters.dataset[k * clusters.leading_dim + j] = newCentroids[k * clusters.leading_dim + j] / (double) clusters.members[k];

    }
  }
}
 MPI_Bcast( clusters.dataset, clusters.leading_dim*clusters.secondary_dim, MPI_DOUBLE, 0,
          MPI_COMM_WORLD);



    if(fabs(SumOfDist - new_SumOfDist)<threshold){
      break;
    }
    SumOfDist = new_SumOfDist;

   
    

 }
  



//End of calculations!




free (dp);
free (cp);
free (newC);
free(newCentroids);



//Show results

if (SelfTID==0 ){

      
  end = clock();     // end of execution timer
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  printf("Finished after %d iterations.\n", iter);
  printf("The total distance between points and the cluster centers is %f. \n",new_SumOfDist );

  for (i=0 ; i<numClusters ;i++){
  printf("cluster [ %d ]->members: %u\n", i, clusters.members[i] );
  }
  


  printf("Time elapsed: %f seconds. \n", time_spent);  


  /*========save data============*/
  save(&clusters, file1_0, file1_1);
  save(&data_in, file2_0, file2_1);

  /*============clean memory===========*/
  clean(&data_in);
  clean(&clusters);
}



MPI_Finalize();
return 0;
}
