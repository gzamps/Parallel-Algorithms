//=============================================================//
//              PageRank Algorithm with Pthreads               //
//                                                             //
//                       3 - 10 - 2014                         //                
//                    Zampokas Giorgos                         //   
//=============================================================//
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/types.h>
#include <math.h>
#include <pthread.h>
#include  <time.h>
#include <sys/time.h>
#define MAX_THREAD 1000
#define threshold 0.01;
#define  d = 0.85;
#define  od = 0.15;







void error_message(){

	char *help = "Error using PageRank: Three arguments required\n"
  "First: size of Graph\n"
  "Second: number of maximum connections\n"
  "Third: numder of threads\n";

  printf("%s", help);
}




typedef struct {
  int tid;
  int part;
  int elems_remained;
  int points_done;
  int cur_step;
} Worker;

  static int *graph;   // grafos me tis sundeseis
  static float *P;   // to dianysma arxikwn sinthikwn
  static float *E;   // dianysma E
  static int *done;  // dianysma me to an exei sugklinei i timi 
  static int *connections;   // dianysma me arithmo sundesewn 
  static float *PageRank;    // dianysma me timi PageRank
  static float *oldPageRank; // dianysma me tis proigoumenes times PageRank
  static unsigned int max_connections;      // max sundeseis
  static unsigned long int points;    // shmeia grafou
  static int *conv_step; //vhma sugklisis
  static int *points_done; //posa shmeia sugklinan

pthread_t *threads;
pthread_mutex_t mux1;
//pthread_barrier_t barrier; // barrier synchronization object


Worker *workers; // domi me ta dedomena


void *calcPageRank(void* worker){


  int i,j;
  float newP;
  float sum;
  int link;
  int elements = ((Worker *) worker)->part;
  int id =((Worker *) worker)->tid;
  int steps= ((Worker *) worker)->cur_step;
 

  for (i=0 ; i<elements ; i++){    //gia ola ta shmeia 
  newP=0.0;
  if (done[id*elements + i] == 0 ){       //an den exoun sugklinei
      sum = 0.0;
     
    for ( j=0 ; j< connections[id * elements +  i] ; j++){
      link = graph[id*elements*max_connections + j];
     
      if (connections[link] != 0){
      sum = sum + oldPageRank[link]/connections[link];         //to pagerank metavaletai !!
      
      }
    }

  newP = 0.85*sum; + 0.15*E[ id * elements + i];
  
  if (fabsf(oldPageRank[id*elements + i] - newP) < 0.0001 ){         //sugklisi

    done[id*elements + i]=1;
    
    
    points_done[ id ] ++;
    
    if (points_done[ id ] == elements)
    {
      conv_step[id]=steps;
      break;
          }
   // PageRank[i]=newP;
  }else{                                //oxi sugklisi
    PageRank[ id * elements + i]=newP;
  }
  
  }

  


}

}

int main(int argc, char **argv){

 


if(argc<4){
    error_message();
    return 0;
    //printf("Error using kmeans: Three arguments required\n");
}
  

// console input assignments
  points = atoi(argv[1]);
  max_connections = atoi(argv[2]);
  int i = 0;
  int j = 0;
  unsigned long int N = points;
  unsigned int c = max_connections;
  int pick = 0;
  unsigned long int nt=atoi(argv[3]);


  clock_t begin, end;
  double time_spent;


if((nt<1 || nt>MAX_THREAD)){
  printf("Enter valid number of threads, between 1 and %d !\n", MAX_THREAD);
  exit (1);
}
threads = (pthread_t *)malloc(nt*sizeof(pthread_t));

/*=======Memory Allocation=========*/  //  


  
  graph = (int*)malloc(points*max_connections*sizeof(int));
  P = (float*)malloc(points*sizeof(float));
  E = (float*)malloc(points*sizeof(float));
  PageRank= (float*)malloc(points*sizeof(float));
  oldPageRank= (float*)malloc(points*sizeof(float));
  done = (int*)malloc(points*sizeof(int));
  connections = (int*)malloc(points*sizeof(int));
  conv_step = (int*)malloc(nt*sizeof(int));
  points_done = (int*)malloc(nt*sizeof(int));



// OPEN FILE GRAPH
  FILE *fp;
  int* data = (int*)malloc(N*c*sizeof(int));

  fp = fopen("G1000000.bin", "rb");
  if(fp==NULL){printf("Error opening the file\n");}

  size_t w = fread(data, sizeof(int), N*c, fp);
  if(w!=N*c){printf("Error reading the data from G\n");}

  fclose(fp);


  for(i=0; i<points; i++){
    for(j=0; j<max_connections; j++){
      graph[i*max_connections + j] = data[i*max_connections + j]; 
    }
  }

//OPEN FILE P
  
  double* data2 = (double*)malloc(N*sizeof(double));

  fp = fopen("P1000000.bin", "rb");
  if(fp==NULL){printf("Error opening the file\n");}

   w = fread(data2, sizeof(double), N, fp);
  if(w!=N){printf("Error reading the data from P\n");}

  fclose(fp);


  for(i=0; i<points; i++){
    
      P[i] = data2[i] ;
    
  }

//OPEN FILE E
  
  double* data3 = (double*)malloc(N*sizeof(double));

  fp = fopen("E1000000.bin", "rb");
  if(fp==NULL){printf("Error opening the file\n");}

   w = fread(data3, sizeof(double), N, fp);
  if(w!=N){printf("Error reading the data from E\n");}

  fclose(fp);


  for(i=0; i<points; i++){
    
      E[i] = data3[i] ;
    
  }


//initialize points_done
for ( i = 0; i < nt; i++)
{
  points_done[i]=0;
}

// intialize done
  for (i = 0; i < points; i++)
  {
  	done[i]=0;
  }

//initialize PageRank
  for (i = 0; i < points; i++)
  {
  
  	oldPageRank[i]=P[i];
  	
  }


//divide dataset
  int part = points/nt; //posa shmeia antistoixoun se kathe processs
  int elems_remained = points - nt* part ; //posa shmeia perissevoun extra sto telos gia na ta parei h last process
  



begin = clock();


workers = malloc(nt*sizeof(Worker));
if (workers==NULL)
{
  printf("error me tin workers\n");
}
for (i=0;i<nt;i++){

  if (i!=nt){
    workers[i].tid=i ;
    workers[i].part=part;
   printf("%d\n", workers[i].part);

  
  }else{
     workers[i].tid=i ;
    workers[i].elems_remained=elems_remained;
  }
}
int steps;
int PointsConv=0;
for ( steps= 0; steps < 50; steps++){
  
  if(steps>0){
    for(j=0 ; j<points ; j++){

      oldPageRank[j]=PageRank[j];
    }
  }
  PointsConv=0;
for (i=0;i<nt;i++){
  PointsConv=PointsConv + points_done[i];
}
  if (PointsConv==points)
  {
    break;
  }
 printf("sunolika sugkl shmeia %d\n",PointsConv );
for (i=0;i<nt;i++){
  workers[i].cur_step = steps;
  pthread_create(&threads[i],NULL,calcPageRank,(void *)&workers[i]);
  
  pthread_join(threads[i],NULL);

}

}
for (i = 0; i < points; ++i){
  
}

int max_conv_step=0;
for (i=0 ;i<nt ; i++){
  if ( conv_step[i]>max_conv_step){
    max_conv_step=conv_step[i];
  }
}
  printf("sugklisi se %d vhmata\n",max_conv_step+1 );

  end = clock();     // end of execution timer
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;


  printf("Time elapsed: %f seconds. \n", time_spent);  


return 1;
}
