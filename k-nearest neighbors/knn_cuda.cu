/* knn cuda */


#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "cuda.h"
#include <time.h>
#include <sys/time.h>



/*gia compiler 
export PATH=/Developer/NVIDIA/CUDA-5.5/bin:$PATH
 kai
export DYLD_LIBRARY_PATH=/Developer/NVIDIA/CUDA-5.5/lib:$DYLD_LIBRARY_PATH
*/



typedef struct{
  float *dataset;
  int leading_dim;
  int secondary_dim;
} knn_struct;


//#define BlockSize 512
//#define NUMBER_OF_SUBMATRICES 128

void error_message_fewer(){
	char *help = "Entered less than four arguments";

	char *help2 = "Error using knns: Four arguments required\n"
  "First: number of elements\n"
  "Second: number of attributes (dimensions)\n"
  "Third: numder of queries\n"
  "Fourth: number of nearest neighbours\n";

  printf("\e[1;34m %s \e[0m", help);
  	printf("\n");
  printf("\e[1;34m %s \e[0m", help2);
}

void error_message_more(){
	char *help = "Entered more than four arguments";

	char *help2 = "Error using knns: Four arguments required\n"
  "First: number of elements\n"
  "Second: number of attributes (dimensions)\n"
  "Third: numder of queries\n"
  "Fourth: number of nearest neighbours\n";

  printf("\e[1;34m %s \e[0m", help);
  	printf("\n");
  printf("\e[1;34m %s \e[0m", help2);
}

char* choose_data_file(int n){
	
	if (n==524288){
		return "base524288.bin";
	}else if(n==786432){
		return "base786432.bin";
	} else{
		return"base1048576.bin";
	}
	
	
}

char* choose_query_file(int q){
	
	
	if (q==1){
		return "query1.bin";
	}else if(q==100){
		return "query100.bin";
	}else if(q==200){
		return "query200.bin";
	}else if(q==300){
		return "query300.bin";
	}else if(q==400){
		return "query400.bin";
	}else if(q==500){
		return "query500.bin";
	}else if(q==600){
		return "query600.bin";
	}else if(q==700){
		return "query700.bin";
	}else if(q==800){
		return "query800.bin";
	}else if(q==900){
		return "query900.bin";
	}else {
		return "query1000.bin";
	}		
		
		
	
}

void save_distances(float* tmp_dataset, char *filename1,int n,int k){

  
  FILE *outfile;
  //int n = data2save->leading_dim;
  //int m = data2save->secondary_dim;
  //double *tmp_dataset = data2save->dataset;
  //unsigned int *tmp_members = data2save->members;

  printf("Saving data to files: "); printf(filename1);  printf("\n");

  /*===========Save to file 1===========*/
  if((outfile=fopen(filename1, "wb")) == NULL){
    printf("Can't open output file\n");
  }

  fwrite(tmp_dataset, sizeof(float), n*k, outfile);

  fclose(outfile);

}

void save_indexes(int* tmp_dataset, char *filename1,int n,int k){

  
  FILE *outfile;
  //int n = data2save->leading_dim;
  //int m = data2save->secondary_dim;
  //double *tmp_dataset = data2save->dataset;
  //unsigned int *tmp_members = data2save->members;

  printf("Saving data to files: "); printf(filename1); printf("\n");

  /*===========Save to file 1===========*/
  if((outfile=fopen(filename1, "wb")) == NULL){
    printf("Can't open output file\n");
  }

  fwrite(tmp_dataset, sizeof(int), n*k, outfile);

  fclose(outfile);

}


void cleanDevice(knn_struct *data){

  cudaFree(data->dataset);

}
/*
__global__ void selection_2(float* distances,float* NNdist,int* NNidx,int numObjects,int numQueries,int k){
	
	extern __shared__ float data[];
	
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int current_elem = tid;
	int current_block = blockIdx.x;
	
	__syncthreads();
	data[threadIdx.x]= distances[tid] 
	
}
*/

__global__ void select_kernel_last(float* next_distances_2, int*  next_indexes_2,float* NNdist,int* NNidx,int k, int q,int bs,float* last_distances,int* last_indexes){

	extern __shared__ float shared_data[];
	
	float* sdata = (float*)shared_data;
	int* sindexes = (int*)&sdata[512];
	
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	//int current_elem = tid;
	//int current_block = blockIdx.x;
	float temp_dist;
	int temp_idx;
	if (tid<bs){
	
	float kdist[8];
	int kidx[8];
	__syncthreads();
	//load elements in shared memory
	sdata[threadIdx.x] = next_distances_2[ tid ];
	sindexes[threadIdx.x] = next_indexes_2[ tid];
	
	__syncthreads();

	for (int neighbour = 0 ; neighbour < k ; neighbour ++){

	
		
		for (unsigned int s=1; s < blockDim.x; s *= 2) {
			
			int index = 2 * s * threadIdx.x;
			if (index < blockDim.x) {
				
				if (sdata[index]>sdata[index + s]){

					temp_dist=sdata[index];
					sdata[index]=sdata[index + s];
					sdata[index + s]=temp_dist;
			
					temp_idx=sindexes[index];
					sindexes[index]=sindexes[index + s];
					sindexes[index + s]=temp_idx;
				}
	
			}

__syncthreads();

		}

		if (threadIdx.x==0){
		kdist[neighbour]=sdata[0];
		kidx[neighbour]=sindexes[0];
		sdata[0]=FLT_MAX;
		}
	__syncthreads();	
	}


	if(threadIdx.x==0){
		for (int neighbour = 0 ; neighbour < k ; neighbour ++){
			NNdist[  q*k + neighbour ] = kdist[neighbour];
			last_distances[ neighbour ]=kdist[neighbour];
		//distances[q*numObjects + neighbour]=sdata[neighbour];
			NNidx[  q*k + neighbour ] = kidx[neighbour];
			last_indexes[neighbour ]=kidx[neighbour];
		}
	}
	}

}
__global__ void select_kernel_2(float* next_distances,int* next_indexes,float* next_distances_2,int* next_indexes_2,int numQueries,int k ){
	
	extern __shared__ float shared_data[];
	
	float* sdata = (float*)shared_data;
	int* sindexes = (int*)&sdata[512];
	
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	//int current_elem = tid;
	//int current_block = blockIdx.x;
	float temp_dist;
	int temp_idx;
	
	
	float kdist[8];
	int kidx[8];

	__syncthreads();
	//load elements in shared memory
	sdata[threadIdx.x] = next_distances[ tid ];
	sindexes[threadIdx.x] = next_indexes[ tid ];
	
	__syncthreads();
	
	for (int neighbour = 0 ; neighbour < k ; neighbour ++){
		
		for (unsigned int s=1; s < blockDim.x; s *= 2) {
			
			int index = 2 * s * threadIdx.x;
			if (index < blockDim.x) {
				
				if (sdata[index]>sdata[index + s]){

					temp_dist=sdata[index];
					sdata[index]=sdata[index + s];
					sdata[index + s]=temp_dist;
			
					temp_idx=sindexes[index];
					sindexes[index]=sindexes[index + s];
					sindexes[index + s]=temp_idx;
				}
	
			}

		__syncthreads();

		}
		if (threadIdx.x==0){
		kdist[neighbour]=sdata[0];
		kidx[neighbour]=sindexes[0];
		sdata[0]=FLT_MAX;
		}
		__syncthreads();
	}


	if(threadIdx.x==0){
		for (int neighbour = 0 ; neighbour < k ; neighbour ++){
			next_distances_2[  k*blockIdx.x + neighbour ] = kdist[neighbour];
		//distances[q*numObjects + neighbour]=sdata[neighbour];
			next_indexes_2[ k*blockIdx.x + neighbour ] = kidx[neighbour];
		}
	}
	
	
}




__global__ void select_kernel(float* distances,float* next_distances,int* next_indexes,int numObjects,int k){
	
	extern __shared__ float shared_data[];
	
	float* sdata = (float*)shared_data; 
	 int* sindexes = (int*)&sdata[512]; 
	  
	
	
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	//int current_elem = tid;
	//int current_block = blockIdx.x;
	float temp_dist;
	int temp_idx;
	
	
	float kdist[8];
	int kidx[8];

	
	
	__syncthreads();
	
	//load elements in shared memory
	sdata[threadIdx.x] = distances[ tid ];
	sindexes[threadIdx.x] = tid;
	//sdata[0]=FLT_MAX;
	
	
	__syncthreads();
	
	
	for (int neighbour = 0 ; neighbour < k ; neighbour ++){
		



		for (unsigned int s=1; s < blockDim.x; s *= 2) {
			
			int index = 2 * s * threadIdx.x;
			if (index < blockDim.x) {
				
				if (sdata[index]>sdata[index + s]){

					temp_dist=sdata[index];
					sdata[index]=sdata[index + s];
					sdata[index + s]=temp_dist;
			
					temp_idx=sindexes[index];
					sindexes[index]=sindexes[index + s];
					sindexes[index + s]=temp_idx;
				}
	
			}



		}
__syncthreads();
		if (threadIdx.x==0){
		kdist[neighbour]=sdata[0];
		kidx[neighbour]=sindexes[0];
		sdata[0]=FLT_MAX;
		


		}
		
__syncthreads();




	}


	if(threadIdx.x==0){
		for (int neighbour = 0 ; neighbour < k ; neighbour ++){
			next_distances[ k*blockIdx.x + neighbour ] = kdist[neighbour];
		//distances[q*numObjects + neighbour]=sdata[neighbour];
			next_indexes[ k*blockIdx.x + neighbour ] = kidx[neighbour];
		}
	}
	
	
}


__device__ float euclidean_distance_gpu(float *v1, float *v2, int attributes, int numObjects){

  float dist = 0;
  
#pragma unroll 2
  for( int i = 0; i < attributes; i++ ){
    float tmp = v2[i*numObjects] - v1[i];
    dist += tmp * tmp;
  }
  return dist;
}

__device__ float my_euclidean_distance_gpu(float *v1, float* v2, int attributes,int numObjects){

	float dist = 0;
//#pragma unroll 2
	for(int i=0; i<attributes; i++){
		float tmp= v2[i]-v1[i];
		dist += tmp * tmp;
	}
	return dist;
}


__global__ void calculate_distances_seperately(float* dataset ,float* queries ,float* distances,int numObjects,int numAttributes,int numQueries,int k,int q){
	
	extern __shared__ float querymeans[];
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int current_elem = tid;
	float tmp_dist;
	
	if (tid<numObjects){
		
			tmp_dist=0;
			
			__syncthreads();
			if ( threadIdx.x < numAttributes){
				querymeans[threadIdx.x]= queries[q*numAttributes + threadIdx.x];
				}
			__syncthreads();
			//tmp_dist = euclidean_distance_gpu( querymeans, dataset + current_elem, numAttributes, numObjects);
			tmp_dist = my_euclidean_distance_gpu( querymeans, dataset + current_elem*numAttributes , numAttributes, numObjects);
			distances[tid]=tmp_dist;
		
		
	}
	

	
	
}
__global__ void clear_distances(float* distances, int numObjects){

	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid<numObjects){
		distances[tid]=0.0;
	}
}

void knns(knn_struct* d_training_set, knn_struct* d_query_set, float* d_NNdist, int* d_NNidx, int tk,float* distances){
	
	float* dataset=d_training_set->dataset;
	float* queries=d_query_set->dataset;
	int k=tk;
	float* NNdist=d_NNdist;
	int* NNidx=d_NNidx;
	int numObjects=d_training_set->secondary_dim;
	int numAttributes=d_training_set->leading_dim;
	int numQueries=d_query_set->secondary_dim;
	//size_t memory_free, memory_total;
	//int i;
	int BlockSize=512;

	
#ifdef TIMEONLY
  float elapsedTime_kernel;
  cudaEvent_t start_kernel, stop_kernel;
  cudaEventCreate(&start_kernel);
  cudaEventCreate(&stop_kernel);
#endif
	
  
printf("tha diairesw %d dia  %d = %d",numObjects , BlockSize , numObjects/BlockSize);
	int tmp_grid_size = numObjects/BlockSize;
  int tmp_block_size = numObjects<BlockSize ? numObjects:BlockSize; 

  

  printf("tmp_grid_size = %d \n",tmp_grid_size );
  dim3 grid(tmp_grid_size,1);
  dim3 threads(tmp_block_size, 1);
  printf(" \n");
 
  
  
#ifdef TIMEONLY
  cudaEventRecord(start_kernel, 0);
#endif

  
  
int blocks_reducing=numObjects/BlockSize; 
int blocks_reducing_2=blocks_reducing/BlockSize; 
  float* next_distances;
  float* next_distances_2;
  int* next_indexes;
  int* next_indexes_2;
/*
  float* test_nd=(float *)malloc(k*blocks_reducing*sizeof(float));
  int* test_ni=(int *)malloc(k*blocks_reducing*sizeof(int));
  float* test_nd_2=(float *)malloc(k*k*blocks_reducing_2*sizeof(float));
  int* test_ni_2=(int *)malloc(k*k*blocks_reducing_2*sizeof(int));
  */
  
  cudaMalloc((void**)&next_indexes, k*blocks_reducing*sizeof(int));
  cudaMalloc((void**)&next_distances, k*blocks_reducing*sizeof(float));
/*
int* next_indexes_s;
float* next_distances_s;
  cudaMalloc((void**)&next_indexes_s, k*blocks_reducing*sizeof(int));
  cudaMalloc((void**)&next_distances_s, k*blocks_reducing*sizeof(float));
*/
  cudaMalloc((void**)&next_indexes_2, k*k*blocks_reducing_2*sizeof(int));
  cudaMalloc((void**)&next_distances_2, k*k*blocks_reducing_2*sizeof(float));


  float* last_distances;
  int* last_indexes;
  cudaMalloc((void**)&last_distances, k*sizeof(float));
  cudaMalloc((void**)&last_indexes, k*sizeof(int));
  //float* testld=(float *)malloc(k*sizeof(float));
 // int* testid=(int *)malloc(k*sizeof(int));
  printf("knn starts \n");




//test
//  float *te=(float *)malloc(1*sizeof(float));
//float *te1=(float *)malloc(1*sizeof(float));
//float *te2=(float *)malloc(1*sizeof(float));
//float *t=(float *)malloc((numObjects-1000000)*sizeof(float));

//====== for each query =====//
for (int q=0; q<numQueries; q++){
 
 
  clear_distances<<<grid ,threads >>>(distances, numObjects);
  
  calculate_distances_seperately<<<grid , threads , numAttributes*sizeof(float)>>>(dataset, queries, distances, numObjects, numAttributes, numQueries, k, q);
 
/*
 cudaMemcpy(t, distances + 1000000, (numObjects-1000000)*sizeof(float), cudaMemcpyDeviceToHost);
  for(i=0;i< numObjects-1000000 -48000 ;i++){
	printf("stoixio[%d] = %f ,,,,, ", i+1000000 , t[ i ]);
}
  /*
  //cudaMemcpy(testdist, distances, numObjects*sizeof(float), cudaMemcpyDeviceToHost);
  for (i=500000;i<500500;i++){
  	printf("| %d-> %f |.",i, testdist[i]);
  }
  */
  //	cudaMemcpy(te, distances + 1000021, 1*sizeof(float), cudaMemcpyDeviceToHost);
  	//printf("meta to cal dist, te:%f\n",te[0]);
  select_kernel<<<  grid, threads, 2*BlockSize*sizeof(float) >>>(distances, next_distances, next_indexes, numObjects, k);
	


//cudaMemcpy(te1, next_distances + k*blocks_reducing- 5, 1*sizeof(float), cudaMemcpyDeviceToHost);
//printf("meta to cal dist, te1:%f\n",te1[0]);

/*
cudaMemcpy(test_ni, next_indexes, k*blocks_reducing*sizeof(int), cudaMemcpyDeviceToHost);
cudaMemcpy(test_nd, next_distances, k*blocks_reducing*sizeof(float), cudaMemcpyDeviceToHost);


	 for(i=0;i<k*blocks_reducing;i++){
	 	printf("next_indexes[%d]= %d me dist =%f\n", i,test_ni[i],test_nd[i]);
	 }
*/
	 dim3 newgrid((int)k*blocks_reducing/BlockSize,1);
	 dim3 newthreads((int)BlockSize,1);
	 
	 select_kernel_2<<<  newgrid, newthreads, 2*BlockSize*sizeof(float) >>>( next_distances, next_indexes, next_distances_2, next_indexes_2, numQueries, k);
	
	//cudaMemcpy(te2, distances, 1*sizeof(float), cudaMemcpyDeviceToHost);

	 //cudaMemcpy(test_ni_2, next_indexes_2, k*k*blocks_reducing_2*sizeof(int), cudaMemcpyDeviceToHost);
     //cudaMemcpy(test_nd_2, next_distances_2, k*k*blocks_reducing_2*sizeof(float), cudaMemcpyDeviceToHost);
/*
	for(i=0;i<k*k*blocks_reducing_2;i++){
	 	printf("next_indexes_2[%d]= %d me dist =%f \n", i,test_ni_2[i],test_nd_2[i]);
	 }
*/

	 dim3 lastgrid(1,1);
	 dim3 lastthreads((int)BlockSize,1); //2 -> 64 

	 select_kernel_last<<<  lastgrid, lastthreads, 2*BlockSize*sizeof(float) >>>( next_distances_2, next_indexes_2, NNdist, NNidx, k ,q, k*k*blocks_reducing_2, last_distances,last_indexes);
	 //cudaMemcpy(testld, last_distances, k*sizeof(float), cudaMemcpyDeviceToHost);
	// cudaMemcpy(testid, last_indexes, k*sizeof(int), cudaMemcpyDeviceToHost);
	 /*
	 for (int i = 0; i < k; i++)
	 {
	 	printf("last_distances= %f apo index= %d \n",testld[i],testid[i] );
	 }
	 */
  }
  
 


  
  printf("Done with knns\n");
 

  
#ifdef TIMEONLY
  cudaEventRecord(stop_kernel, 0);  
  cudaEventSynchronize(stop_kernel);
#endif
  
  
  
#ifdef TIMEONLY
  cudaEventElapsedTime(&elapsedTime_kernel, start_kernel, stop_kernel);
  printf("Time elapsed for kernel execution: %f ms\n", elapsedTime_kernel);
#endif
	
	
  cudaFree(distances);
  cudaFree(next_distances);
  cudaFree(next_indexes);

#ifdef TIMEONLY
  cudaEventDestroy(start_kernel);
  cudaEventDestroy(stop_kernel);
#endif

}	
		





int main(int argc, char **argv){


  struct timeval first, second, lapsed;
  struct timezone tzp;
  size_t memory_free, memory_total;
  
  
cuMemGetInfo_v2(&memory_free, &memory_total);
  printf("Totel memory: %zd, free memory: %zd\n", memory_total, memory_free);

  
 

  if(argc<5){
    error_message_fewer();
    return 0;
  }
  if (argc>5){
  	error_message_more();
  }


  int numObjects = atoi(argv[1]);
  if (numObjects<524288 || numObjects>1048576){
    printf("invalind number of objects\n");
    return 1;
  }
  int numAttributes = atoi(argv[2]);
  if (numAttributes<128 || numAttributes>128){
    printf("invalid number of attributes\n");
    return 2;
  }
  int numQueries = atoi(argv[3]);
  if ( numQueries<1 && (numQueries>=100) && ( (numQueries % 100) !=0 ) && ( numQueries >1000 ) ) {
    printf("invalid number of queries\n");
    return 3;
  }
  int k = atoi(argv[4]);
  if (k<1 || k>8){
    printf("invalid number of nearest neighbours\n");
    return 4;
  }

  char *dataset_file = choose_data_file(numObjects);
  char *query_file = choose_query_file(numQueries);



  printf("name file dataset: %s \n", dataset_file);
  printf("name file queries: %s \n",query_file );
  
  printf("objects: %d\n", numObjects);
  printf("dimentions: %d\n", numAttributes);
  printf("queries: %d\n", numQueries);
  printf("k: %d\n", k);
  
  char *file1 = "NN_distances_cuda.bin";
  char *file2 = "NN_indexes_cuda.bin";

  knn_struct training_set;
  knn_struct query_set;
  float *NNdist;
  int *NNidx;

  training_set.leading_dim = numAttributes;
  training_set.secondary_dim = numObjects;
  query_set.leading_dim = numAttributes;
  query_set.secondary_dim = numQueries;

  /*======== Memory allocation ======*/
  training_set.dataset = (float*)malloc(numObjects*numAttributes*sizeof(float));
  query_set.dataset = (float*)malloc(numQueries*numAttributes*sizeof(float));
  NNdist = (float*)malloc(numQueries*k*sizeof(float));
  NNidx = (int*)malloc(numQueries*k*sizeof(int));



  int i,j;


//===file data====//

  FILE *fp;
  float* data = (float*)malloc(numObjects*numAttributes*sizeof(float));

  fp = fopen( dataset_file , "rb");
  if(fp==NULL){printf("Error opening the file\n");}

  int w = fread(data, sizeof(float), numObjects*numAttributes, fp);
  if(w!=numObjects*numAttributes){printf("Error reading the data\n");}
  
  
 
  for (i=0 ; i<numObjects ; i++){
  	for (j=0 ; j< numAttributes ; j++ ){
  		//z= data[ i*numAttributes + j];
  		training_set.dataset[ i*numAttributes + j] = data[ i*numAttributes + j];
  	}
  	
  }
  
  fclose(fp);
  free(data);



//===file query ====//
FILE *fpq;
  float* queries = (float*)malloc(numQueries*numAttributes*sizeof(float));

  fpq = fopen(query_file, "rb");
  if(fpq==NULL){printf("Error opening the file\n");}

  int w2 = fread(queries, sizeof(float), numQueries*numAttributes, fpq);
  if(w2!=numQueries*numAttributes){printf("Error reading the data\n");}

  for (i=0 ; i<numQueries ; i++){
  		  for (j=0 ; j<numAttributes ; j++){
		query_set.dataset[ i*numAttributes + j] = queries[ i*numAttributes + j];
	  	
	  }
  }

  fclose(fpq);
  free(queries);
  /*===== Cuda Events===*/

  float elapsedTime;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  
  // ==== device stuff ==== //
  cuMemGetInfo_v2(&memory_free, &memory_total);
  printf("Total memory: %zd, free memory: %zd\n", memory_total, memory_free);
  
  knn_struct d_training_set;
  knn_struct d_query_set;
  float *d_NNdist;
  int *d_NNidx;
  float *d_distances;

  d_training_set.leading_dim = numAttributes;
  d_training_set.secondary_dim = numObjects;
  d_query_set.leading_dim = numAttributes;
  d_query_set.secondary_dim = numQueries;
  
 
  /*========= device memory allocation======*/

  cudaMalloc((void**)&d_training_set.dataset, d_training_set.leading_dim*d_training_set.secondary_dim*sizeof(float));

  cudaMalloc((void**)&d_query_set.dataset, d_query_set.leading_dim*d_query_set.secondary_dim*sizeof(float));

  cudaMalloc((void**)&d_NNdist , d_query_set.secondary_dim*k*sizeof(float));
  
  cudaMalloc((void**)&d_NNidx , d_query_set.secondary_dim*k*sizeof(int));
  
  
  cudaMalloc((void**)&d_distances, d_training_set.secondary_dim*sizeof(float));
  
 

  cudaMemcpy(d_training_set.dataset, training_set.dataset, training_set.leading_dim*training_set.secondary_dim*sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(d_query_set.dataset, query_set.dataset, query_set.leading_dim*query_set.secondary_dim*sizeof(float), cudaMemcpyHostToDevice);
  


 cudaEventRecord(start, 0);
  

  knns(&d_training_set, &d_query_set, d_NNdist, d_NNidx, k, d_distances);
  
  
  cudaMemcpy(NNdist, d_NNdist, k*query_set.secondary_dim*sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(NNidx, d_NNidx, k*query_set.secondary_dim*sizeof(int), cudaMemcpyDeviceToHost);
  
 
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);


  cudaEventElapsedTime(&elapsedTime, start, stop);

#ifdef TIMEONLY
  printf("Time elapsed: %f ms\n", elapsedTime);
#endif

   
  printf("Time elapsed: %f ms\n", elapsedTime);
  
  
  gettimeofday(&second, &tzp);


  if(first.tv_usec>second.tv_usec){
    second.tv_usec += 1000000;
    second.tv_sec--;
  }
  
  lapsed.tv_usec = second.tv_usec - first.tv_usec;
  lapsed.tv_sec = second.tv_sec - first.tv_sec;

  printf("Time elapsed: %d.%06dsec\n", lapsed.tv_sec, lapsed.tv_usec); 

/*
 for (i=0;i<numQueries;i++){
	for (j=0; j<k; j++){
		printf("apostash %d query apo ton geirona %d = %f\n",i,NNidx[i*k + j],NNdist[i*k + j] );
	}
}
*/
  /*========save data============*/
  save_distances(NNdist, file1, numQueries, k);
  save_indexes(NNidx, file2, numQueries, k);
  

  /*==== clean device===*/
  cleanDevice(&d_training_set);
  cleanDevice(&d_query_set);


  cudaEventDestroy(start);
  cudaEventDestroy(stop);


  cudaDeviceReset();
  
  return 1;
}
  
  
  