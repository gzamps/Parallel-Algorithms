#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/types.h>
#include <math.h>
#include <pthread.h>
#define MAX_THREAD 1000

//@@@ structure with box info
typedef struct Box {
int level, boxid, parent, child[8],
n, start,
colleague[26];
char type;
int tid,id;
float akmi;
float limits[6]; // mporei kai 8 ?_- !![]_[]
float center[3];
} Box;

Box box[]; // array of all boxes
Box leaf[]; // array of non-null leaf boxes
float **A; //array with all the elements
float **B; // array of elements (sorted) after the execution

static int boxcounter=0;
static int leafcounter;
static int NUMBEROFPOINTS;
static int S;
static int NextBstarts;
static int levelWeWork;
static int l;


//@@@ dimiourgia kuvou
void *cube( boxParent){
	int localcounter;
	//mesa se ena mutex!
	boxcounter++;
	localcounter=boxcounter;
	Box newbox;
	box[boxcounter]=newbox;
	// telos mutex
	//newbox.boxid=localcounter;
	//box[localcounter].boxid=localcounter; //boxid .. mporw na to balw kai pio katw gia na mh to ksanaallaksw
	box[localcounter].parent=boxParent.boxid;   //parent
	int i,j,k;
	float limitsParent[6];
	for (i=0; i<5 ;i++){
		limitsParent[i]=boxParent.limits[i];
	}
	/*
	int type=box.tid;		
	box[localcounter].type=type;					//type
	float akmi=box.akmi/2; 
	box[localcounter].akmi=akmi; 					//akmi
	float *lim=findLimits(limits,type,akmi);
	*/
	
	int type=localcounter % 8;
	box[localcounter].type=type;
	float akmi=boxParent.akmi/2;
	float* lim=findLimits(limitsParent,type,akmi);

	float *centerPoints=findCenter(lim[0],lim[2],lim[4],akmi);
	box[localcounter].center[0]=centerPoints[0];       
	box[localcounter].center[1]=centerPoints[1];				//center[3]
	box[localcounter].center[2]=centerPoints[2];
	free (centerpoints);
	for (j=0;j<5;j++){
		box[localcounter].limits[j]=lim[j];  //limits[6]
	}
	/*
	xmin=lim[0];
	xmax=lim[1];
	ymin=lim[2];
	ymax=lim[3];
	zmin=lim[4];
	zmax=lim[5];
	*/
	if ((localcounter % 8 )==1){
		startpoint=boxParent.start;
		//ti endpoint
	}else{
		startpoint=box[localcounter-1].start + box[localcounter-1].n;
	}
	

	if (boxcounter<9){
	//box[localcounter].n=searchForPoints(lim[0],lim[1],lim[2],lim[3],lim[4],lim[5],A);
	int *indexes=searchForPoints(lim[0],lim[1],lim[2],lim[3],lim[4],lim[5],A);
	}else {
	//box[localcounter].n=searchForPoints(lim[0],lim[1],lim[2],lim[3],lim[4],lim[5],B,startpoint);//endpoint);
	int *indexes=searchForPoints(lim[0],lim[1],lim[2],lim[3],lim[4],lim[5],B,startpoint);
	}

	free(lim);

	
	//gia leaf adeio
	if (box[localcounter].n==0){
		box[localcounter].boxid=0;
		boxParent.child[localcounter%8]=0; //miden to boxid pou antistoixei sto paidi ;)
		box[localcounter].start=startpoint;


		//gia leaf non-empty
	}else if (box[localcounter].n<=S){

		/*
		//setarisma tou non-empty leaf
		//se mutex2
			box[localcounter]->start=NextBstarts;   //start
			ThisBstarts=NextBstarts;
			NextBstarts=NextBstarts+n;
		//telosmutex
		*/
		boxParent.child[localcounter%8]=box[localcounter].boxid; //perasma to boxid pou antistoixei sto paidi ;)
		box[localcounter].start=startpoint;
		allocatePoints(indexes,  box[localcounter].start,  box[localcounter].n);
		free(indexes);
		/*
		for (k=0;k<n;k++){
			B[0][k+ThisBstarts]= A[0][points[k]];
			B[1][k+ThisBstarts]= A[1][points[k]];
			B[2][k+ThisBstarts]= A[2][points[k]];
		}
	*/
		
		l
	}else{ //  gia n > S
		/*//setarisma 8 newn childZ
		int childId;
		levelOfParent=boxParent.level;
		
		for (j=1;j<8;j++){
		childId=8*levelOfParent+j;
		
        */

       boxParent.child[localcounter%8]=box[localcounter].boxid; //perasma to boxid pou antistoixei sto paidi ;)
        
        }

	
	return;
}


//@@@ epistrefei tis suntetagmenes twn koryfwn tou kuvou 8 h 6 ? :S
 float * findLimits(float limitsOfParent[],int type,float akmi){
		float * limitsOfChild=malloc(6*sizeof(float));
	switch (type){

		case '1': //xyz
				limitsOfChild[0]=limitsOfParent[0]+akmi/2; //xmin
				limitsOfChild[1]=limitsOfParent[1];		   //xmax
				limitsOfChild[2]=limitsOfParent[2]+akmi/2; //ymin
				limitsOfChild[3]=limitsOfParent[3];		   //ymax
				limitsOfChild[4]=limitsOfParent[4]+akmi/2; //zmin
				limitsOfChild[5]=limitsOfParent[5];	  	  //zmax
		case '2': //-xyz
				limitsOfChild[0]=limitsOfParent[0]; 		//xmin
				limitsOfChild[1]=limitsOfParent[1]-akmi/2;		//xmax
				limitsOfChild[2]=limitsOfParent[2]+akmi/2; //ymin
				limitsOfChild[3]=limitsOfParent[3];		   //ymax
				limitsOfChild[4]=limitsOfParent[4]+akmi/2; //zmin
				limitsOfChild[5]=limitsOfParent[5];	  	   //zmax
		case '3': //-x-yz
				limitsOfChild[0]=limitsOfParent[0]; 		 //xmin
				limitsOfChild[1]=limitsOfParent[1]-akmi/2;		 //xmax
				limitsOfChild[2]=limitsOfParent[2]; 		 //ymin
				limitsOfChild[3]=limitsOfParent[3]-akmi/2;		 //ymax
				limitsOfChild[4]=limitsOfParent[4]+akmi/2;   //zmin
				limitsOfChild[5]=limitsOfParent[5];	  		 //zmax
		case '4': //-xy-z
				limitsOfChild[0]=limitsOfParent[0]; 		//xmin
				limitsOfChild[1]=limitsOfParent[1]-akmi/2;	   //xmax
				limitsOfChild[2]=limitsOfParent[2]+akmi/2; //ymin
				limitsOfChild[3]=limitsOfParent[3];		   //ymax
				limitsOfChild[4]=limitsOfParent[4]; 		//zmin
				limitsOfChild[5]=limitsOfParent[5]-akmi/2;	   //zmax
		case '5': //x-yz
				limitsOfChild[0]=limitsOfParent[0]+akmi/2; //xmin
				limitsOfChild[1]=limitsOfParent[1];		   //xmax
				limitsOfChild[2]=limitsOfParent[2]; 		//ymin
				limitsOfChild[3]=limitsOfParent[3]-akmi/2;	   //ymax
				limitsOfChild[4]=limitsOfParent[4]+akmi/2; //zmin
				limitsOfChild[5]=limitsOfParent[5];	      //zmax
		case '6': //x-y-z
				limitsOfChild[0]=limitsOfParent[0]+akmi/2; //xmin
				limitsOfChild[1]=limitsOfParent[1];		   //xmax
				limitsOfChild[2]=limitsOfParent[2]; 		//ymin
				limitsOfChild[3]=limitsOfParent[3]-akmi/2;	    //ymax
				limitsOfChild[4]=limitsOfParent[4]; 		//zmin
				limitsOfChild[5]=limitsOfParent[5]-akmi/2;	   //zmax
		case '7': //xy-z
				limitsOfChild[0]=limitsOfParent[0]+akmi/2; //xmin
				limitsOfChild[1]=limitsOfParent[1];		   //xmax
				limitsOfChild[2]=limitsOfParent[2]+akmi/2; //ymin
				limitsOfChild[3]=limitsOfParent[3];		   //ymax
				limitsOfChild[4]=limitsOfParent[4]; 		//zmin
				limitsOfChild[5]=limitsOfParent[5]-akmi/2;	  	 //zmax
		case '8': //-x-y-z
				limitsOfChild[0]=limitsOfParent[0]; 		//xmin
				limitsOfChild[1]=limitsOfParent[1]-akmi/2;	    //xmax
				limitsOfChild[2]=limitsOfParent[2];		    //ymin
				limitsOfChild[3]=limitsOfParent[3]-akmi/2;	    //ymax
				limitsOfChild[4]=limitsOfParent[4]; 		//zmin
				limitsOfChild[5]=limitsOfParent[5]-akmi/2;	   //zmax

	}
	return limitsOfChild;
}

float *findCenter(float xmin,float ymin,float zmin,float akmi){
	float *c;
	c[0]=xmin+akmi/2;
	c[1]=ymin+akmi/2;
	c[2]=zmin+akmi/2;

	return c;
}
		/*
//@@@ psaxnei gia points kai ta apothikeuei ean dimiourgithei o kuvos
//args megethosA,A,To S pou dinoume apo konsola
float searchForPoints(float xmin,float xmax,float ymin,float ymax,float zmin,float zmax,float **A){

	float *indexes;
	int i;
	int numOfPointsInCube=0;

	for (i=0;i<NUM;i++ ){
		if ( (A[0][i] > xmin )&&( A[0][i] < xmax) &&(A[1][i] > ymin)&&(A[1][i] < ymax)&&(A[2][i] > zmin)&&(A[2][i] < zmax){
			numOfPointsInCube++;
			indexes[i]=i;
		}
	}
	indexes[i+1]=numOfPointsInCube;
	return indexes;

}

*/
int* searchForPoints(float xmin,float xmax,float ymin,float ymax,float zmin,float zmax,float **M,int startofpoints, int numberofpts){

	int *indexes=malloc(numberofpts*sizeof(int));
	int i;
	int numOfPointsInCube=0;

	for (i=startofpoints;i<numberofpts;i++ ){
		if ( (M[0][i] > xmin )&&( M[0][i] < xmax) &&(M[1][i] > ymin) && (M[1][i] < ymax) && (M[2][i] > zmin) && (M[2][i] < zmax)){
			numOfPointsInCube++;
			indexes[i]=i;
			
		}
	}
	box[boxcounter].n=numOfPointsInCube;
	return indexes;

}
/*
//@@@  root creation
void createRootsChilren(void *rootData){

	int i;
	for (i=1;i<8;i++){
	Box[i]= struct Box { .tid=i };
	//what kid dld//!
	pthread_create(&thread[i],NULL,cube,(void *)rootData);
    pthread_join(thread[i],NULL);
    
}
*/
//@@@
void allocatePoints(int *index,int startpoint,int n){
	int i;
	int j=0;
	//int *indexes = malloc(n * sizeof(int));
	for (i=startpoint;i<startpoint+n;i++){
		B[0][i]=B[0][index[j]];
		B[1][i]=B[1][index[j]];
		B[2][i]=B[2][index[j]];
		j++;

	}
	return;
}
//@@@ generates ta shmeia panw sti unitary sphere
float *generatePointsOnSphere(){

static float p[3];
float x,y,z;



x=(float)rand()/RAND_MAX;
y=(float)rand()/RAND_MAX;
z=(float)rand()/RAND_MAX;


float length = (float)sqrt((x * x) + (y * y) + (z * z)) ;

p[0]=x/length;
p[1]=y/length;
p[2]=z/length;

return p;

}
//@@@ telos_generatePointsOnSphere




//@@@ MAIN
void main(int argc, char* argv[]) {  //na doso kai S kai nt gia dokimes


unsigned long int n,i,seedz,nt;
n=atoi(argv[1]);
NUMBEROFPOINTS=n;
S=atoi(argv[2]);
nt=atoi(argv[3]);

/*pthread_t *threads;
if((nt<1 || nt>MAX_THREAD)){
	printf("Enter valid number of threads, between 1 and %d !\n", MAX_THREAD);
	exit (1);
}
threads = (pthread_t *)malloc(nt*sizeof(pthread_t));
*/
i=0;
seedz=0;
srand( (unsigned)time( NULL ));


//@@@ 1A) dilosh tou A
float **A = malloc(3*sizeof(float));

for (i = 0; i < 3; i++) {
  A[i] = malloc(n * sizeof(float));
}
//@@@ 1) telos_1A


//@@@ 1B) dilosh B
float **B = malloc(3*sizeof(float));

for (i = 0; i < 3; i++) {
  B[i] = malloc(n * sizeof(float));
}
//@@@ 1B) telos_1B


//@@@ 2) gemisma tou A
for (i = 0; i < n; i++) {



//printf(" done epanalipsi gia fores = %lu \n", i );

float* points=generatePointsOnSphere();
//??float as=points[0];

A[0][i]=points[0];
A[1][i]=points[1];
A[2][i]=points[2];

printf( "points[0] = %f\n", A[0][i]);
printf( "points[1] = %f\n", A[1][i]);
printf( "points[2] = %f\n", A[2][i]);


}
//@@@ telos_2


//@@@ 3) setup rizas
struct Box root;
root.boxid=1;
root.level=0;
root.n=n;
root.parent=0;
root.center[0]=0.5;
root.center[1]=0.5;
root.center[2]=0.5;
root.akmi=1;
root.tid=0;

for (i=1; i<8 ;i++)
{
	root.child[i]=i+1;
}
/*
for (i=0 ; i<27 ;i++)
{
	root->colleague[i]=0;     //den exei colleague h riza
}
*/

box[0]=root; 				  //1o antikeimeno to root

//@@@ telos_3

//pthread_create(&threads[0],NULL,createRootsChilren,(void *)&box[0]);


//@@@ 4) dimiourgia dentrou





//serial

for (i=1;i<10;i++){
	
	

	if (i<=8){
		cube(box[0]);
		}else if(i>8 && i<=64){
			cube(box[(i-1)/8]);
		}else if(i>64 && i<=8*64){
			cube(box[(i-1)/8]);
			}else if (i>8*64 && i<=64*64){
				cube(box[(i-1)/8]);
				}else if(i>64*64 && i<=64*64*8){
					cube(box[(i-1)/8]);
					}else if(i>64*64*8 && i<=64*64*64){
						cube(box[(i-1)/8]);
						}else if(i>64*64*64 && i<=64*64*64*8){
							cube(box[(i-1)/8]);
							}else {
								cube(box[(i-1)/8]);
							}
						}
					
				

			
			
		











free(A);
free(B);
printf("TELOC omg\n");
exit(1);
}


