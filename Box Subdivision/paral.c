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
int tid;
float akmi;
float limits[6]; // mporei kai 8 ?_- !![]_[]
float center[3];
} Box;

Box box[]; // array of all boxes
Box leaf[]; // array of non-null leaf boxes
float **A; //array with all the elements
float **B; // array of elements (sorted) after the execution
static int boxcounter;
static int leafcounter;
static int NUMBEROFPOINTS;
static int S;
static int NextBstarts;

/*
#ifdef RAND_MAX
#undef RAND_MAX
#define RAND_MAX [10000]
#endif
*/ 


//#define MAX_THREAD 1000


/*
	Produces a random int x, min <= x <= max 
	following a uniform distribution
*/
/*float randomFloat(int min, int max)
{
 float a;

a=min + rand() % (max - min + 1);
printf("%f \n",a);
return a; 
}
*/


//@@@ dimiourgia kuvou
void *cube( void *boxParent){
	int lcounter;
	//mesa se ena mutex!
	boxcounter++;
	lcounter=boxcounter;
	Box newbox;
	Box[boxcounter]=newbox;
	// telos mutex
	Box[lcounter]->boxid=lcounter; //boxid .. mporw na to balw kai pio katw gia na mh to ksanaallaksw
	Box[lcounter]->parent=boxParent->boxid;   //parent
	int i,j,k;
	float limitsParent[6];
	for (i=0; i<5 ;i++){
		limitsParent[i]=((Box *) boxParent )->limits[i];
	}
	int type=((Box *) box)->tid;		
	Box[lcounter]->type=type;					//type
	float akmi=((Box *) box )->akmi/2; 
	Box[lcounter]->akmi=akmi; 					//akmi
	float *lim=findLimits(limits,type,akmi);
	
	float *centerPoints=findCenter(lim[0],lim[2],lim[4],akmi);
	Box[lcounter]->center[0]=centerPoints[0];       
	Box[lcounter]->center[1]=centerPoints[1];				//center[3]
	Box[lcounter]->center[2]=centerPoints[2];

	for (j=0;j<5;j++){
		Box[lcounter]->limits[j]=lim[j];  //limits[6]
	}
	/*
	xmin=lim[0];
	xmax=lim[1];
	ymin=lim[2];
	ymax=lim[3];
	zmin=lim[4];
	zmax=lim[5];
	*/
	float *points;
	points = searchForPoints(lim[0],lim[1],lim[2],lim[3],lim[4],lim[5],NUMBEROFPOINTS,A,S);
	Box[lcounter]->n= points[NUMBEROFPOINTS+1]; //n  passarw ws teleutaio argument tou pointer to n!
	

	
	//gia leaf adeio
	if (n==0){
		box[lcounter]->boxid=0;

		//gia leaf non-empty
	}else if (n<=S){


		//setarisma tou non-empty leaf
		//se mutex2
			Box[lcounter]->start=NextBstarts;   //start
			ThisBstarts=NextBstarts;
			NextBstarts=NextBstarts+n;
		//telosmutex



		for (k=0;k<n;k++){
			B[0][k+ThisBstarts]= A[0][points[k]];
			B[1][k+ThisBstarts]= A[1][points[k]];
			B[2][k+ThisBstarts]= A[2][points[k]];
		}

		
		l
	}else{ //  gia n > S
		//setarisma 8 newn childZ
		int childId;
		levelOfParent=((Box *) box )->level;
		
		for (j=1;j<8;j++){
		childId=8*levelOfParent+j;
		pthread_create(&thread[childId],NULL,cube,//dat[j*8+1]);
        ptread_join(thread[childId],NULL);
        
        }

	}

}

//@@@ epistrefei tis suntetagmenes twn koryfwn tou kuvou 8 h 6 ? :S
 float * findLimits(float limitsOfParent[],int type,float akmi){
		float  limitsOfChild[6];
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

//@@@ psaxnei gia points kai ta apothikeuei ean dimiourgithei o kuvos
//args megethosA,A,To S pou dinoume apo konsola
float *searchForPoints(float xmin,float xmax,float ymin,float ymax,float zmin,float zmax,int NUM,float **A,int S){

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

//@@@  root creation
void createRootsChilren(void *rootData){

	int i;
	for (i=1;i<8;i++){
	Box[i]= struct Box { .tid=i };
	//what kid dld//!
	pthread_create(&thread[i],NULL,cube,(void *)rootData);
    pthread_join(thread[i],NULL);
    
}


//@@@ generates ta shmeia panw sti unitary sphere
float * generatePointsOnSphere(){

static float p[3];
float x,y,z;
float *points;


x=(float)rand()/RAND_MAX;
y=(float)rand()/RAND_MAX;
z=(float)rand()/RAND_MAX;


float length = (float)sqrt((x * x) + (y * y) + (z * z)) ;

p[0]=x/length;
p[1]=y/length;
p[2]=z/length;

return points;

}
//@@@ telos_generatePointsOnSphere




//@@@ MAIN
int main(int argc, char* argv[]) {  //na doso kai S kai nt gia dokimes


unsigned long int n,i,seedz,nt;
n=atoi(argv[1]);
NUMBEROFPOINTS=n;
S=atoi(argv[2]);
nt=atoi(argv[3]);

pthread_t *threads;
if((nt<1 || nt>MAX_THREAD)){
	printf("Enter valid number of threads, between 1 and %d !\n", MAX_THREAD);
	exit (1);
}
threads = (pthread_t *)malloc(nt*sizeof(pthread_t));

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
struct Box root { .boxid=1 , .level=0 , .n=n , .center[0]=0.5 , .center[1]=0.5 , .center[2]=0.5 , .akmi=1.0 , .tid =0 };


/*root->boxid=1;
root->level=0;
root->n=n;
root->parent=0;
root->center[0]=0.5;
root->center[1]=0.5;
root->center[2]=0.5;
root->akmi=1;
root->tid=0;
*/

for (i=0; i<7 ;i++)
{
	root->child[i]=i+2;
}
/*
for (i=0 ; i<27 ;i++)
{
	root->colleague[i]=0;     //den exei colleague h riza
}
*/

box[0]=root; 				  //1o antikeimeno to root

//@@@ telos_3

pthread_create(&threads[0],NULL,createRootsChilren,(void *)&box[0]);


//@@@ 4) dimiourgia dentrou
















}

free(A);
printf("TELOC omg\n");
exit(1);
return 1;
}


