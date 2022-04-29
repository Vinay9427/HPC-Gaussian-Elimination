#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/times.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h>
#define MAXN 3000// Maximum Dimension for Matrix 

void backSubstitution();
void displayMat();
void printAnswer();
void gaussian_mpi(int N);

int proc,id,N;
static double A[MAXN][MAXN],B[MAXN],X[MAXN];


void ReadInputAndInitializeMatrix(int argc, char **argv)
{
    FILE *inputfp = NULL;
    inputfp = fopen(argv[1], "r");
    int row, col;

    /* read number of equations */
    fscanf(inputfp, "%d", &N);
    for (row = 0; row < N; row++)
    {
        for (col = 0; col < N; col++)
        {
            fscanf(inputfp, "%lf", &A[row][col]);
        }
        fscanf(inputfp, "%lf", &B[row]);
    }
    fclose(inputfp);
}

/* Displays the matrix which has being initialized */
void displayMat()
{	
	int i,j;
	if (N <= 10)
    {
		printf("Displaying Initial Matrix.\n");
		for (i=0;i<N;i++)
		{
			printf("| ");
			for(j=0;j<N; j++)
			{	
				printf("%lf ",A[i][j]);
			}
			printf("| | %lf |\n",B[i]);
		}
	}
}

/* This function performs the backsubstitution */
void backSubstitution()
{
	int i,j;
	for (i=N-1;i>=0;i--)
	{
		int count=0;
		X[i] = B[i];
		for (j=i+1;j<N;j++)
		{
			X[i]-=A[i][j]*X[j];
			
		}
		X[i] = X[i]/A[i][i];
	}
}

/* This function performs gaussian elimination with MPI implementation through static interleave */
void gaussian_mpi(int N)
{	
	double wp_time,wa_time=0;
	MPI_Request request;
	MPI_Status status;
	int p,k,i,j;
	float mp;	

	MPI_Barrier(MPI_COMM_WORLD);// waiting for all processors	
	if(id==0)// Processors starts the MPI Timer i.e MPI_Wtime()
	{
		wa_time = MPI_Wtime();
	}	

	for (k=0;k<N-1;k++)
 	{	
		//Broadcsting X's and Y's matrix from 0th rank processor to all other processors.
		MPI_Bcast(&A[k][0],N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&B[k],1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
		if(id==0)
		{
			for (p=1;p<proc;p++)
			{
		  		for (i=k+1+p;i<N;i+=proc)
		  		{
				/* Sending X and y matrix from oth to all other processors using non blocking send*/
				   MPI_Send(&A[i], N, MPI_DOUBLE, p, 0, MPI_COMM_WORLD);
				   MPI_Send(&B[i], 1, MPI_DOUBLE, p, 0, MPI_COMM_WORLD);
		  		}
			}
			// implementing gaussian elimination 
			for (i=k+1 ; i<N ; i += proc)
			{
	  			mp = A[i][k] / A[k][k];
	  			for (j = k; j < N; j++)
	 			{
	   				A[i][j] -= A[k][j] * mp;
	 			}
	   			B[i] -= B[k] * mp;
			}
			// Receiving all the values that are send by 0th processor.
			for (p = 1; p < proc; p++)
			{
			  for (i = k + 1 + p; i < N; i += proc)
			  {
			    MPI_Recv(&A[i], N, MPI_DOUBLE, p, 1, MPI_COMM_WORLD, &status);
			    MPI_Recv(&B[i], 1, MPI_DOUBLE, p, 1, MPI_COMM_WORLD, &status);
			  }
			}
			//Stopping the MPI_Timer
			if (k == N - 2)
			{
 				wp_time = MPI_Wtime();
 				printf("elapsed time = %f\n", wp_time-wa_time);
			}
		}
		
		
		else
		{
			for (i = k + 1 + id; i < N; i += proc)
			{
				MPI_Recv(&A[i], N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);		
				MPI_Recv(&B[i], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
				mp = A[i][k] / A[k][k];
				for (j = k; j < N; j++)
				{
				    A[i][j] -= A[k][j] * mp;
				}
				B[i] -= B[k] * mp;
				MPI_Send(&A[i], N, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);						    
				MPI_Send(&B[i], 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
			}
		}
		 //MPI_Barrier(MPI_COMM_WORLD);//Waiting for all processors
	}
}


/*  Printing Solution Matrix X   */
void printAnswer()
{
	int i;
	if(N <=20){
		printf("\nSolution Vector (x):\n\n");
		for (i=0;i<N;i++)
		{
			printf("|%lf|\n", X[i]);
		}
	}	
}

int main(int argc,char *argv[])
{
	double start, end;
	MPI_Init(&argc,&argv);//Initiating MPI
	MPI_Comm_rank(MPI_COMM_WORLD,&id);//Getting rank of current processor.
	MPI_Comm_size(MPI_COMM_WORLD,&proc);//Getting number of processor in MPI_COMM_WORLD
		
	
	if (argc >= 2) 
	{
	    N = atoi(argv[2]);//getting matrix dimension from command line argument
   	}


    start = MPI_Wtime();

	if(id==0)
	{
		//initializeMat();//initializing matrix
        ReadInputAndInitializeMatrix(argc, argv);
		displayMat();//displaying the matrix
		/* Start Clock */
 		printf("\nStarting clock.\n");
        start = MPI_Wtime();				
	}
	
	gaussian_mpi(N);//implementing the gaussian elimination
	printf("Size of matrix = %d", N);
	if(id==0)   
	{
		backSubstitution();
        end = MPI_Wtime();	
  		printf("Stopped clock.\n");
		printAnswer();

        printf("\nTime for Execution = %f s.\n", end-start);

	}
	MPI_Finalize(); //Finalizing the MPI
  	return 0;
}