#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* Program Parameters */
#define BILLION 1000000000L
#define MAXN 2000 /* Max value of N */
int N;            /* Matrix size */
int seed = 0;     /* Random seed */

/* Matrices and vectors */
double A[MAXN][MAXN], B[MAXN], X[MAXN];

void ReadInputAndInitializeMatrix(int argc, char **argv)
{
    FILE *inputfp = NULL;
    inputfp = fopen(argv[1], "r");
    int row, col;

    /* read number of equations */
    fscanf(inputfp, "%d", &N);
    for (row = 0; row < N; row++)
    {
        for (col = 0; col < N+1; col++)
        {
            fscanf(inputfp, "%lf", &A[row][col]);
        }
        //fscanf(inputfp, "%d", &B[row]);
        X[row] = 0;
    }
    fclose(inputfp);
}

void RandMatrixVectorGen(int argc, char **argv)
{
    srand(time(NULL));
    N = atoi(argv[1]);
    int row, col;

    for (col = 0; col < N; col++)
    {
        for (row = 0; row < N; row++)
        {
            A[row][col] = rand() % N + 1;
        }
        B[col] = rand() % N + 1;
        X[col] = 0;
    }
}

/* Print input matrices */
void print_inputs()
{
    int row, col;

    if (N < 20)
    {
        printf("\nA =\n\t");
        for (row = 0; row < N; row++)
        {
            for (col = 0; col <=N; col++)
            {
                printf("%lf%s", A[row][col], (col < N ) ? ", " : ";\n\t");
            }
        }
        printf("\nB = [");
        for (col = 0; col < N; col++)
        {
            printf("%lf%s", B[col], (col < N - 1) ? "; " : "]\n");
        }
    }
}

void print_X()
{
    int row;
    if (N < 100)
    {
        printf("\nX = [");
        for (row = 0; row < N; row++)
        {
            printf("%lf%s", X[row], (row < N - 1) ? "; " : "]\n");
        }
    }
}

void gauss()
{
    int i, j, k;
    int temp;

    /* Gaussian elimination */
    
    for(i=0;i<N-1;i++){
        //Partial Pivoting
        for(k=i+1;k<N;k++){
            //If diagonal element(absolute vallue) is smaller than any of the terms below it
            if(fabs(A[i][i])<fabs(A[k][i])){
                //Swap the rows
                for(j=0;j<N+1;j++){                
                    double temp;
                    temp=A[i][j];
                    A[i][j]=A[k][j];
                    A[k][j]=temp;
                }
            }
        }
        //Begin Gauss Elimination
        for(k=i+1;k<N;k++){
            register double  tmp=A[k][i];
            for(j=0;j<N+1;j++){
                A[k][j]=A[k][j]-(tmp/ A[i][i])*A[i][j];
            }
        }
         
    }

    printf("\nA =\n\t");
        for (i = 0; i < N; i++)
        {
            for (j = 0; j <=N; j++)
            {
                printf("%lf%s", A[i][j], (j < N ) ? ", " : ";\n\t");
            }
        }
        
    //Begin Back-substitution
    for(i=N-1;i>=0;i--){
        X[i]=A[i][N];
        for(j=i+1;j<N;j++){
            X[i]=X[i]-A[i][j]*X[j];
        }
        X[i]=X[i]/A[i][i];
    }
}

int main(int argc, char *argv[])
{
    struct timespec start, end;
    double elapsed_time;

    // RandMatrixVectorGen(argc, argv);
    ReadInputAndInitializeMatrix(argc, argv);
    print_inputs();
    clock_gettime(CLOCK_MONOTONIC, &start);
    gauss();
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_time = (end.tv_sec - start.tv_sec) * BILLION + end.tv_nsec - start.tv_nsec;
    print_X();

    printf("gauss runtime: %llu nanoseconds\n", (long long unsigned int)elapsed_time);

    return 0;
}

