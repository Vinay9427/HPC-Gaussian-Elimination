#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <time.h>
#include <omp.h>

/* Program Parameters */
#define MAXN 3000 /* Max value of N */
int N;            /* Matrix size */

/* Matrices and vectors */
double A[MAXN][MAXN], B[MAXN], X[MAXN];

void gauss(); 

void ReadInputAndInitializeMatrix(int argc, char **argv)
{
    FILE *inputfp = NULL;
    inputfp = fopen(argv[2], "r");
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
        X[row] = 0;
    }
    fclose(inputfp);
}

/* Print input matrices */
void print_inputs()
{
    int row, col;

    if (N < 10)
    {
        printf("\nA =\n\t");
        for (row = 0; row < N; row++)
        {
            for (col = 0; col < N; col++)
            {
                printf("%lf%s", A[row][col], (col < N - 1) ? ", " : ";\n\t");
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

    if (N <= 10)
    {
        printf("\nX = [");
        for (row = 0; row < N; row++)
        {
            printf("%lf%s", X[row], (row < N - 1) ? "; " : "]\n");
        }
    }
}

int main(int argc, char **argv)
{
    /* Timing variables */
    double start, end;
    int num_threads = atoi(argv[1]);
    omp_set_num_threads(num_threads);

    ReadInputAndInitializeMatrix(argc, argv);
    print_inputs();

    start = omp_get_wtime();
    gauss();
    end = omp_get_wtime();

    print_X();

    /* Display timing results */
    printf("\nElapsed time = %f s.\n", end-start);

    exit(0);
}


void gauss()
{
    int norm, row, col; /* Normalization row, and zeroing
                         * element row and col */
    double multiplier;

    printf("Computing using OpenMP.\n");

    /* Gaussian elimination */

    for (norm = 0; norm < N - 1; norm++)
    {
#pragma omp parallel for shared(A, B) private(multiplier, row, col)
        for (row = norm + 1; row < N; row++)
        {
            multiplier = A[row][norm] / A[norm][norm];
            for (col = norm; col < N; col++)
            {
                A[row][col] -= A[norm][col] * multiplier;
            }
            B[row] -= B[norm] * multiplier;
        }
    }
    /* (Diagonal elements are not normalized to 1.  This is treated in back
     * substitution.)
     */

    /* Back substitution */
    for (row = N - 1; row >= 0; row--)
    {
        X[row] = B[row];
        for (col = N - 1; col > row; col--)
        {
            X[row] -= A[row][col] * X[col];
        }
        X[row] /= A[row][row];
    }
}