#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* Program Parameters */
#define BILLION 1000000000L
#define MAXN 3000 /* Max value of N */
int N;            /* Matrix size */

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
        for (col = 0; col < N; col++)
        {
            fscanf(inputfp, "%lf", &A[row][col]);
        }
        fscanf(inputfp, "%lf", &B[row]);
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
    double temp;

    /* Gaussian elimination */
    for (i = 0; i < N - 1; i++)
    {
        for (j = i + 1; j < N; j++)
        {
            temp = A[j][i] / A[i][i];
            for (k = i + 1; k < N; k++)
            {
                A[j][k] = A[j][k] - A[i][k] * temp;
            }
            B[j] = B[j] - B[i] * temp;
        }
    }

    /* Back substitution */
    for (j = N - 1; j >= 0; j--)
    {
        X[j] = B[j];
        for (k = N - 1; k > j; k--)
        {
            X[j] = X[j] - A[j][k] * X[k];
        }
        X[j] = X[j] / A[j][j];
    }
}

int main(int argc, char *argv[])
{
    struct timespec start, end;
    double elapsed_time, elapsed_time_us, elapsed_time_ms, elapsed_time_s;

    // RandMatrixVectorGen(argc, argv);
    ReadInputAndInitializeMatrix(argc, argv);
    print_inputs();
    clock_gettime(CLOCK_MONOTONIC, &start);
    gauss();
    clock_gettime(CLOCK_MONOTONIC, &end);
   // elapsed_time = (end.tv_sec - start.tv_sec) * BILLION + end.tv_nsec - start.tv_nsec;
    elapsed_time_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
    elapsed_time_ms = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_nsec - start.tv_nsec) / 1000000;
    elapsed_time_s = (end.tv_sec - start.tv_sec)  + (end.tv_nsec - start.tv_nsec) / 1000000000;
    print_X();

    //printf("gauss runtime: %llu nanoseconds\n", (long long unsigned int)elapsed_time);
    printf("the input matrix size: %d\n", N);
    printf("gauss runtime in microseconds: %llu us\n", (long long unsigned int)elapsed_time_us);
    printf("gauss runtime in milliseconds: %llu ms\n", (long long unsigned int)elapsed_time_ms);
    printf("gauss runtime in seconds: %llu s\n", (long long unsigned int)elapsed_time_s);

    return 0;
}

