#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define BILLION 1000000000L

int A[3][3] = {
    {1, -2, 3},
    {-1, 3, -1},
    {2, -5, 5}};
int B[3] = {9, -6, 17};
int X[3] = {0, 0, 0};
int N = 3;

int RandMatrixGen(double *M, int n)
{
    int i;
    srand(time(NULL));
    for (i = 0; i < n * n; i++)
        M[i] = rand() % n + 1;
    return 0;
}

int RandVectorGen(double *M, int n)
{
    int i;
    srand(time(NULL));
    for (i = 0; i < n; i++)
        M[i] = rand() % n + 1;
    return 0;
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
                printf("%d%s", A[row][col], (col < N - 1) ? ", " : ";\n\t");
            }
        }
        printf("\nB = [");
        for (col = 0; col < N; col++)
        {
            printf("%d%s", B[col], (col < N - 1) ? "; " : "]\n");
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
            printf("%d%s", X[row], (row < N - 1) ? "; " : "]\n");
        }
    }
}

void gauss()
{
    int i, j, k;
    int temp;
    printf("Computing Serially.\n");

    /* Gaussian elimination */
    for (i = 0; i < N - 1; i++)
    {
        for (j = i + 1; j < N; j++)
        {
            temp = A[j][i] / A[i][i];
            for (k = i; k < N; k++)
            {
                A[j][k] -= A[i][k] * temp;
            }
            B[j] -= B[i] * temp;
        }
    }

    printf("\nA =\n\t");
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            printf("%d%s", A[i][j], (j < N - 1) ? ", " : ";\n\t");
        }
    }

    /* Back substitution */
    for (j = N - 1; j >= 0; j--)
    {
        X[j] = B[j];
        for (k = N - 1; k > j; k--)
        {
            X[j] -= A[j][k] * X[k];
        }
        X[j] /= A[j][j];
    }
}

int main(int argc, char *argv[])
{
    struct timespec start, end;
    double elapsed_time, elapsed_time_us;

    print_inputs();
    clock_gettime(CLOCK_MONOTONIC, &start);
    gauss();
    clock_gettime(CLOCK_MONOTONIC, &end);
    //elapsed_time = (end.tv_sec - start.tv_sec) * BILLION + end.tv_nsec - start.tv_nsec;
    elapsed_time_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
    print_X();

    //printf("gauss runtime: %llu nanoseconds\n", (long long unsigned int)elapsed_time);
    printf("gauss runtime: %llu microseconds\n", (long long unsigned int)elapsed_time_us);

    return 0;
}
