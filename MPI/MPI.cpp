#include "mpi.h"
#include <iostream>

#include "matrix.h"

int main(int argc, char argv[])
{
	MPI_Init(NULL, NULL);
	
	int rank;
	int size;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	double *MA = NULL;
	double *MB = NULL;
	double *MC = NULL;

	if (rank == 0)
	{
		int rows;
		int cols;

		setMatrices(MA, MB, MC, rows, cols);
		printMatrix(MA, rows, cols);
		
		shuffleMatrices(MA, MB, rows, cols, size);

		printMatrix(MA, rows, cols);
	}
	else
	{
		
	}

	int test;
	scanf("%d", &test);
	
	MPI_Finalize();

	return 0;
}

