#include "mpi.h"
#include <cmath>
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

	int rows;
	int cols;

	if (rank == 0)
	{
		setMatrices(MA, MB, MC, rows, cols);
		shuffleMatrices(MA, MB, rows, cols, size);

		printMatrix(MA, rows, cols);

		initialSendMatrices(MA, MB, rows, cols, size);
	}
	else
	{
		receiveAndSetDimensions(MA, MB, MC, rows, cols);
		initialReceiveMatrices(MA, MB, rows, cols, rank, size);
	}

	delete[] MA;
	delete[] MB;
	delete[] MC;

	int test;
	scanf("%d", &test);
	
	MPI_Finalize();

	return 0;
}

