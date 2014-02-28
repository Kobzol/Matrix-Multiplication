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

	if (rank == 0)
	{
		
	}
	else
	{
		
	}

	
	MPI_Finalize();

	return 0;
}

