﻿#include <mpi.h>
#include <cmath>
#include <iostream>

#include "matrix.h"

using std::cout;
using std::endl;

const bool PRINT_MATRICES = true;

int main(int argc, char **argv)
{
	MPI_Init(NULL, NULL);
	
	int rank;
	int size;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	double *MA = NULL;
	double *MB = NULL;
	double *MC = NULL;

	int sub_rows;
	int sub_cols;
	int subMatrixElems;

	int procWidth = (int) sqrt(size);

	if (rank == 0)
	{
		double *FullMatrixA = NULL;
		double *FullMatrixB = NULL;
		double *FullMatrixC = NULL;

		int full_rows;
		int full_cols;
		
		// initialize given and result matrices
		initializeMatrices(FullMatrixA, FullMatrixB, FullMatrixC, full_rows, full_cols);

		if (PRINT_MATRICES)
		{
			printMatrix(FullMatrixA, full_rows, full_cols);
			printMatrix(FullMatrixB, full_rows, full_cols);
		}

		shuffleMatrices(FullMatrixA, FullMatrixB, full_rows, full_cols, size);

		// initialize submatrices
		sub_rows = sub_cols = (int) sqrt((full_rows * full_cols) / size);
		subMatrixElems = sub_rows * sub_cols;
		
		setMatrices(MA, MB, MC, sub_rows, sub_rows);

		initialSendMatrices(FullMatrixA, FullMatrixB, full_rows, full_cols, size);

		int subMatrixIndex = 0;

		for (int i = 0; i < sub_rows; i++)
		{
			for (int j = 0; j < sub_cols; j++)
			{
				MA[subMatrixIndex] = FullMatrixA[i * full_cols + j];
				MB[subMatrixIndex++] = FullMatrixB[i * full_cols + j];
			}
		}

		// counting phase
		for (int i = 0; i < procWidth; i++)
		{
			multiplyMatrices(MA, MB, MC, sub_rows, sub_cols);

			moveMatrix(MA, subMatrixElems, DIRECTION_LEFT, TAG_MATRIX_A, rank, size);
			moveMatrix(MB, subMatrixElems, DIRECTION_UP, TAG_MATRIX_B, rank, size);

			receiveMatrices(MA, MB, subMatrixElems, rank, size);
		}

		subMatrixIndex = 0;

		for (int i = 0; i < sub_rows; i++)
		{
			for (int j = 0; j < sub_cols; j++)
			{
				FullMatrixC[i * full_cols + j] = MC[subMatrixIndex++];
			}
		}

		for (int i = 0; i < size - 1; i++)
		{
			collectResult(FullMatrixC, full_rows, full_cols, size);
		}

		if (PRINT_MATRICES)
		{
			printMatrix(FullMatrixC, full_rows, full_cols);
		}

		delete[] FullMatrixA;
		delete[] FullMatrixB;
		delete[] FullMatrixC;
	}
	else
	{
		receiveDimensions(sub_rows, sub_cols);

		setMatrices(MA, MB, MC, sub_rows, sub_cols);	// bottleneck?

		initialReceiveMatrices(MA, MB, sub_rows * sub_cols, rank, size);

		subMatrixElems = sub_rows * sub_cols;

		for (int i = 0; i < procWidth; i++)
		{
			multiplyMatrices(MA, MB, MC, sub_rows, sub_cols);

			moveMatrix(MA, subMatrixElems, DIRECTION_LEFT, TAG_MATRIX_A, rank, size);
			moveMatrix(MB, subMatrixElems, DIRECTION_UP, TAG_MATRIX_B, rank, size);

			receiveMatrices(MA, MB, subMatrixElems, rank, size);
		}

		sendResult(MC, subMatrixElems, 0);
	}

	delete[] MA;
	delete[] MB;
	delete[] MC;
	
	MPI_Finalize();

	return 0;
}

