#include "matrix.h"

void setMatrices(double * &MA, double * &MB, double * &MC, int &rows, int &cols)
{
	rows = 3;
	cols = 3;

	MA = new double[rows * cols];
	MB = new double[rows * cols];
	MC = new double[rows * cols];

	memset(MC, 0, sizeof(*MA) * rows * cols);

	int count = rows * cols;

	for (int i = 0; i < count; i++)
	{
		MA[i] = i + 1;
		MB[i] = i % 3;
	}
}
void shuffleMatrices(double *&MA, double *&MB, int rows, int cols, int size)
{
	int subMatrixElems = (rows * cols) / size;

	// horizontal swap
	for (int i = 1; i < rows; i++)
	{
		int pivot = i * cols;

		for (int j = 0; j < cols; j++)
		{
			int targetCol = (j - i + cols) % cols;

			std::swap(MA[pivot], MA[i * cols + targetCol]);
			std::swap(MB[pivot], MB[i * cols + targetCol]);
		}
	}
}

void printMatrix(double *matrix, int rows, int cols)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			std::cout << matrix[i * cols + j] << " ";
		}

		std::cout << std::endl;
	}
}