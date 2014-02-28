#include "matrix.h"

void setMatrices(double * &MA, double * &MB, double * &MC, int &rows, int &cols)
{
	rows = 4;
	cols = 4;

	MA = new double[rows * cols];
	MB = new double[rows * cols];
	MC = new double[rows * cols];

	int count = rows * cols;

	for (int i = 0; i < count; i++)
	{
		MA[i] = i + 1;
		MB[i] = i % 3;
	}
}
void shuffleMatrices(double *MA, double *MB, int rows, int cols, int size)
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

	// vertical swap
	for (int j = 1; j < cols; j++)
	{
		int pivot = j;

		for (int i = 0; i < rows; i++)
		{
			int targetRow = (i - j + rows) % rows;

			std::swap(MA[pivot], MA[targetRow * cols + j]);
			std::swap(MB[pivot], MB[targetRow * cols + j]);
		}
	}
}

void initialSendMatrices(const double * MA, const double * MB, int rows, int cols, int size)
{
	int subMatrixElems = (rows * cols) / size;	// number of elements in submatrix
	int subMatrixDim = subMatrixElems / 2;		// dimension of submatrix
	
	double *bufferMA = new double[subMatrixElems];
	double *bufferMB = new double[subMatrixElems];
	int bufferPos = 0;

	int dimensions[2] = { rows, cols };

	int procWidth = sqrt(size);	// dimension of grid

	int i;	// processor x coord in grid
	int j;	// processor y coord in grid

	for (int p = 1; p < size; p++)
	{
		MPI_Send(dimensions, 2, MPI_INT, p, TAG_MATRIX_DIMENSIONS, MPI_COMM_WORLD);

		bufferPos = 0;

		i = p / procWidth; 
		j = p % procWidth; 

		int xStart = i * subMatrixDim;
		int xEnd = xStart + subMatrixDim;

		int yStart = j * subMatrixDim;
		int yEnd = yStart + subMatrixDim;		

		for (int x = xStart; x < xEnd; x++)
		{
			for (int y = yStart; y < yEnd; y++)
			{
				bufferMA[bufferPos] = MA[x * cols + y];
				bufferMB[bufferPos++] = MB[x * cols + y];
			}
		}

		MPI_Send(bufferMA, subMatrixElems, MPI_DOUBLE, p, TAG_MATRIX_A, MPI_COMM_WORLD);
		MPI_Send(bufferMB, subMatrixElems, MPI_DOUBLE, p, TAG_MATRIX_B, MPI_COMM_WORLD);
	}

	delete[] bufferMA;
	delete[] bufferMB;
}
void receiveAndSetDimensions(double * &MA, double * &MB, double * &MC, int &rows, int &cols, int size)
{
	int dimensions[2];

	MPI_Recv(dimensions, 2, MPI_INT, 0, TAG_MATRIX_DIMENSIONS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	rows = dimensions[0];
	cols = dimensions[1];

	int subMatrixElems = (rows * cols) / size;	// number of elements in submatrix

	MA = new double[subMatrixElems];
	MB = new double[subMatrixElems];
	MC = new double[subMatrixElems];
}
void initialReceiveMatrices(double * MA, double * MB, int rows, int cols, int rank, int size)
{
	int subMatrixElems = (rows * cols) / size;	// number of elements in submatrix

	MPI_Status status;

	double *buffer = new double[subMatrixElems];

	for (int i = 0; i < 2; i++)
	{
		MPI_Recv(buffer, subMatrixElems, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		if (status.MPI_TAG == TAG_MATRIX_A)
		{
			memcpy(MA, buffer, sizeof(*buffer) * subMatrixElems);
		}
		else memcpy(MB, buffer, sizeof(*buffer) * subMatrixElems);
	}

	delete[] buffer;
}

void multiplyMatrices(const double *MA, const double *MB, double *MC, int rows, int cols)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			int position = i * cols + j;
			MC[position] = 0;

			for (int k = 0; k < cols; k++)
			{
				MC[position] += MA[i * cols + k] * MB[k * cols + j];
			}
		}
	}
}

void printMatrix(double *matrix, int rows, int cols)
{
	std::cout << std::endl;

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			std::cout << matrix[i * cols + j] << " ";
		}

		std::cout << std::endl;
	}
}