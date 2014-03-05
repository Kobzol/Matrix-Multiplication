#include "matrix.h"

void initializeMatrices(double * &MA, double * &MB, double * &MC, int &rows, int &cols)
{
	rows = 4;
	cols = 4;

	setMatrices(MA, MB, MC, rows, cols);

	int count = rows * cols;

	/* Fill matrix with some data */
	for (int i = 0; i < count; i++)
	{
		MA[i] = (i + 1) % 3;
		MB[i] = i % 3;
	}
}

void setMatrices(double * &MA, double * &MB, double * &MC, int rows, int cols)
{
	int elems = rows * cols;

	MA = new double[elems];
	MB = new double[elems];
	MC = new double[elems];

	memset(MC, 0, sizeof(*MC) * elems);
}

void shuffleMatrices(double *MA, double *MB, int rows, int cols, int size)
{
	int sub_rows = (int) sqrt((rows * cols) / size);
	int sub_cols = sub_rows;
	
	int subMatrixElems = sub_rows * sub_cols;

	double *rowBuffer = new double[cols];

	// horizontal swap
	unsigned int gap = 0;

	for (int i = sub_rows; i < rows; i++)
	{
		if (i % sub_rows == 0)
		{
			gap++;
		}

		for (int j = 0; j < cols; j++)
		{
			int targetCol = (j - (gap * sub_cols) + cols) % cols;	// find target column

			rowBuffer[targetCol] = MA[i * cols + j];
		}

		memcpy(MA + i * cols, rowBuffer, sizeof(*rowBuffer) * cols);
	}

	delete[] rowBuffer;

	double *colBuffer = new double[rows];

	// vertical swap
	gap = 0;

	for (int j = sub_cols; j < cols; j++)
	{
		if (j % sub_rows == 0)
		{
			gap++;
		}

		for (int i = 0; i < rows; i++)
		{
			int targetRow = (i - (gap * sub_rows) + rows) % rows;	// find target row

			colBuffer[targetRow] = MB[i * cols + j];
		}

		for (int i = 0; i < rows; i++)
		{
			MB[i * cols + j] = colBuffer[i];
		}
	}

	delete[] colBuffer;
}

void initialSendMatrices(const double * MA, const double * MB, int rows, int cols, int size)
{
	int subMatrixElems = (rows * cols) / size;	// number of elements in submatrix
	int subMatrixDim = (int) sqrt(subMatrixElems);		// dimension of submatrix
	
	double *bufferMA = new double[subMatrixElems];
	double *bufferMB = new double[subMatrixElems];
	int bufferPos = 0;

	int dimensions[2] = { subMatrixDim, subMatrixDim };

	int procWidth = (int) sqrt(size);	// dimension of grid

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

void receiveDimensions(int &rows, int &cols)
{
	int * dimensions = new int[2];

	MPI_Recv(dimensions, 2, MPI_INT, 0, TAG_MATRIX_DIMENSIONS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	rows = dimensions[0];
	cols = dimensions[1];

	delete[] dimensions;
}

void initialReceiveMatrices(double * MA, double * MB, int subMatrixElems, int rank, int size)
{
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

void moveMatrix(const double *M, int elems, int direction, int tag, int rank, int size)
{
	int procWidth = (int) sqrt(size);

	int location[2] = {
		rank / procWidth
		,rank % procWidth
	};

	int targetX = location[0];
	int targetY = location[1];

	if (direction == DIRECTION_LEFT)
	{
		targetY = (targetY - 1 + procWidth) % procWidth;
	}
	else if (direction == DIRECTION_UP)
	{
		targetX = (targetX - 1 + procWidth) % procWidth;
	}

	int target = targetX * procWidth + targetY;

	MPI_Send(M, elems, MPI_DOUBLE, target, tag, MPI_COMM_WORLD);
}

void receiveMatrices(double *MA, double *MB, int elems, int rank, int size)
{
	double *buffer = new double[elems];

	MPI_Status status;

	/*for (int i = 0; i < 2; i++)
	{
		MPI_Recv(buffer, elems, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		if (status.MPI_TAG == TAG_MATRIX_A)
		{
			memcpy(MA, buffer, sizeof(*buffer) * elems);
		}
		else memcpy(MB, buffer, sizeof(*buffer) * elems);
	}*/

	int procWidth = (int) sqrt(size);

	int location[2] = {
		rank / procWidth
		, rank % procWidth
	};

	int bottomX = (location[0] + 1 + procWidth) % procWidth;
	int bottomY = location[1];

	int rightX = location[0];
	int rightY = (location[1] + 1 + procWidth) % procWidth;

	int bottom = (bottomX * procWidth) + bottomY;
	int right = (rightX * procWidth) + rightY;

	MPI_Recv(MA, elems, MPI_DOUBLE, right, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	MPI_Recv(MB, elems, MPI_DOUBLE, bottom, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

	delete[] buffer;
}

void sendResult(const double *M, int elems, int target)
{
	MPI_Send(M, elems, MPI_DOUBLE, target, TAG_MATRIX_RESULT, MPI_COMM_WORLD);
}

void collectResult(double *M, int rows, int cols, int size)
{
	int procWidth = (int) sqrt(size);

	int sub_rows = (int) sqrt((rows * cols) / size);
	int sub_cols = sub_rows;

	int elems = sub_rows * sub_cols;

	double *buffer = new double[elems];
	int bufferPos = 0;

	MPI_Status status;

	MPI_Recv(buffer, elems, MPI_DOUBLE, MPI_ANY_SOURCE, TAG_MATRIX_RESULT, MPI_COMM_WORLD, &status);

	int rank = status.MPI_SOURCE;

	int location[2] = {
		rank / procWidth
		,rank % procWidth
	};

	int xStart = location[0] * sub_rows;
	int xEnd = xStart + sub_rows;

	int yStart = location[1] * sub_cols;
	int yEnd = yStart + sub_cols;

	for (int x = xStart; x < xEnd; x++)
	{
		for (int y = yStart; y < yEnd; y++)
		{
			M[x * cols + y] = buffer[bufferPos++];
		}
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