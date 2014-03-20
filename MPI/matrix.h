#ifndef MATRIX_H_
#define MATRIX_H_

#include <algorithm>
#include <iostream>
#include <cstring>
#include <cmath>

#include "mpi.h"

const int TAG_MATRIX_A			= 0x00000001;
const int TAG_MATRIX_B			= 0x00000002;
const int TAG_MATRIX_DIMENSIONS = 0x00000003;
const int TAG_MATRIX_RESULT		= 0x00000004;

const int DIRECTION_LEFT		= 0x00000001;
const int DIRECTION_UP			= 0x00000002;

/* Receive (or set) initial matrices for the computation */
void initializeMatrices(double * &MA, double * &MB, double * &MC, int &rows, int &cols);

/* Allocates three matrices and sets the matrix C to zero */
void setMatrices(double *& matrixA, double *& matrixB, double *& matrixC, int rows, int cols);

/* Initial shuffling of submatrices */
void shuffleMatrices(double *matrixA, double *matrixB, int rows, int cols, int size);

/* Initial emission of submatrices to slaves */
void initialSendMatrices(const double * matrixA, const double * matrixB, int rows, int cols, int size);

/* Receive submatrix dimensions from master */
void receiveDimensions(int &rows, int &cols);

/* Receive initial submatrix from master */
void initialReceiveMatrices(double *matrixA, double *matrixB, int subMatrixElems, int rank, int size);

/* Multiply two matrices */
void multiplyMatrices(const double *matrixA, const double *matrixB, double *matrixC, int rows, int cols);

/* Send submatrix to a neighbour process, indicated by direction */
void moveMatrix(double *M, int elems, int direction, int tag, int rank, int size);

/* Receive submatrices from neighbour processes */
void receiveMatrices(double *MA, double *MB, int elems, int rank, int size);

/* Send complete submatrix of matrix C to master */
void sendResult(double *M, int elems, int target);

/* Receive complete submatrix of matrix C from slave */
void collectResult(double *M, int rows, int cols, int size);

/* [Debug] Print a matrix */
void printMatrix(double * matrix, int rows, int cols);

#endif