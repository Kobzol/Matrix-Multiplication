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

void initializeMatrices(double * &MA, double * &MB, double * &MC, int &rows, int &cols);

void setMatrices(double *& matrixA, double *& matrixB, double *& matrixC, int rows, int cols);

void shuffleMatrices(double *matrixA, double *matrixB, int rows, int cols, int size);
void initialSendMatrices(const double * matrixA, const double * matrixB, int rows, int cols, int size);
void receiveDimensions(int &rows, int &cols);
void initialReceiveMatrices(double *matrixA, double *matrixB, int subMatrixElems, int rank, int size);

void multiplyMatrices(const double *matrixA, const double *matrixB, double *matrixC, int rows, int cols);

void moveMatrix(const double *M, int elems, int direction, int tag, int rank, int size);
void receiveMatrices(double *MA, double *MB, int elems, int rank, int size);

void sendResult(const double *M, int elems, int target);
void collectResult(double *M, int rows, int cols, int size);

void printMatrix(double * matrix, int rows, int cols);

#endif