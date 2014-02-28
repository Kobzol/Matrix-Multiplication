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

void setMatrices(double *& matrixA, double *& matrixB, double *& matrixC, int &rows, int &cols);

void shuffleMatrices(double *matrixA, double *matrixB, int rows, int cols, int size);
void initialSendMatrices(const double * matrixA, const double * matrixB, int rows, int cols, int size);
void receiveAndSetDimensions(double * &matrixA, double * &matrixB, double * &matrixC, int &rows, int &cols, int size);
void initialReceiveMatrices(double *matrixA, double *matrixB, int rows, int cols, int rank, int size);

void resendMatrices();
void multiplyMatrices(const double *matrixA, const double *matrixB, double *matrixC, int rows, int cols);

void returnMatrix();

void composeMatrix();

void printMatrix(double * matrix, int rows, int cols);

#endif