#ifndef MATRIX_H_
#define MATRIX_H_

#include <algorithm>
#include <iostream>
#include <cstring>

void setMatrices(double *& matrixA, double *& matrixB, double *& matrixC, int &rows, int &cols);

void shuffleMatrices(double *& matrixA, double *& matrixB, int rows, int cols, int size);
void initialSendMatrices();
void initialReceiveMatrices();

void resendMatrices();
void computeProduct();

void returnMatrix();

void composeMatrix();

void printMatrix(double * matrix, int rows, int cols);

#endif