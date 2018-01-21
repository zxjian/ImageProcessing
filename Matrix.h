#pragma once
#include <iostream>
#include <cmath>
using namespace std;
class Matrix
{
private:
	int row;
	int column;
	double* matrix;
public:
	Matrix();
	Matrix(int _row, int _column);
	Matrix Correlation(const Matrix &kernel);
	Matrix Convolution(const Matrix &kernel);
	Matrix Derivative(const Matrix &derivative_mask, int mask_length);
	int Row();
	int Column();
	void SetValue(double value, int i, int j);
	double GetValue(int i, int j);
	void PrintMatrix();
};