#pragma once
#include "Matrix.h"

Matrix::Matrix() : row(0), column(0) {}
Matrix::Matrix(int _row, int _column)
{
	row = _row;
	column = _column;
	int n = row*column;
	matrix = new double[n];
	for (int i = 0; i < n; ++i)
		matrix[i] = 0;
}

int Matrix::Row()
{
	return row;
}

int Matrix::Column()
{
	return column;
}

Matrix Matrix::Correlation(const Matrix &kernel)
{
	Matrix new_image = Matrix(row, column);
	int kernel_half_length = kernel.row / 2;
	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < column; ++j)
		{
			if (i - kernel_half_length >= 0 && i + kernel_half_length < row && j - kernel_half_length >= 0 && j + kernel_half_length < column)
			{
				double temp_value = 0;
				for (int k = 0; k < kernel.row; ++k)
				{
					for (int l = 0; l < kernel.column; ++l)
					{
						temp_value += matrix[(i+k-kernel_half_length)*column + (j+l-kernel_half_length)] * kernel.matrix[k*kernel.column+l];
					}
				}
				new_image.matrix[i*column + j] = temp_value;
			}			
		}
	}
	return new_image;
}

Matrix Matrix::Convolution(const Matrix &kernel)
{
	Matrix new_image = Matrix(row, column);
	int kernel_half_length = kernel.row / 2;
	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < column; ++j)
		{
			if (i - kernel_half_length >= 0 && i + kernel_half_length < row && j - kernel_half_length >= 0 && j + kernel_half_length < column)
			{
				double temp_value = 0;
				for (int k = kernel.row-1; k >= 0; --k)
				{
					for (int l = kernel.column-1; l >= 0; --l)
					{
						temp_value += matrix[(i - k + kernel_half_length)*column + (j - l + kernel_half_length)] * kernel.matrix[k + kernel.column + l];
					}
				}
				new_image.matrix[i*column + j] = temp_value;
			}
		}
	}
	return new_image;
}

Matrix Matrix::Derivative(const Matrix &derivate_mask, int mask_length)
{
	Matrix new_image = Matrix(row, column);
	int kernel_half_length = derivate_mask.row / 2;
	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < column; ++j)
		{
			if (i - kernel_half_length >= 0 && i + kernel_half_length < row && j - kernel_half_length >= 0 && j + kernel_half_length < column)
			{
				double temp_value = 0;
				for (int k = 0; k < derivate_mask.row; ++k)
				{
					for (int l = 0; l < derivate_mask.column; ++l)
					{
						temp_value += matrix[(i + k - kernel_half_length)*column + (j + l - kernel_half_length)] * derivate_mask.matrix[k*derivate_mask.column + l];
					}
				}
				new_image.matrix[i*column + j] = temp_value/ mask_length;
			}
		}
	}
	return new_image;
}

void Matrix::SetValue(double value, int i, int j)
{
	if (i >= 0 && i < row && j >= 0 && j < column)
	{
		matrix[i*column + j] = value;
	}
}

double Matrix::GetValue(int i, int j)
{
	if (i >= 0 && i < row && j >= 0 && j < column)
		return matrix[i*column + j];
	else
		return NULL;
}

void Matrix::PrintMatrix()
{
	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < column; ++j)
		{
			cout << matrix[i*column + j] << " ";
		}
		cout << endl;
	}
}