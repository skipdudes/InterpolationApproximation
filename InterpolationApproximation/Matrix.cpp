#include "Matrix.h"

Matrix::Matrix(int N)
{
	tab = new double* [N];
	size = N;

	for (int i = 0; i < N; i++)
		tab[i] = new double[N] {};
}

Matrix::Matrix(const Matrix& A)
{
	tab = new double* [A.size];
	size = A.size;

	for (int i = 0; i < A.size; i++)
		tab[i] = new double[A.size] {};

	// Przepisanie wartoœci
	for (int i = 0; i < A.size; i++)
		for (int j = 0; j < A.size; j++)
			tab[i][j] = A.tab[i][j];
}

Matrix::~Matrix()
{
	for (int i = 0; i < size; i++)
		delete[] tab[i];
	delete[] tab;
}

std::vector<double> Matrix::operator*(const std::vector<double>& v)
{
	std::vector<double> u(size, 0);

	for (int i = 0; i < size; i++)
	{
		double sum = 0;
		for (int j = 0; j < size; j++)
			sum += tab[i][j] * v[j];
		u[i] = sum;
	}

	return u;
}