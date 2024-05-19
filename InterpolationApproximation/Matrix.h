#pragma once
#include <vector>

class Matrix
{
public:
	// Konstruktory
	Matrix(int N);
	Matrix(const Matrix& A);

	// Destruktor
	~Matrix();

	// Operator mno¿enia (poprawiony: memory leak fix)
	std::vector<double> operator*(const std::vector<double>& v);

	//private:
	double** tab = nullptr;
	int size = 0;
};