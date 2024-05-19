#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <chrono>
#include "Matrix.h"
#include "Sample.h"
using namespace std;

// Sta³e
const int LAGRANGE = 0;
const int SPLINES = 1;

// Liczba próbek z których zostan¹ wybrane wêz³y interpolacyjne
const int sampleNumber = 512;

// Nazwy folderów z plikami
const string outputFolder = "Output";
const string inputFolder = "Input";
const string methodLName = "Lagrange";
const string methodFSName = "Funkcje sklejane";

// Nazwy plików z danymi wejœciowymi i ich liczba
const string inputSamples[] = { "genoa_rapallo.txt", "MountEverest.csv", "SpacerniakGdansk.csv", "WielkiKanionKolorado.csv" };
const int inputNumber = sizeof(inputSamples) / sizeof(string);

// Liczba próbek które zostan¹ pominiête podczas wybierania wêz³ów interpolacyjnych
const int skippedSamples[] = { 16,32,64,128 };
const int skippedSampleNumber = sizeof(skippedSamples) / sizeof(int);

// Pivoting
void pivotFunc(Matrix* U, Matrix* L, Matrix* P, int i)
{
	double pivot = abs(U->tab[i][i]);
	int pivotIndex = i;

	for (int j = i + 1; j < U->size; j++)
	{
		if (abs(U->tab[j][i]) > pivot)
		{
			pivot = abs(U->tab[j][i]);
			pivotIndex = j;
		}
	}

	// Macierz osobliwa
	if (!U->tab[pivotIndex][i])
		return;

	if (pivotIndex != i)
	{
		for (int j = 0; j < U->size; j++)
		{
			if (j >= i)
				swap(U->tab[i][j], U->tab[pivotIndex][j]);
			else
				swap(L->tab[i][j], L->tab[pivotIndex][j]);

			swap(P->tab[i][j], P->tab[pivotIndex][j]);
		}
	}
}

// Metoda faktoryzacji LU
void LU(const int size, Matrix* A, std::vector<double>& x, std::vector<double>& b)
{
	// Rozbicie macierzy
	Matrix* U = new Matrix(*A);
	Matrix* L = new Matrix(size);
	Matrix* P = new Matrix(size);

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			if (i == j)
			{
				L->tab[i][j] = 1;
				P->tab[i][j] = 1;
			}

	// A = L * U
	for (int i = 0; i < size - 1; i++)
	{
		pivotFunc(U, L, P, i);

		for (int j = i + 1; j < size; j++)
		{
			L->tab[j][i] = U->tab[j][i] / U->tab[i][i];

			for (int k = i; k < size; k++)
			{
				U->tab[j][k] = U->tab[j][k] - L->tab[j][i] * U->tab[i][k];
			}
		}
	}

	// L * y = b
	b = *P * b;
	vector<double> y(size, 0);

	for (int i = 0; i < size; i++)
	{
		double sum = 0;
		for (int j = 0; j < i; j++)
			sum += L->tab[i][j] * y[j];
		y[i] = (b[i] - sum) / L->tab[i][i];
	}

	// U * x = y
	for (int i = size - 1; i >= 0; i--)
	{
		double sum = 0;
		for (int j = i + 1; j < size; j++)
			sum += U->tab[i][j] * x[j];
		x[i] = (y[i] - sum) / U->tab[i][i];
	}

	delete U;
	delete L;
	delete P;
}

// Metoda Lagrange'a
double Lagrange(Sample* samples, int nodeNumber, double distance)
{
	double elevation = 0.00f;

	for (int i = 0; i < nodeNumber; i++)
	{
		double a = 1.00f;

		for (int j = 0; j < nodeNumber; j++)
			if (i != j)
				a *= (distance - samples[j].distance) / (samples[i].distance - samples[j].distance);

		elevation += a * samples[i].elevation;
	}

	return elevation;
}

// Metoda funkcji sklejanych
double Spline(Sample* samples, int nodeNumber, double distance)
{
	int N = 4 * (nodeNumber - 1);
	Matrix* M = new Matrix(N);
	vector<double> x(N, 1);
	vector<double> b(N, 0);

	M->tab[0][0] = 1;
	b[0] = samples[0].elevation;

	double x_dist = samples[1].distance - samples[0].distance;
	M->tab[1][0] = 1;
	M->tab[1][1] = x_dist;
	M->tab[1][2] = x_dist * x_dist;
	M->tab[1][3] = x_dist * x_dist * x_dist;
	b[1] = samples[1].elevation;

	M->tab[2][2] = 1;
	b[2] = 0;

	x_dist = samples[nodeNumber - 1].distance - samples[nodeNumber - 2].distance;
	M->tab[3][4 * (nodeNumber - 2) + 2] = 2;
	M->tab[3][4 * (nodeNumber - 2) + 3] = 6 * x_dist;
	b[3] = 0;

	for (int i = 1; i < nodeNumber - 1; i++)
	{
		x_dist = samples[i].distance - samples[i - 1].distance;

		M->tab[4 * i][4 * i] = 1;
		b[4 * i] = samples[i].elevation;

		M->tab[4 * i + 1][4 * i] = 1;
		M->tab[4 * i + 1][4 * i + 1] = x_dist;
		M->tab[4 * i + 1][4 * i + 2] = x_dist * x_dist;
		M->tab[4 * i + 1][4 * i + 3] = x_dist * x_dist * x_dist;
		b[4 * i + 1] = samples[i + 1].elevation;

		M->tab[4 * i + 2][4 * (i - 1) + 1] = 1;
		M->tab[4 * i + 2][4 * (i - 1) + 2] = 2 * x_dist;
		M->tab[4 * i + 2][4 * (i - 1) + 3] = 3 * x_dist * x_dist;
		M->tab[4 * i + 2][4 * i + 1] = -1;
		b[4 * i + 2] = 0;

		M->tab[4 * i + 3][4 * (i - 1) + 2] = 2;
		M->tab[4 * i + 3][4 * (i - 1) + 3] = 6 * x_dist;
		M->tab[4 * i + 3][4 * i + 2] = -2;
		b[4 * i + 3] = 0;
	}

	LU(N, M, x, b);
	double elevation = 0.00f;

	for (int i = 0; i < nodeNumber - 1; i++)
	{
		elevation = 0.00f;

		if (distance >= samples[i].distance && distance <= samples[i + 1].distance)
		{
			for (int j = 0; j < 4; j++)
			{
				double temp_dist = distance - samples[i].distance;
				elevation += x[4 * i + j] * pow(temp_dist, j);
			}

			break;
		}
	}

	delete M;

	return elevation;
}

// G³ówna funkcja interpoluj¹ca
int interpolation(Sample* samples, Sample* nodes, int nodeNumber, int method, string filename)
{
	// Plik z wynikami
	ofstream output;
	string methodName = (method == LAGRANGE) ? methodLName : methodFSName;
	string directory = outputFolder + '/' + methodName + '/' + filename;
	string nodeName = '_' + to_string(nodeNumber) + '.';
	directory.replace(directory.find('.'), 1, nodeName);
	output.open(directory.c_str());

	if (!output.good())
	{
		cout << "Nie mozna otworzyc pliku " << directory << endl;
		return 0;
	}

	double result = 0.00f;
	chrono::steady_clock::time_point t_start = chrono::steady_clock::now();

	double tempDist = nodes[0].distance;
	while (tempDist <= nodes[nodeNumber - 1].distance)
	{
		bool interp = true;
		for (int i = 0; i < nodeNumber; i++)
		{
			if ((int)nodes[i].distance == tempDist)
			{
				// Nie trzeba przeprowadzaæ interpolacji
				interp = false;
				result = nodes[i].elevation;
				break;
			}
		}

		// Przeprowadzenie interpolacji
		if (interp)
			result = (method == LAGRANGE) ? Lagrange(nodes, nodeNumber, tempDist) : Spline(nodes, nodeNumber, tempDist);

		// Zapisanie wyniku
		output << tempDist << ' ' << result << endl;
		tempDist += 8;
	}

	output.close();
	chrono::steady_clock::time_point t_end = chrono::steady_clock::now();
	int time = (int)chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count();

	// Zapisanie wêz³ów interpolacyjnych (abc.txt -> abc_node.txt)
	directory.replace(directory.find('.'), 1, "_nodes.");
	output.open(directory.c_str());

	if (!output.good())
	{
		cout << "Nie mozna otworzyc pliku " << directory << endl;
		return time;
	}

	for (int i = 0; i < nodeNumber; i++)
		output << nodes[i].distance << ' ' << nodes[i].elevation << endl;

	output.close();
	return time;
}

// Wypisanie œredniego czasu wykonania w zale¿noœci od metody
void printAverageTime(int* tab)
{
	for (int i = 0; i < skippedSampleNumber; i++)
	{
		int nodeNumber = sampleNumber / skippedSamples[i];
		int elapsed = tab[i] / inputNumber;
		cout << nodeNumber << " wezlow: " << elapsed << " ms" << endl;
	}
}

int main()
{
	// Wektory przechowuj¹ce œrednie czasy wykonania
	int* averageTimeL = new int[skippedSampleNumber] {};
	int* averageTimeFS = new int[skippedSampleNumber] {};

	for (int i = 0; i < inputNumber; i++)
	{
		// Dane wejœciowe
		ifstream input;
		string directory = inputFolder + '/' + inputSamples[i];
		input.open(directory.c_str());

		if (!input.good())
		{
			cout << "Nie mozna otworzyc pliku " << directory << endl;
			return 1;
		}

		// Wczytanie danych
		Sample* samples = new Sample[sampleNumber];
		int load_index = 0;

		while (!input.eof() && load_index < sampleNumber)
		{
			input >> samples[load_index].distance;
			input >> samples[load_index].elevation;
			load_index++;
		}

		// Wczytano dane
		input.close();

		if (load_index < sampleNumber)
		{
			cout << "Za mala liczba probek w pliku " << directory << endl;
			return 1;
		}

		// Ró¿na liczba wêz³ów interpolacji
		for (int j = 0; j < skippedSampleNumber; j++)
		{
			int nodeNumber = sampleNumber / skippedSamples[j];
			Sample* nodes = new Sample[nodeNumber];

			// Przepisanie wêz³ów
			for (int k = 0, l = 0; l < nodeNumber; k += skippedSamples[j], l++)
			{
				nodes[l].distance = samples[k].distance;
				nodes[l].elevation = samples[k].elevation;
			}

			// Funkcja interpolacyjna
			averageTimeL[j] += interpolation(samples, nodes, nodeNumber, LAGRANGE, inputSamples[i]);
			averageTimeFS[j] += interpolation(samples, nodes, nodeNumber, SPLINES, inputSamples[i]);

			delete[] nodes;
		}

		// Dealokacja
		delete[] samples;
	}

	// Czas metod¹ Lagrange'a
	cout << "Sredni czas wykonania - metoda Lagrange'a" << endl;
	printAverageTime(averageTimeL);

	// Czas metod¹ funkcji sklejanych
	cout << "Sredni czas wykonania - metoda funkcji sklejanych" << endl;
	printAverageTime(averageTimeFS);

	delete[] averageTimeL;
	delete[] averageTimeFS;

	return 0;
}