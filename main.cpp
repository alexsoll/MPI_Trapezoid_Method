#define _CRT_SECURE_NO_WARNINGS
#include <mpi.h>
#include <thread>
#include <iostream>
#include <ctime>
#include <time.h>
#include "parser.h"
#include "C:\Users\alsol\source\repos\MPI_Trapezoid_Method\MPI_Trapezoid_Method\parser.cpp"
using namespace std;

double val_func(double x1, char* formula) {
	TParser parser;
	static double x[10];
	parser.SetX((double *)&x);
	x[0] = x1;
	parser.Compile(formula);
	parser.Evaluate();
	return parser.GetResult();
}

void integral(const double a1, const double b1, const double h, double *res, char* func) {
	double sum = 0.0;
	double h_left, h_right;
	for (double i = a1; i < b1; i += h) {
		h_left = val_func(i, func);
		h_right = val_func(i + h, func);
		sum += ((h_right + h_left) / 2) * h;
	}
	*res = sum;
}


int main(int argc, char **argv)
{
	double a1, b1, a2, b2, h;
	a1 = 0.0;
	b1 = 4.0;
	h = 0.00001;
	const char* func = "x[0] * x[0]";
	double step, left, right, sum = 0.0;
	double starttime, endtime;
	int rank, process_num;
	double res = 0.0;

	for (int i = 1; i < argc; ++i) {
		if ((!strcmp(argv[i], "-a1")) && (i + 1 < argc)) {
			a1 = atof(argv[i + 1]);
		}
		if ((!strcmp(argv[i], "-b1")) && (i + 1 < argc)) {
			b1 = atof(argv[i + 1]);
		}
		if ((!strcmp(argv[i], "-h")) && (i + 1 < argc)) {
			h = atof(argv[i + 1]);
		}
		if ((!strcmp(argv[i], "-f(x)")) && (i + 1 < argc)) {
			func = argv[i + 1];
		}
	}

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &process_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
		starttime = MPI_Wtime();

	step = (b1 - a1) / process_num;  // Each process counts its right and left borders
	left = rank * step;
	right = left + step;

	std::cout <<"hello. i'm "<< rank << " rank. " << "left: " << left << " right: " << right << '\n';
	
	integral(left, right, h, &res, const_cast<char *>(func)); // Сounting the value of the integral

	std::cout << "res: " << res << " \n";

	MPI_Reduce(&res, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		endtime = MPI_Wtime();
		std::cout << "Integral value of the function " << func << " = " << sum << "\n";
		std::cout << "Parallel time for " << process_num << ": " << endtime - starttime << "\n";
	}


	MPI_Finalize(); //заканчивем работу MPI
	return 0;
}