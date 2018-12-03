#define _CRT_SECURE_NO_WARNINGS
#include <mpi.h>
#include <thread>
#include <iostream>
#include <ctime>
#include <time.h>
#include "parser.h"
#include "C:\Users\alsol\source\repos\MPI_Trapezoid_Method\MPI_Trapezoid_Method\parser.cpp"
using namespace std;


double val_func(double x1, double x2, char* formula, TParser *parser, double *x) {
	//TParser parser;
	//static double x[10];
	//parser.SetX((double *)&x);
	x[0] = x1;
	x[1] = x2;
	//parser->Compile(formula);
	parser->Evaluate();
	return parser->GetResult();
	//return 1;
}

void two_dimensional_integral(const double a1, const double b1, const double a2,
	const double b2, const double h, double *res, char* func, TParser *parser, double *x)
{	
	double sum;
	double midval;
	sum = 0.0;

	for (double i = a1; i < b1; i += h)
	{
		for (double j = a2; j < b2; j += h)
		{
			midval = (val_func(i, j, func, parser, x) + val_func(i + h, j, func, parser, x) + val_func(i, j + h, func, parser, x) + val_func(i + h, j + h, func, parser, x)) / 4;
			sum += midval * h * h;

		}
	}
	*res = sum;
}


void one_dimensional_integral(const double a1, const double b1, const double h, double *res, char* func, TParser *parser, double *x) {
	double sum = 0.0;
	double h_left, h_right;
	for (double i = a1; i < b1; i += h) {
		h_left = val_func(i, 0, func, parser, x);
		h_right = val_func(i + h, 0, func, parser, x);
		sum += ((h_right + h_left) / 2) * h;
	}
	*res = sum;
}


int main(int argc, char **argv)
{
	double a1, b1, a2 = INFINITY, b2 = INFINITY, h;
	a1 = 0.0;
	b1 = 4.0;
	a2 = 0.0;
	b2 = 4.0;
	h = 0.001;
	const char* func = "x[0] * x[1]";
	double stepX, leftX, rightX, sum = 0.0;
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
		if ((!strcmp(argv[i], "-a2")) && (i + 1 < argc)) {
			a2 = atof(argv[i + 1]);
		}
		if ((!strcmp(argv[i], "-b2")) && (i + 1 < argc)) {
			b2 = atof(argv[i + 1]);
		}
	}

	TParser parser;
	static double x[10];
	parser.SetX((double *)&x);
	parser.Compile(const_cast<char *>(func));

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &process_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		starttime = MPI_Wtime();
	}
	stepX = (b1 - a1) / process_num;  // Each process counts its right and left borders
	leftX = rank * stepX;
	rightX = leftX + stepX;

	std::cout <<"hello. i'm "<< rank << " rank. " << "leftX: " << leftX << " rightX: " << rightX <<'\n';
	if (a2 == INFINITY || b2 == INFINITY) {
		one_dimensional_integral(leftX, rightX, h, &res, const_cast<char *>(func), &parser, x); // Ñounting the value of the one dimensional integral
	}
	else {
		two_dimensional_integral(leftX, rightX, a2, b2, h, &res, const_cast<char *>(func), &parser, x); // Ñounting the value of the two dimensional integral
	}

	std::cout << "res: " << res << " \n";

	MPI_Reduce(&res, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		endtime = MPI_Wtime();
		std::cout << "Integral value of the function " << func << " = " << sum << "\n";
		std::cout << "Parallel time for " << process_num << ": " << endtime - starttime << "\n";
	}


	MPI_Finalize(); //çàêàí÷èâåì ðàáîòó MPI
	return 0;
}