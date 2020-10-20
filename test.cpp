#include "./src/linearSystems.hpp"
#include "./src/matrix.hpp"
#include "./src/matrix_full.hpp"
#include <iostream>
#include <random>
#include <ctime>
#include <vector>

int main()
{
	std::srand(std::time(NULL));

	int n = 4;

	Matrix_full<double> M(n, n, 0);
	for (int i=0; i<n; ++i)
		for (int j=0; j<n; ++j)
			M.set(i, j, double(std::rand() % 1000)/1000);

	std::vector<double> v(4, 1);

	std::vector<double> x = linearSystems::GaussJordan(M, v);

	std::cout << "x[0] = " << x[0] <<
	std::endl << "x[1] = " << x[1] <<
	std::endl << "x[2] = " << x[2] <<
	std::endl << "x[3] = " << x[3];

	return 0;
}