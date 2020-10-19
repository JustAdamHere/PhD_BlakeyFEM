#include "./src/linearSystems.hpp"
#include "./src/matrix.hpp"
#include "./src/matrix_full.hpp"
#include <iostream>
#include <vector>

int main()
{
	Matrix_full<double> M(3, 3, 0);
	M.set(0, 0, 1);
	M.set(1, 1, 1);
	M.set(2, 2, 1);

	std::vector<double> v(3, 1);

	std::vector<double> x = linearSystems::GaussJordan(M, v);

	std::cout << "x[0] = " << x[0] <<
	std::endl << "x[1] = " << x[1] <<
	std::endl << "x[2] = " << x[2];

	return 0;
}