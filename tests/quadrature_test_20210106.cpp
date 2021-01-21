#include "../src/common.hpp"
#include "../src/quadrature.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <functional>

int main()
{
	int    n = 9;
	double h = double(2)/(n-1);

	for (int i=0; i<n; ++i)
	{
		std::cout << quadrature::legendrePolynomial(5, 3)(-1 + i*h) << std::endl;
	}

	return 0;
}