#include "../src/common.hpp"
#include "../src/element.hpp"
#include "../src/matrix.hpp"
#include "../src/mesh.hpp"
#include "../src/refinement.hpp"
#include "../src/solution.hpp"
#include "../src/solution_nonlinear.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <functional>

double bratu(double x, double u)
{
	return -exp(u);
}

double sinpi(double x)
{
	return sin(M_PI*x);
}

double zero(double x, double u)
{
	return 0;
}

int main()
{
	// Sets up problem.
	int n = 20;
	Mesh*               myMesh     = new Mesh(n);
	Solution_nonlinear* mySolution = new Solution_nonlinear(myMesh, bratu, bratu, 1);

	// Solves the new problem, and then outputs solution and mesh to files.
	double a = 1;
	std::vector<double> u0(n+1);
	for (int i=0; i<u0.size(); ++i)
		u0[i] = a*sinpi(double(i)/n);
	mySolution->Solve(1e-15, 1e-12, u0);
	mySolution->output_solution();
	mySolution->output_mesh();
	
	delete mySolution;
	delete myMesh;

	return 0;
}