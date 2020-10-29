#include "../src/common.hpp"
#include "../src/element.hpp"
#include "../src/matrix.hpp"
#include "../src/mesh.hpp"
#include "../src/refinement.hpp"
#include "../src/solution.hpp"
#include "../src/solution_linear.hpp"
#include "../src/solution_dg_linear.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <functional>

double zero(double x)
{
	return 0;
}

double one(double x)
{
	return 1;
}

// Should yield u(x) = x(1-x)e^x
double test(double x)
{
	return x*(x+3)*exp(x);
}

double exact(double x)
{
	double a = 1e-3;

	return -exp(x/sqrt(a))/(exp(double(1)/sqrt(a)) + 1) - (exp(-x/sqrt(a)) * exp(double(1)/sqrt(a)))/(exp(double(1)/sqrt(a)) + 1) + 1;
}

double exact_(double x)
{
	double a = 1e-3;

	return -exp(x/sqrt(a))/(exp(double(1)/sqrt(a)) + 1)/sqrt(a) + (exp(-x/sqrt(a)) * exp(double(1)/sqrt(a)))/(exp(double(1)/sqrt(a)) + 1)/sqrt(a);
}

int main()
{
	// Sets up problem.
	Mesh*     myMesh     = new Mesh(4);
	Solution_dg_linear* mySolution = new Solution_dg_linear(myMesh, test);	

	// Solves the new problem, and then outputs solution and mesh to files.
	mySolution->Solve(1e-10);
	mySolution->output_solution(exact);

	delete mySolution;
	delete myMesh;

	return 0;
}