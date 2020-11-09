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

// Should yield u(x) = sin(pi*x)
double f(double x)
{
	return pow(M_PI, 2)*sin(M_PI*x);
}

double exact(double x)
{
	return sin(M_PI*x);
}

int main()
{
	// Sets up problem.
	Mesh*     myMesh     = new Mesh(8);
	Solution_dg_linear* mySolution = new Solution_dg_linear(myMesh, f);	

	// Solves the new problem, and then outputs solution and mesh to files.
	mySolution->Solve(1e-10);
	mySolution->output_solution(exact);

	delete mySolution;
	delete myMesh;

	return 0;
}