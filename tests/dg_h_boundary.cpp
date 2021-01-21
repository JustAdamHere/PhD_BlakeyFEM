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

double exact_(double x)
{
	return M_PI*cos(M_PI*x);
}

int main()
{
	// Sets up problem.
	Mesh*               myMesh     = new Mesh(4);
	Solution_dg_linear* mySolution = new Solution_dg_linear(myMesh, f);

	// Refinement variables.
	Mesh*               myNewMesh;
	Solution_dg_linear* myNewSolution_type;
	Solution*           myNewSolution = myNewSolution_type;

	refinement::refinement(myMesh, &myNewMesh, mySolution, &myNewSolution, 1e-15, 1e-3, 15, true, false, true, exact, exact_);

	// Solves the new problem, and then outputs solution and mesh to files.
	myNewSolution->output_solution(exact);
	myNewSolution->Solve(1e-15);
	myNewSolution->output_mesh();

	//delete myNewSolution; // DEFFO a memory leak somewhere.
	delete mySolution;
	//delete myNewMesh;
	delete myMesh;

	return 0;
}