#include "../src/common.hpp"
#include "../src/element.hpp"
#include "../src/matrix.hpp"
#include "../src/mesh.hpp"
#include "../src/refinement.hpp"
#include "../src/solution.hpp"
#include "../src/solution_linear.hpp"
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
	Mesh*     myMesh     = new Mesh(4);
	Solution_linear* mySolution = new Solution_linear(myMesh, f, 1, zero);

	// Adaptivity variables.
	Mesh*            myNewMesh;
	Solution_linear* myNewSolution_linear;
	Solution*        myNewSolution = myNewSolution_linear;

	// Performs the refinement with the correct type of adaptivity.
	refinement::refinement(myMesh, &myNewMesh, mySolution, &myNewSolution, 1e-15, 1e-3, 3, true, false, true, exact, exact_);

	// Solves the new problem, and then outputs solution and mesh to files.
	myNewSolution->output_solution(exact);
	myNewSolution->Solve(1e-15);
	myNewSolution->output_mesh();

	delete mySolution;
	delete myMesh;

	return 0;
}