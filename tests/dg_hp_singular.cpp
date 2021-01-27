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

double exact(double x)
{
	if (x <= double(1)/3)
		return 0;
	else
		return pow(x - double(1)/3, 2.5)*(1-x);
}

double exact_(double x)
{
	if (x <= double(1)/3)
		return 0;
	else
		return pow(x - double(1)/3, 1.5)*(17-21*x)/(6*sqrt(3));
}

double exact__(double x)
{
	if (x <= double(1)/3)
		return 0;
	else
		return pow(x - double(1)/3, 0.5)*(1-x)*double(15)/4 - 5*pow(x-double(1)/3, 1.5);
}

double f(double x)
{
	return -exact__(x);
}

int main()
{
	// Sets up problem.
	Mesh*               myMesh     = new Mesh(2);
	Solution_dg_linear* mySolution = new Solution_dg_linear(myMesh, f, 1, zero);

	// Refinement variables.
	Mesh*               myNewMesh;
	Solution_dg_linear* myNewSolution_type;
	Solution*           myNewSolution = myNewSolution_type;

	refinement::refinement(myMesh, &myNewMesh, mySolution, &myNewSolution, 1e-15, 0, 10, true, true, true, exact, exact_);

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