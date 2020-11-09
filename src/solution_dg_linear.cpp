/******************************************************************************
 * @details This is a file containing definitions of [Solution].
 * 
 * @author     Adam Matthew Blakey
 * @date       2020/10/09
 ******************************************************************************/
#include "common.hpp"
#include "element.hpp"
#include "linearSystems.hpp"
#include "matrix.hpp"
#include "matrix_full.hpp"
#include "mesh.hpp"
#include "quadrature.hpp"
#include "solution.hpp"
#include "solution_dg_linear.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

/******************************************************************************
 * __Solution_dg_linear__
 * 
 * @details 	The Mesh constructor, taking 1 argument for the mesh.
 * 
 * @param[in] a_mesh 		The mesh the solution is defined on.
 ******************************************************************************/
Solution_dg_linear::Solution_dg_linear(Mesh* const &a_mesh, f_double const &a_f)
{
	this->noElements		= a_mesh->get_noElements();
	this->solution 			.resize(a_mesh->get_noNodes());
	this->boundaryConditions.resize(2);
	this->mesh 				= a_mesh;
	this->f       = a_f; 
	this->linear  = true;
}

Solution_dg_linear::Solution_dg_linear(Mesh* const &a_mesh, Solution_dg_linear* const &a_solution)
{
	this->noElements		= a_mesh->get_noElements();
	this->solution 			.resize(a_mesh->get_noNodes());
	this->boundaryConditions.resize(2);
	this->mesh 				= a_mesh;
	this->f       = a_solution->get_f();
	this->linear  = true;
}

double Solution_dg_linear::l(Element* currentElement, f_double &basis)
{
	double J = currentElement->get_Jacobian();
	double integral = 0;

	std::vector<double> coordinates;
	std::vector<double> weights;
	currentElement->get_elementQuadrature(coordinates, weights);

	for (int i=0; i<coordinates.size(); ++i)
	{
		double f_value = this->f(currentElement->mapLocalToGlobal(coordinates[i]));
		double b_value = basis(coordinates[i]);
		
		integral += b_value*f_value*weights[i]*J;
	}

	return integral;
}

double Solution_dg_linear::a(Element* currentElement, f_double &basis1_, f_double &basis2_)
{
	double J = currentElement->get_Jacobian();
	double integral = 0;

	std::vector<double> coordinates;
	std::vector<double> weights;
	currentElement->get_elementQuadrature(coordinates, weights);

	for (int i=0; i<coordinates.size(); ++i)
	{
		double b_value = basis1_(coordinates[i]) * basis2_(coordinates[i]);

		integral += b_value*weights[i]/J;
	}

	return integral;
}

/******************************************************************************
 * __Solve__
 * 
 * @details 	Uses the stored data to calculate and populate the value in
 * 					local variable 'solution'.
 ******************************************************************************/
void Solution_dg_linear::Solve(const double &a_cgTolerance)
{	
	double A = 0;
	double B = 0;

	int n = this->mesh->elements->get_dg_DoF();

	Elements* elements = this->mesh->elements;

	Matrix_full<double> stiffnessMatrix(n, n, 0);
	std::vector<double> loadVector(n, 0);

	// Element loop.
	for (int elementNo = 0; elementNo<this->noElements; ++elementNo)
	{
		Element* currentElement = (*(this->mesh->elements))[elementNo];

		std::vector<int> elementDoFs = elements->get_dg_elementDoFs(elementNo);
		for (int a=0; a<elementDoFs.size(); ++a)
		{
			int j = elementDoFs[a];
			f_double basis = currentElement->basisLegendre(a, 0);

			loadVector[j] += this->l(currentElement, basis);

			for (int b=0; b<elementDoFs.size(); ++b)
			{
				int i = elementDoFs[b];
				f_double basis1  = currentElement->basisLegendre(b, 0);
				f_double basis2  = currentElement->basisLegendre(a, 0);
				f_double basis1_ = currentElement->basisLegendre(b, 1);
				f_double basis2_ = currentElement->basisLegendre(a, 1);

				double value = stiffnessMatrix(i, j);
				stiffnessMatrix.set(i, j, value + this->a(currentElement, basis1_, basis2_));
			}
		}
	}

	// Flux parameters.
	double theta = -1;
	double J = (*(this->mesh->elements))[0]->get_Jacobian();
	double h = 2*J;
	double sigma = 10/h;

	// Loop over faces.
	for (int faceNo=0; faceNo<this->mesh->get_noNodes(); ++faceNo)
	{
		double currentFace = this->mesh->elements->get_nodeCoordinates()[faceNo];

		// Left boundary face.
		if (faceNo == 0)
		{
			int elementNo = faceNo;
			Element* currentElement = (*(this->mesh->elements))[elementNo];
			std::vector<int> elementDoFs = elements->get_dg_elementDoFs(elementNo);

			for (int a=0; a<elementDoFs.size(); ++a)
			{
				int j = elementDoFs[a];

				// DOES SOMETHING HAVE TO BE ADDED TO THE RHS HERE?
				// Yes, it's the boundary conditions, which are zero in this test problem.

				for (int b=0; b<elementDoFs.size(); ++b)
				{
					int i = elementDoFs[b];

					double u  = currentElement->basisLegendre(b, 0)(-1);
					double v  = currentElement->basisLegendre(a, 0)(-1);
					double u_ = currentElement->basisLegendre(b, 1)(-1)/J;
					double v_ = currentElement->basisLegendre(a, 1)(-1)/J;

					double b_value = -(u_)*(-v) + theta*(v_)*(-u) + sigma*(-u)*(-v);

					double value = stiffnessMatrix(i, j);
					stiffnessMatrix.set(i, j, value + b_value);
				}
			}	
		}
		// Right boundary face.
		else if (faceNo == this->mesh->get_noNodes()-1)
		{
			int elementNo = faceNo-1;
			Element* currentElement = (*(this->mesh->elements))[elementNo];
			std::vector<int> elementDoFs = elements->get_dg_elementDoFs(elementNo);

			for (int a=0; a<elementDoFs.size(); ++a)
			{
				int j = elementDoFs[a];

				// DOES SOMETHING HAVE TO BE ADDED TO THE RHS HERE?
				// Yes, it's the boundary conditions, which are zero in this test problem.

				for (int b=0; b<elementDoFs.size(); ++b)
				{
					int i = elementDoFs[b];

					double u  = currentElement->basisLegendre(b, 0)(1);
					double v  = currentElement->basisLegendre(a, 0)(1);
					double u_ = currentElement->basisLegendre(b, 1)(1)/J;
					double v_ = currentElement->basisLegendre(a, 1)(1)/J;

					double b_value = -(u_)*(v) + theta*(v_)*(u) + sigma*(u)*(v);

					double value = stiffnessMatrix(i, j);
					stiffnessMatrix.set(i, j, value + b_value);
				}
			}	
		}
		// Interior faces.
		else
		{
			int prevElementNo = faceNo-1;
			int nextElementNo = faceNo;
			Element* prevElement = (*(this->mesh->elements))[prevElementNo];
			Element* nextElement = (*(this->mesh->elements))[nextElementNo];
			std::vector<int> prevElementDoFs = elements->get_dg_elementDoFs(prevElementNo);
			std::vector<int> nextElementDoFs = elements->get_dg_elementDoFs(nextElementNo);

			for (int a=0; a<prevElementDoFs.size(); ++a)
			{
				int j = prevElementDoFs[a];

				for (int b=0; b<prevElementDoFs.size(); ++b)
				{
					int i = prevElementDoFs[b];

					// --
					double um  = prevElement->basisLegendre(b, 0)(1);
					double vm  = prevElement->basisLegendre(a, 0)(1);
					double um_ = prevElement->basisLegendre(b, 1)(1)/J;
					double vm_ = prevElement->basisLegendre(a, 1)(1)/J;

					double b_value = -(um_)*(vm)/2 + theta*(vm_)*(um)/2 + sigma*(um)*(vm);

					double value = stiffnessMatrix(i, j);
					stiffnessMatrix.set(i, j, value + b_value);
				}

				for (int b=0; b<nextElementDoFs.size(); ++b)
				{
					int i = nextElementDoFs[b];

					// +-
					double up  = nextElement->basisLegendre(b, 0)(-1);
					double vm  = prevElement->basisLegendre(a, 0)( 1);
					double up_ = nextElement->basisLegendre(b, 1)(-1)/J;
					double vm_ = prevElement->basisLegendre(a, 1)( 1)/J;

					double b_value = -(up_)*(vm)/2 + theta*(vm_)*(-up)/2 + sigma*(-up)*(vm);

					double value = stiffnessMatrix(i, j);
					stiffnessMatrix.set(i, j, value + b_value);
				}
			}

			for (int a=0; a<nextElementDoFs.size(); ++a)
			{
				int j = nextElementDoFs[a];

				for (int b=0; b<prevElementDoFs.size(); ++b)
				{
					int i = prevElementDoFs[b];

					// -+
					double um  = prevElement->basisLegendre(b, 0)( 1);
					double vp  = nextElement->basisLegendre(a, 0)(-1);
					double um_ = prevElement->basisLegendre(b, 1)( 1)/J;
					double vp_ = nextElement->basisLegendre(a, 1)(-1)/J;

					double b_value = -(um_)*(-vp)/2 + theta*(vp_)*(um)/2 + sigma*(um)*(-vp);

					double value = stiffnessMatrix(i, j);
					stiffnessMatrix.set(i, j, value + b_value);
				}

				for (int b=0; b<nextElementDoFs.size(); ++b)
				{
					int i = nextElementDoFs[b];

					// ++
					double up  = nextElement->basisLegendre(b, 0)(-1);
					double vp  = nextElement->basisLegendre(a, 0)(-1);
					double up_ = nextElement->basisLegendre(b, 1)(-1)/J;
					double vp_ = nextElement->basisLegendre(a, 1)(-1)/J;

					double b_value = -(up_)*(-vp)/2 + theta*(vp_)*(-up)/2 + sigma*(-up)*(-vp);

					double value = stiffnessMatrix(i, j);
					stiffnessMatrix.set(i, j, value + b_value);
				}
			}
		}
	}

	this->solution = linearSystems::GaussJordan(stiffnessMatrix, loadVector);
}

f_double Solution_dg_linear::get_f() const
{
	return this->f;
}

double Solution_dg_linear::compute_residual(const double &a_uh, const double &a_uh_2, const double &a_x) const
{
	return 0;
}

// Not entirely sure that this is correct...
double Solution_dg_linear::compute_energyNormDifference2(f_double const &a_u, f_double const &a_u_1) const
{
	int n = this->mesh->get_noElements();

	double norm = 0;

	for (int i=0; i<n; ++i)
	{
		// Gets the current element.
		Element* currentElement = (*(this->mesh->elements))[i];

		// Retrieves quadrature information.
		std::vector<double> coordinates;
		std::vector<double> weights;
		currentElement->get_elementQuadrature(coordinates, weights);

		for (int j=0; j<coordinates.size(); ++j)
		{			
			// Actual and approximate solution at coordinates.
			double uh   = compute_uh(i, coordinates[j], 0);
			double uh_1 = compute_uh(i, coordinates[j], 1);
			double u    = a_u  (currentElement->mapLocalToGlobal(coordinates[j]));
			double u_1  = a_u_1(currentElement->mapLocalToGlobal(coordinates[j]));

			double Jacobian = currentElement->get_Jacobian();

			norm += pow((u_1 - uh_1), 2)*weights[j]*Jacobian
				 +  pow(sqrt(this->f(currentElement->mapLocalToGlobal(coordinates[j])))*(u - uh), 2)*weights[j]*Jacobian;
		}
	}

	return norm;
}

double Solution_dg_linear::compute_errorIndicator(const double &a_i) const
{
	return 0;
}

double Solution_dg_linear::compute_uh(const int &a_i, const double &a_xi, const int &a_n) const
{
	Element* currentElement = (*(this->mesh->elements))[a_i];
	double J = currentElement->get_Jacobian(); // Needs to be inverse transpose of Jacobi in dimensions higher than 1.

	double result = 0;

	std::vector<int> elementDoFs = this->mesh->elements->get_dg_elementDoFs(a_i);
	for (int j=0; j<elementDoFs.size(); ++j)
	{
		f_double basis = (*(this->mesh->elements))[a_i]->basisLegendre(j, a_n);

		result += this->solution[elementDoFs[j]] * basis(a_xi);
	}

	return result / pow(J, a_n);
}