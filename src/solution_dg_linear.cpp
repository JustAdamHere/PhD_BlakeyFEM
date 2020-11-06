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
	// double J = currentElement->get_Jacobian();
	// double integral = 0;
	
	// std::vector<double> coordinates;
	// std::vector<double> weights;
	// currentElement->get_elementQuadrature(coordinates, weights);

	// for (int k=0; k<coordinates.size(); ++k)
	// {
	// 	double b_value = basis(coordinates[k]);

	// 	double f_value = this->f(currentElement->mapLocalToGlobal(coordinates[k]));
	// 	integral += b_value*f_value*weights[k]*J;
	// }

	// return integral;

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
	// double J = currentElement->get_Jacobian();
	// double integral = 0;
	
	// std::vector<double> coordinates;
	// std::vector<double> weights;
	// currentElement->get_elementQuadrature(coordinates, weights);

	// for (int k=0; k<coordinates.size(); ++k)
	// {
	// 	double b_value = basis1_(coordinates[k]) * basis2_(coordinates[k]);

	// 	integral += this->epsilon*b_value*weights[k]/J;
	// }

	// for (int k=0; k<coordinates.size(); ++k)
	// {
	// 	double b_value = basis1(coordinates[k]) * basis2(coordinates[k]);

	// 	double c_value = this->c(currentElement->mapLocalToGlobal(coordinates[k]));

	// 	integral += c_value*b_value*weights[k]*J;
	// }

	// return integral;

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

// See Andrea book on page 20 for maybe a better way of doing all this...
double Solution_dg_linear::b(const int &faceNo, f_double &leftBasis1, f_double &leftBasis2, f_double &leftBasis1_, f_double &leftBasis2_, f_double &rightBasis1, f_double &rightBasis2, f_double &rightBasis1_, f_double &rightBasis2_)
{
	// TEMPORARY: How do I get the correct Jacobian here?
	double J = (*(this->mesh->elements))[0]->get_Jacobian();

	double currentFace = this->mesh->elements->get_nodeCoordinates()[faceNo];
	double integral = 0;

	// Stability parameters.
	double theta = -1;
	double sigma = 10/J; // Should there be a dependance of O(p^2/h) here?

	// Variables for jumps and averages.
	double u_0_jump;
	double v_0_jump;
	double u_1_average;
	double v_1_average;

	// Left boundary.
	if (faceNo == 0)
	{
		u_0_jump    = leftBasis1 (currentFace);
		v_0_jump    = leftBasis2 (currentFace);
		u_1_average = leftBasis1_(currentFace);
		v_1_average = leftBasis2_(currentFace);
	}
	// Right boundary.
	else if (faceNo == this->mesh->get_noNodes())
	{
		u_0_jump    = -leftBasis1 (currentFace);
		v_0_jump    = -leftBasis2 (currentFace);
		u_1_average =  leftBasis1_(currentFace);
		v_1_average =  leftBasis2_(currentFace);
	}
	// Interior element boundary.
	else
	{
		// I think maybe we'll actually need 8 different basis functions?!
		u_0_jump    = jump   (leftBasis1,  rightBasis1,  currentFace); // Needs modifying.
		v_0_jump    = jump   (leftBasis2,  rightBasis2,  currentFace); // Needs modifying.
		u_1_average = average(leftBasis1_, rightBasis1_, currentFace); // Needs modifying.
		v_1_average = average(leftBasis2_, rightBasis2_, currentFace); // Needs modifying.
	}
	
	integral = - u_1_average*v_0_jump + theta*v_1_average*u_0_jump + sigma*u_0_jump*v_0_jump;

	integral /= J;

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

	// GENERAL PROBLEM
	// - (au')' + (bu)' + cu = f,    in omega
	// u = g_d,                      on omega Dirichlet boundary
	// n * (au') = g_N,              on omega Neumann boundary

	// POISSON
	// -u'' = f,    in omega
	//  u   = 0,    on omega boundary

	// int(u'v') + 

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

					if ((i == 6) && (j == 7))
						std::cout << "WOW! " << "left" << std::endl;
					/*std::cout << i << " " << j << " " << b_value << std::endl;*/
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

					if ((i == 6) && (j == 7))
						std::cout << "WOW! " << "right" << std::endl;
					/*std::cout << i << " " << j << " " << b_value << std::endl;*/
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

					if ((i == 6) && (j == 7))
						std::cout << "WOW! " << "mm" << std::endl;
					/*std::cout << i << " " << j << " " << b_value << std::endl;*/
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

					if ((i == 6) && (j == 7))
						std::cout << "WOW! " << "pm" << std::endl;
					/*std::cout << i << " " << j << " " << b_value << std::endl;*/
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

					if ((i == 6) && (j == 7))
						std::cout << "WOW! " << "mp" << std::endl;
					/*std::cout << i << " " << j << " " << b_value << std::endl;*/
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

					if ((i == 6) && (j == 7))
						std::cout << "WOW! " << "pp" << std::endl;
					/*std::cout << i << " " << j << " " << b_value << std::endl;*/
				}
			}
		}
	}

























	

	/*Matrix_full<double> ppMatrix(2, 2, 0); // Assumes linear everywhere.
	Matrix_full<double> pmMatrix(2, 2, 0);
	Matrix_full<double> mpMatrix(2, 2, 0);
	Matrix_full<double> mmMatrix(2, 2, 0);
	for (int faceNo = 0; faceNo<this->mesh->get_noNodes(); ++ faceNo)
	{
		double currentFace = this->mesh->elements->get_nodeCoordinates()[faceNo];

		if (faceNo == 0)
		{
			int elementNo = faceNo;
			Element* currentElement = (*(this->mesh->elements))[elementNo];
			std::vector<int> elementDoFs = elements->get_dg_elementDoFs(elementNo);

			for (int a=0; a<elementDoFs.size(); ++a)
			{
				for (int b=0; b<elementDoFs.size(); ++b)
				{
					int i = elementDoFs[b];
					int j = elementDoFs[a];

					double u  = currentElement->basisLegendre(a, 0)(currentFace);
					double v  = currentElement->basisLegendre(b, 0)(currentFace);
					double u_ = currentElement->basisLegendre(a, 1)(currentFace);
					double v_ = currentElement->basisLegendre(b, 1)(currentFace);

					double b_value = -(u_)*(v)/2 + theta*(v_)*(u)/2 + sigma*(u)*(v);

					double value = stiffnessMatrix(i, j);
					stiffnessMatrix.set(i, j, value + b_value);
				}
			}
		}
		else if (faceNo == this->mesh->get_noNodes()-1)
		{
			int elementNo = faceNo-1;
			Element* currentElement = (*(this->mesh->elements))[elementNo];
			std::vector<int> elementDoFs = elements->get_dg_elementDoFs(elementNo);

			for (int a=0; a<elementDoFs.size(); ++a)
			{
				for (int b=0; b<elementDoFs.size(); ++b)
				{
					int i = elementDoFs[b];
					int j = elementDoFs[a];

					double u  = currentElement->basisLegendre(a, 0)(currentFace);
					double v  = currentElement->basisLegendre(b, 0)(currentFace);
					double u_ = currentElement->basisLegendre(a, 1)(currentFace);
					double v_ = currentElement->basisLegendre(b, 1)(currentFace);

					double b_value = -(u_)*(-v)/2 + theta*(v_)*(-u)/2 + sigma*(-u)*(-v);

					double value = stiffnessMatrix(i, j);
					stiffnessMatrix.set(i, j, value + b_value);
				}
			}
		}
		else
		{
			int prevElementNo = faceNo;
			int nextElementNo = faceNo+1;
			Element* prevElement = (*(this->mesh->elements))[prevElementNo];
			Element* nextElement = (*(this->mesh->elements))[nextElementNo];
			std::vector<int> prevElementDoFs = elements->get_dg_elementDoFs(prevElementNo);
			std::vector<int> nextElementDoFs = elements->get_dg_elementDoFs(nextElementNo);

			for (int a=0; a<prevElementDoFs.size(); ++a)
			{
				for (int b=0; b<prevElementDoFs.size(); ++b)
				{
					for (int c=0; c<prevElementDoFs.size(); ++c)
					{
						for (int d=0; d<prevElementDoFs.size(); ++d)
						{
							int i = prevElementDoFs[b];
							int j = prevElementDoFs[a];
							int k = nextElementDoFs[d];
							int l = nextElementDoFs[c];

							double um  = prevElement->basisLegendre(a, 0)(currentFace);
							double vm  = prevElement->basisLegendre(b, 0)(currentFace);
							double up  = nextElement->basisLegendre(c, 0)(currentFace);
							double vp  = nextElement->basisLegendre(d, 0)(currentFace);
							double um_ = prevElement->basisLegendre(a, 1)(currentFace);
							double vm_ = prevElement->basisLegendre(b, 1)(currentFace);
							double up_ = nextElement->basisLegendre(c, 1)(currentFace);
							double vp_ = nextElement->basisLegendre(d, 1)(currentFace);

							std::cout
                                << um  << std::endl
								<< vm  << std::endl
								<< up  << std::endl
								<< vp  << std::endl
								<< um_ << std::endl
								<< vm_ << std::endl
								<< up_ << std::endl
								<< vp_ << std::endl
							<< std::endl;

							double b_value = -(up_ + um_)*(vm - vp)/2 + theta*(vp_ + vm_)*(um - up)/2 + sigma*(um - up)*(vm - vp);

							std::cout << b_value << std::endl;

							double value = stiffnessMatrix(i, j);
							stiffnessMatrix.set(i, j, value + b_value);
						}
					}
				}
			}
		}
	}*/




	// Interior face loop.
	/*for (int faceNo = 1; faceNo<this->mesh->get_noNodes()-1; ++faceNo) // Only works in 1D
	{
		double currentFace = this->mesh->elements->get_nodeCoordinates()[faceNo]; // May need modifying

		int prevElementNo = faceNo-1;
		int nextElementNo = faceNo;

		Element* prevElement = (*(this->mesh->elements))[prevElementNo];
		Element* nextElement = (*(this->mesh->elements))[nextElementNo];
		
		std::vector<int> prevElementDoFs = elements->get_dg_elementDoFs(prevElementNo);
		std::vector<int> nextElementDoFs = elements->get_dg_elementDoFs(nextElementNo);

		for (int a=0; a<nextElementDoFs.size(); ++a)
		{
			int j = nextElementDoFs[a];
			
			for (int b=0; b<prevElementDoFs.size(); ++b)
			{
				int i = prevElementDoFs[b];
				f_double leftBasis1   = prevElement->basisLegendre(b, 0);
				f_double leftBasis2   = prevElement->basisLegendre(a, 0);
				f_double leftBasis1_  = prevElement->basisLegendre(b, 1);
				f_double leftBasis2_  = prevElement->basisLegendre(a, 1);
				f_double rightBasis1  = nextElement->basisLegendre(b, 0);
				f_double rightBasis2  = nextElement->basisLegendre(a, 0);
				f_double rightBasis1_ = nextElement->basisLegendre(b, 1);
				f_double rightBasis2_ = nextElement->basisLegendre(a, 1);

				double value = stiffnessMatrix(i, j); // Bit messy...
				stiffnessMatrix.set(i, j, value + this->b(faceNo, leftBasis1, leftBasis2, leftBasis1_, leftBasis2_, rightBasis1, rightBasis2, rightBasis1_, rightBasis2_)); // Call to this will change soon...
			}
		}
	}

	// Exterior face loop.
	std::vector<int> exteriorFaces = {0, this->mesh->get_noNodes()-1};
	for (int faceCounter=0; faceCounter<exteriorFaces.size(); ++faceCounter)
	{
		int elementNo;
		if (faceCounter == 0)
			elementNo = 0;
		else if(faceCounter == 1)
			elementNo = this->mesh->get_noElements()-1;

		int faceNo = exteriorFaces[faceCounter];

		Element* currentElement = (*(this->mesh->elements))[elementNo];
		std::vector<int> elementDoFs = elements->get_dg_elementDoFs(elementNo);

		// I'M VERY UNSURE IF THIS IS CORRECT. A LITTLE UNSURE ON PLACEMENT OF EXTERIOR FACE DOFS.
		for (int a=0; a<elementDoFs.size(); ++a)
		{
			int j = elementDoFs[a];

			for (int b=0; b<elementDoFs.size(); ++b)
			{
				int i = elementDoFs[b];
				f_double basis1  = currentElement->basisLegendre(b, 0);
				f_double basis2  = currentElement->basisLegendre(a, 0);
				f_double basis1_ = currentElement->basisLegendre(b, 1);
				f_double basis2_ = currentElement->basisLegendre(a, 1);

				double value = stiffnessMatrix(i, j); // Bit messy...
				stiffnessMatrix.set(i, j, value + this->b(faceNo, basis1, basis2, basis1_, basis2_, basis1, basis1, basis1, basis1)); // Final 4 arguments are not important on exterior boundaries.
				this->b(faceNo, basis1, basis2, basis1_, basis2_, basis1, basis1, basis1, basis1);
			}
		}
	}*/

	for (int i=0; i<stiffnessMatrix.get_noColumns(); ++i)
	{
		for (int j=0; j<stiffnessMatrix.get_noRows(); ++j)
			std::cout << std::setw(8) << stiffnessMatrix(i, j);
		std::cout << std::endl;
	}

	std::cout << std::endl << std::endl;

	for (int i=0; i<loadVector.size(); ++i)
		std::cout << loadVector[i] << std::endl;



	//this->solution = linearSystems::conjugateGradient(stiffnessMatrix, loadVector, 1e-5);
	this->solution = linearSystems::GaussJordan(stiffnessMatrix, loadVector);
	//this->solution.resize(n, 0);



	// int n = this->mesh->elements->get_DoF();//this->noElements + 1; // Number of nodes.

	// Elements* elements = this->mesh->elements;

	// Matrix_full<double> stiffnessMatrix(n, n, 0);
	// std::vector<double> loadVector(n, 0);

	// for (int elementCounter=0; elementCounter<this->noElements; ++elementCounter)
	// {
	// 	Element* currentElement = (*(this->mesh->elements))[elementCounter];
	// 	int polynomialDegree = currentElement->get_polynomialDegree();

	// 	double elementLeft  = currentElement->get_nodeCoordinates()[0];
	// 	double elementRight = currentElement->get_nodeCoordinates()[1];

	// 	std::vector<int> elementDoFs = elements->get_elementDoFs(elementCounter);
	// 	for (int a=0; a<elementDoFs.size(); ++a)
	// 	{
	// 		int j = elementDoFs[a];
	// 		f_double basis = currentElement->basisFunction(a, 0);

	// 		loadVector[j] += this->l(currentElement, basis);

	// 		for (int b=0; b<elementDoFs.size(); ++b)
	// 		{
	// 			int i = elementDoFs[b];
	// 			f_double basis1  = currentElement->basisFunction(b, 0);
	// 			f_double basis2  = currentElement->basisFunction(a, 0);
	// 			f_double basis1_ = currentElement->basisFunction(b, 1);
	// 			f_double basis2_ = currentElement->basisFunction(a, 1);

	// 			double value = stiffnessMatrix(i, j); // Bit messy...
	// 			stiffnessMatrix.set(i, j, value + this->a(currentElement, basis1, basis2, basis1_, basis2_));
	// 		}
	// 	}
	// }

	// std::vector<double> F_(n);
	// std::vector<double> u0(n, 0);

	// int m = this->mesh->elements->get_noElements(); // Only works in 1D!
	
	// for (int i=0; i<stiffnessMatrix.get_noRows(); ++i)
	// 	stiffnessMatrix.set(0, i, 0);
	// for (int j=0; j<stiffnessMatrix.get_noColumns(); ++j)
	// 	stiffnessMatrix.set(j, 0, 0);
	// loadVector[0] = 0;

	// for (int i=0; i<stiffnessMatrix.get_noRows(); ++i)
	// 	stiffnessMatrix.set(m, i, 0);
	// for (int j=0; j<stiffnessMatrix.get_noColumns(); ++j)
	// 	stiffnessMatrix.set(j, m, 0);
	// loadVector[m] = 0;

	// u0[0] = A;
	// u0[m] = B;
	
	// F_ = stiffnessMatrix*u0;
	// for (int i=0; i<n; ++i)
	// 	loadVector[i] -= F_[i];

	// for (int i=0; i<stiffnessMatrix.get_noRows(); ++i)
	// 	stiffnessMatrix.set(0, i, 0);
	// for (int j=0; j<stiffnessMatrix.get_noColumns(); ++j)
	// 	stiffnessMatrix.set(j, 0, 0);
	// stiffnessMatrix.set(0, 0, 1);

	// for (int i=0; i<stiffnessMatrix.get_noRows(); ++i)
	// 	stiffnessMatrix.set(m, i, 0);
	// for (int j=0; j<stiffnessMatrix.get_noColumns(); ++j)
	// 	stiffnessMatrix.set(j, m, 0);
	// stiffnessMatrix.set(m, m, 1);

	// this->solution = linearSystems::conjugateGradient(stiffnessMatrix, loadVector, a_cgTolerance);

	// this->solution[0] = A;
	// this->solution[m] = B;
}

f_double Solution_dg_linear::get_f() const
{
	return this->f;
}

double Solution_dg_linear::compute_residual(const double &a_uh, const double &a_uh_2, const double &a_x) const
{
	// return this->f(a_x) + this->epsilon*a_uh_2 - this->c(a_x)*a_uh;

	return 0;
}

double Solution_dg_linear::compute_energyNormDifference2(f_double const &a_u, f_double const &a_u_1) const
{
	// int n = this->mesh->get_noElements();
	// double sqrt_epsilon = sqrt(this->epsilon);

	// double norm = 0;

	// for (int i=0; i<n; ++i)
	// {
	// 	// Gets the current element.
	// 	Element* currentElement = (*(this->mesh->elements))[i];

	// 	// Retrieves quadrature information.
	// 	std::vector<double> coordinates;
	// 	std::vector<double> weights;
	// 	currentElement->get_elementQuadrature(coordinates, weights);

	// 	for (int j=0; j<coordinates.size(); ++j)
	// 	{			
	// 		// Actual and approximate solution at coordinates.
	// 		double uh   = compute_uh(i, coordinates[j], 0);
	// 		double uh_1 = compute_uh(i, coordinates[j], 1);
	// 		double u    = a_u  (currentElement->mapLocalToGlobal(coordinates[j]));
	// 		double u_1  = a_u_1(currentElement->mapLocalToGlobal(coordinates[j]));

	// 		double Jacobian = currentElement->get_Jacobian();

	// 		norm += pow(sqrt_epsilon*(u_1 - uh_1), 2)*weights[j]*Jacobian
	// 			 +  pow(sqrt(this->c(currentElement->mapLocalToGlobal(coordinates[j])))*(u - uh), 2)*weights[j]*Jacobian;
	// 	}
	// }

	// return norm;

	return 0;
}

double Solution_dg_linear::compute_errorIndicator(const double &a_i) const
{
	/*
	// Gets element and its properties.
	Element* currentElement = (*(this->mesh->elements))[a_i];
	int P = currentElement->get_polynomialDegree();
	double leftNode  = currentElement->get_nodeCoordinates()[0];
	double rightNode = currentElement->get_nodeCoordinates()[1];
	double Jacobian  = currentElement->get_Jacobian();

	// Calculates L2 norm on element with weight and residual.
	double norm_2 = 0;
	std::vector<double> quadratureCoordinates;
	std::vector<double> quadratureWeights;
	currentElement->get_elementQuadrature(quadratureCoordinates, quadratureWeights);

	// Loops over quadrature coordinates and weights.
	for (int j=0; j<quadratureCoordinates.size(); ++j)
	{
		double uh   = compute_uh(a_i, quadratureCoordinates[j], 0);
		double uh_2 = compute_uh(a_i, quadratureCoordinates[j], 2);
		double residual = compute_residual(uh, uh_2, currentElement->mapLocalToGlobal(quadratureCoordinates[j]));

		double x = currentElement->mapLocalToGlobal(quadratureCoordinates[j]);
		double weight = (rightNode - x)*(x - leftNode);

		norm_2 += pow(sqrt(weight)*residual, 2)*quadratureWeights[j]*Jacobian;
	}
	
	return double(1)/(P*(P+1)*this->epsilon) * norm_2;
	*/

	return 0;
}

// WHAT ABOUT BOUNDARIES?!
//double Solution_dg_linear::average(const int &a_faceNo, const bool &a_left, f_double &basis)
double Solution_dg_linear::average(f_double &basis1, f_double &basis2, const double &a_face)
{
	return (basis1(a_face) + basis2(a_face))/2;
}

double Solution_dg_linear::jump(f_double &basis1, f_double &basis2, const double &a_face)
{
	return basis2(a_face) - basis1(a_face);
}