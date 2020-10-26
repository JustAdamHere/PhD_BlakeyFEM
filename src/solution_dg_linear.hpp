/******************************************************************************
 * @details This is a file containing declarations of the [Solution] namespace.
 * 
 * @author     Adam Matthew Blakey
 * @date       2020/10/09
 ******************************************************************************/
#ifndef CLASS_SOLUTIONDGLINEAR
#define CLASS_SOLUTIONDGLINEAR

#include "common.hpp"
#include "matrix.hpp"
#include "matrix_full.hpp"
#include "solution.hpp"
#include <vector>

class Solution_dg_linear : public Solution
{
	private:
		// Problem data.
		f_double f;

		// Computes stiffness, load vector, and extra stiffness boundary terms.
		double a(Element* currentElement, f_double &basis1_, f_double &basis2_);
		double l(Element* currentElement, f_double &basis);
		double b(const int &currentFace, f_double &leftBasis1, f_double &leftBasis2, f_double &leftBasis1_, f_double &leftBasis2_, f_double &rightBasis1, f_double &rightBasis2, f_double &rightBasis1_, f_double &rightBasis2_);

		// Computes averages and jumps required for face integrals.
		double average(f_double &basis1, f_double &basis2, const double &a_face);
		double jump   (f_double &basis1, f_double &basis2, const double &a_face);

		// Computers.
		double compute_residual(const double &a_uh, const double &a_uh_2, const double &a_x) const;

	public:
		// Constructors.
		Solution_dg_linear(Mesh* const &a_mesh, Solution_dg_linear* const &a_solution);
		Solution_dg_linear(Mesh* const &a_mesh, f_double const &a_f);

		// Solvers.
		void Solve(const double &a_cgTolerance);

		// Computers.
		double compute_energyNormDifference2(f_double const &a_u, f_double const &a_u_1) const;
		double compute_errorIndicator(const double &a_i) const;

		// Getters.
		f_double get_f() const;
};

#endif