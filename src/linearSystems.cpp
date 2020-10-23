/******************************************************************************
 * LINEARSYSTEMS.CPP
 *
 * This is a file containing functions regarding quadratures.
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/11/4
 ******************************************************************************/

#include "linearSystems.hpp"
#include "matrix.hpp"
#include "matrix_full.hpp"
#include <cassert>
#include <vector>

#include <iostream>

/******************************************************************************
 * thomasInvert
 * 
 * @details    Inverts a tridiagonal matrix.
 *
 * @param[in] n 			Gives the degree of the polynomial.
 * @param[in] x 			The point at which to evaluate the polynomial.
 ******************************************************************************/
namespace linearSystems
{
	void thomasInvert(const std::vector<double> a, const std::vector<double> b, const std::vector<double> c, const std::vector<double> d, std::vector<double> &solution)
	{
		int n = b.size();

		assert(n >= 2);

		std::vector<double> c_(n-1);
		std::vector<double> d_(n);

		c_[0] = c[0]/b[0];
		for (int i=1; i<n-1; ++i)
		{
			c_[i] = c[i]/(b[i] - c_[i-1]*a[i-1]);
		}

		d_[0] = d[0]/b[0];
		for (int i=1; i<n; ++i)
		{
			d_[i] = (d[i] - d_[i-1]*a[i-1])/(b[i] - c_[i-1]*a[i-1]);
		}

		solution[n-1] = d_[n-1];
		for (int i=n-2; i>=0; --i)
		{
			solution[i] = d_[i] - c_[i]*solution[i+1];
		}
	}

	std::vector<double> conjugateGradient(const Matrix<double> &a_M, const std::vector<double> &a_b, const double &a_tolerance)
	{
		std::vector<double> x(a_M.get_noColumns(), 0);

		std::vector<double> r = a_b;
		std::vector<double> p = a_b;

		double r_2 = dotProduct(r, r);
		double errorNorm = sqrt(r_2);

		int noIterations = 0;

		while(errorNorm > a_tolerance)
		{
			std::vector<double> pNew = a_M * p;
			double alpha = r_2/dotProduct(p, pNew);
			x += alpha * p;
			r += -alpha * pNew;
			double r_2New = dotProduct(r, r);
			double beta = r_2New/r_2;

			r_2 = r_2New;
			p = r + beta*p;
			errorNorm = sqrt(r_2);
			++noIterations;
		}

		return x;
	}

	// Modified from:
	// http://www.cplusplus.com/forum/beginner/124448/
	std::vector<double> GaussJordan(const Matrix<double> &a_M, const std::vector<double> &a_b)
	{
		assert(a_M.get_noColumns() == a_M.get_noRows());

		Matrix_full<double> test(a_M);

		int n = a_M.get_noColumns();
		Matrix_full<double> a(n, 2*n, 0);

		// Sets coefficients of matrix correctly.
		for (int i=0; i<n; ++i)
			for (int j=0; j<n; ++j)
				a.set(i, j, a_M(i, j));

		for (int i=0; i<n; ++i)
			a.set(i, i+n, 1);

		// Partial pivoting.
		/*for (int i=n-1; i>0; --i)
		{
			if (a(i-1, 1) < a(i, 1))
				for (int j=0; j<2*n; ++j)
				{
					double d = a(i, j);
					a.set(i,   j, a(i-1, j));
					a.set(i-1, j, d);
				}
		}
		std::cout << "pivoted output: " << std::endl;
		for (int i=0; i<n; ++i)
		{
			for (int j=0; j<2*n; ++j)
				std::cout << a(i, j) << "    ";
			std::cout << std::endl;
		}*/
		
		// Reducing to diagonal matrix.
		for (int i=0; i<n; ++i)
		{
			for (int j=0; j<n; ++j)
				if (j != i)
				{
					double d = a(j, i) / a(i, i);
					for (int k=0; k<2*n; ++k)
						a.set(j, k, a(j, k) - a(i, k)*d);
				}
		}

		// Reducing to unit matrix.
		for (int i=0; i<n; ++i)
		{
			double d = a(i, i);
			for (int j=0; j<2*n; ++j)
				a.set(i, j, a(i, j)/d);
		}

		// Setting the inverse.
		Matrix_full<double> M_(n, n);

		for (int i=0; i<n; ++i)
			for (int j=0; j<n; ++j)
				M_.set(i, j, a(i, n+j));

		std::vector<double> x = M_*a_b;

		return x;
	}

	std::vector<double> direct(const Matrix_full<double> &a_M, const std::vector<double> &a_b)
	{
		int n = a_M.get_noRows();
		Matrix_full<double> M_ = a_M.calculate_adjugate();

		//M_ /= a_M.calculate_determinant();

		std::vector<double> x = M_*a_b;

		return x;
	}

	double dotProduct(const std::vector<double> &a_v1, const std::vector<double> &a_v2)
	{
		double result = 0;

		for (int i=0; i<a_v1.size(); ++i)
			result += a_v1[i] * a_v2[i];

		return result;
	}
}

std::vector<double>& operator+=(std::vector<double> &a_v1, const std::vector<double> &a_v2)
{
	for (int i=0; i<a_v1.size(); ++i)
		a_v1[i] += a_v2[i];

	return a_v1;
}

std::vector<double> operator+(const std::vector<double> &a_v1, const std::vector<double> &a_v2)
{
	std::vector<double> result(a_v1.size());

	for (int i=0; i<a_v1.size(); ++i)
		result[i] = a_v1[i] + a_v2[i];

	return result;
}