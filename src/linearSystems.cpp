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

	// From:
	// https://en.wikipedia.org/wiki/Gaussian_elimination
	// It may not actually work...
	std::vector<double> GaussJordan(const Matrix<double> &a_M, const std::vector<double> &a_b)
	{
		Matrix_full<double> M_ = a_M;

		int m = M_.get_noRows();
		int n = M_.get_noColumns();

		int h = 0;
		int k = 0;
		while ((h < m) && (k < n))
		{
			// Find the kth pivot.
			int i_max = M_.argmax(h, m, k, k+1)[0];

			if (M_(i_max, k) == 0)
			{
				++k;
			}
			else
			{
				M_.swapRows(h, i_max);

				for (int i=h; i<m; ++i)
				{
					double f = M_(i, k) / M_(h, k);

					std::cout << M_(h, k) << std::endl;

					// Fill zeros with lower part of pivot column.
					M_.set(i, k, 0);
					for (int j=k; j<n; ++j)
						M_.set(i, j, M_(i, j) - M_(h, j)*f);
				}

				++h;
				++k;
			}
		}

		//std::vector<double> x(a_M.get_noColumns(), 0);
		std::vector<double> x = M_*a_b;

		/*for (int i=0; i<M_.get_noColumns(); ++i)
		{
			for (int j=0; j<M_.get_noRows(); ++j)
			{
				std::cout << M_(i, j) << " ";
			}
			std::cout << std::endl;
		}*/

		return x;
	}

	std::vector<double> direct(const Matrix<double> &a_M, const std::vector<double> &a_b)
	{

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