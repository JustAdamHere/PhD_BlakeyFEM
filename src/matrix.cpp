/******************************************************************************
 * @details This is a file containing declarations of [Matrix].
 * 
 * @author     Adam Matthew Blakey
 * @date       2020/01/02
 ******************************************************************************/
#ifndef CLASS_SRC_MATRIX
#define CLASS_SRC_MATRIX

#include "matrix.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

/******************************************************************************
 * __get_diagonal__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
std::vector<T> Matrix<T>::get_diagonal() const
{
	int diagonalLength = this->get_noColumns()<this->get_noRows()?this->get_noColumns():this->get_noRows();
	std::vector<T> diagonal(diagonalLength, 0);

	for (int i=0; i<diagonalLength; ++i)
		diagonal[i] = item(i, i);

	return diagonal;
}

/******************************************************************************
 * __operator()__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
const T Matrix<T>::operator()(const int &a_x, const int &a_y) const
{
	return item(a_x, a_y);
}

/******************************************************************************
 * __operator=__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T> &a_RHS)
{
	if (&a_RHS == this)
		return *this;

	this->noRows = a_RHS.get_noRows();
	this->noColumns = a_RHS.get_noColumns();

	this->resize(this->get_noRows() * this->get_noColumns());

	for (int i=0; i<this->get_noColumns(); ++i)
		for (int j=0; j<this->get_noRows(); ++j)
			item(i, j) = a_RHS(i, j);

	return *this;
}

/******************************************************************************
 * __operator+=__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T> &a_RHS)
{
	// Dimensions must be the same.
	if (this->get_noRows() != a_RHS.get_noRows() || this->get_noColumns() != a_RHS.get_noColumns())
	{
		std::cerr << "Matrix dimensions do not match.";
		return *this;
	}

	// Creates new matrix and calculates elements appropriately.
	for (int i=0; i<this->get_noColumns(); ++i)
		for (int j=0; j<this->get_noRows(); ++j)
			item(i, j) += a_RHS(i, j); // I think it's a problem with this LHS -- I don't think you can write to it...

	return *this;
}

/******************************************************************************
 * __operator-=__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T> &a_RHS)
{
	// Dimensions must be the same.
	if (this->get_noRows() != a_RHS.get_noRows() || this->get_noColumns() != a_RHS.get_noColumns())
	{
		std::cerr << "Matrix dimensions do not match.";
		return *this;
	}

	// Creates new matrix and calculates elements appropriately.
	for (int i=0; i<this->get_noColumns(); ++i)
		for (int j=0; j<this->get_noRows(); ++j)
			item(i, j) -= a_RHS(i, j); // I think it's a problem with this LHS -- I don't think you can write to it...

	return *this;
}

/******************************************************************************
 * __operator*=__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T> &a_RHS)
{
	// Matching dimensions.
	if (this->get_noRows() != a_RHS.get_noColumns() || this->get_noRows() != a_RHS.get_noRows())
	{
		std::cerr << "Error: Matrix dimensions do not match." << std::endl;
		return *this;
	}

	Matrix<T> tempMatrix(this->get_noColumns(), this->get_noRows(), 0);

	// Creates new matrix and calculates elements appropriately.
	for (int i=0; i<this->get_noColumns(); ++i)
		for (int j=0; j<this->get_noRows(); ++j)
			for (int k=0; k<this->get_noRows(); ++k)
				tempMatrix(i, j) += item(i, k) * a_RHS(k, j);

	*this = tempMatrix;

	return *this;
}

/******************************************************************************
 * __operator*=__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix<T>& Matrix<T>::operator*=(const T &a_RHS)
{
	// Creates new matrix and calculates elements appropriately.
	for (int i=0; i<this->get_noColumns(); ++i)
		for (int j=0; j<this->get_noRows(); ++j)
			item(i, j) *= a_RHS;

	return *this;
}

/******************************************************************************
 * __operator*__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
std::vector<T> Matrix<T>::operator*(const std::vector<T> &a_RHS) const
{
	if (get_noColumns() != a_RHS.size())
	{
		std::cerr << "Error: Matrix-vector dimensions do not match." << std::endl;
		return a_RHS;
	}

	std::vector<T> tempVector(get_noColumns(), 0);
	for (int i=0; i<get_noColumns(); ++i)
		for (int j=0; j<get_noRows(); ++j)
			tempVector[i] += item(i, j) * a_RHS[j];

	return tempVector;
}

template<class T>
std::vector<T> operator*(const Matrix<T> &a_matrix, const std::vector<T> &a_vector)
{
	return a_matrix*a_vector;
}


/******************************************************************************
 * __operator/=__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix<T>& Matrix<T>::operator/=(const T &a_RHS)
{
	// Creates new matrix and calculates elements appropriately.
	for (int i=0; i<this->get_noColumns(); ++i)
		for (int j=0; j<this->get_noRows(); ++j)
			item(i, j) /= a_RHS;

	return *this;
}

/******************************************************************************
 * __argmax__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
std::vector<int> Matrix<T>::argmax(const int &i_0, const int &i_1, const int &j_0, const int &j_1)
{
	std::vector<int> args  = {i_0, j_0};
	T 				 M_max = item(i_0, j_0);

	for (int i=i_0; i<i_1; ++i)
		for (int j=j_0; j<j_1; ++j)
			if(T value = item(i, j) > M_max)
			{
				M_max = value;
				args[1] = j;
				args[0] = i;
			}

	return args;
}

/******************************************************************************
 * __argmax__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
void Matrix<T>::swapRows(const int &a_row1, const int &a_row2)
{
	std::vector<T> row1(get_noColumns());
	std::vector<T> row2(get_noColumns());

	for (int i=0; i<get_noColumns(); ++i)
	{
		row1[i] = item(i, a_row1);
		row2[i] = item(i, a_row2);
	}

	for (int i=0; i<get_noColumns(); ++i)
	{
		set(i, a_row1, row2[i]);
		set(i, a_row2, row1[i]);
	}
}

#endif