#pragma once
#include <iostream>

#ifndef MATRIX_H
#define MATRIX_H

using namespace std;

template<class T, int m, int n>
class Matrix {
	public:
		T arr[m][n];

		template<class... A>
		Matrix(const A&... args) : arr{ (T)args... } {};
		Matrix(const T (&cpy)[m][n]);

		const T& operator()(int, int) const;
		T& operator()(int, int);
		Matrix<T, m, n> operator+(const Matrix&) const;
		Matrix<T, m, n>& operator+=(const Matrix&);
		Matrix<T, m, n> operator-(const Matrix&) const;
		Matrix<T, m, n>& operator-=(const Matrix&);
		Matrix<T, m, n> operator/ (const T) const;
		Matrix<T, m, n>& operator/= (const T);
		Matrix<T, m, n> operator* (const T) const;
		Matrix<T, m, n>& operator*= (const T);

		template<int x>
		Matrix<T, m, x> operator*(const Matrix<T, n, x>& rhs) const;
};

template<class T, int m, int n>
Matrix<T, m, n>::Matrix(const T (&cpy)[m][n]) {
	for (int row = 0; row < m; row++) {
		for (int col = 0; col < n; col++) {
			arr[row][col] = cpy[row][col];
		}
	}
}

/*template<class T, int m, int n>
typename Matrix<T, m, n>::Row Matrix<T, m, n>::operator[](int index) {
	return Row(*this, index);
}*/

template<class T, int m, int n>
const T& Matrix<T, m, n>::operator() (int row, int col) const{
	return arr[row][col];
}

template<class T, int m, int n>
T& Matrix<T, m, n>::operator() (int row, int col) {
	return arr[row][col];
}

/*template<class T, int m, int n>
T& Matrix<T, m, n>::Row::operator[](int index) {
	return parent.arr[rowNumber][index];
}

template<class T, int m, int n>
const T& Matrix<T, m, n>::Row::operator[](int index) const{
	return parent.arr[rowNumber][index];
}*/

template<class T, int m, int n>
Matrix<T, m, n> Matrix<T, m, n>::operator+(const Matrix& rhs) const{
	Matrix modMatrix;
	for (int row = 0; row < m; row++) {
		for (int col = 0; col < n; col++) {
			modMatrix(row, col) = arr[row][col] + rhs.arr[row][col];
		}
	}
	return modMatrix;
}

template<class T, int m, int n>
Matrix<T, m, n>& Matrix<T, m, n>::operator+=(const Matrix& rhs) {
	for (int row = 0; row < m; row++) {
		for (int col = 0; col < n; col++) {
			arr[row][col] += rhs.arr[row][col];
		}
	}
	return *this;
}

template<class T, int m, int n>
Matrix<T, m, n> Matrix<T, m, n>::operator-(const Matrix& rhs) const {
	Matrix modMatrix;
	for (int row = 0; row < m; row++) {
		for (int col = 0; col < n; col++) {
			modMatrix(row, col) = arr[row][col] - rhs.arr[row][col];
		}
	}
	return modMatrix;
}

template<class T, int m, int n>
Matrix<T, m, n>& Matrix<T, m, n>::operator-=(const Matrix& rhs) {
	for (int row = 0; row < m; row++) {
		for (int col = 0; col < n; col++) {
			arr[row][col] -= rhs.arr[row][col];
		}
	}
	return *this;
}

template<class T, int m, int n>
Matrix<T, m, n> Matrix<T, m, n>::operator/ (const T num) const {
	Matrix modMatrix;
	for (int row = 0; row < m; row++) {
		for (int col = 0; col < n; col++) {
			modMatrix(row, col) = arr[row][col] / num;
		}
	}
	return modMatrix;
}

template<class T, int m, int n>
Matrix<T, m, n>& Matrix<T, m, n>::operator/=(const T num) {
	for (int row = 0; row < m; row++) {
		for (int col = 0; col < n; col++) {
			arr[row][col] /= num;
		}
	}
	return *this;
}

template<class T, int m, int n>
Matrix<T, m, n> Matrix<T, m, n>::operator* (const T num) const {
	Matrix modMatrix;
	for (int row = 0; row < m; row++) {
		for (int col = 0; col < n; col++) {
			modMatrix(row, col) = arr[row][col] * num;
		}
	}
	return modMatrix;
}

template<class T, int m, int n>
Matrix<T, m, n>& Matrix<T, m, n>::operator*=(const T num) {
	for (int row = 0; row < m; row++) {
		for (int col = 0; col < n; col++) {
			arr[row][col] *= num;
		}
	}
	return *this;
}


template<class T, int m, int n>
template<int x> //You need two template declarations because each construct works on a different template argument
Matrix<T, m, x> Matrix<T, m, n>::operator*(const Matrix<T, n, x>& rhs) const{
	Matrix<T, m, x> modMatrix;
	T sum = 0;

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < x; j++) {
			modMatrix(i, j) = 0;
			for (int k = 0; k < n; k++) {
				sum += arr[i][k] * rhs.arr[k][j];
			}
			modMatrix(i, j) = sum;
			sum = 0;
		}
	}
	return modMatrix;
}

template<class T, int m, int n, class W>
Matrix<T, m, n> operator* (const W left, const Matrix<T, m, n>& rhs) {
	return rhs * left;
}

template<class T, int m, int n>
ostream& operator<< (ostream& out, const Matrix<T, m, n>& rhs) {
	for (int row = 0; row < m; row++) {
		out << "  ";
		for (int col = 0; col < n-1; col++) {
			out << rhs.arr[row][col] << "  ";
		}
		out << rhs.arr[row][n - 1];
		if (row != m - 1) {
			out << endl;
		}
	}
	return out;
}
#endif


