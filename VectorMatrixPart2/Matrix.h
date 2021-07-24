#ifndef MATRIX_H
#define MATRIX_H

#include <math.h>
#include <vector>
#include "Vector.h"
//#include <array>

using namespace std;

template<class T, int m, int n>
class Matrix {
public:
	T arr[m][n];

	Matrix() : arr{0} {};

	template<class... A, class = typename std::enable_if<sizeof... (A) == m*n, void>::type>
	explicit Matrix(const A&... args) : arr{ args... } {};


	Matrix(const Matrix&);

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
	Matrix<T, m, n>& operator=(const Matrix&);

	template<int x>
	Matrix<T, m, x> operator*(const Matrix<T, n, x>&) const;
	Matrix<T, m, n> operator+ (const T) const;
	Matrix<T, m, n>& operator+= (const T);
	Matrix<T, m, n> operator- (const T) const;
	Matrix<T, m, n>& operator-= (const T);
	Matrix<T, n, m> t() const;
	Matrix<T, m, m> mk_id();
	void mk_z();

	Vector<T, n> transpose_times(const Vector<T, m>&) const;
	Matrix<T, m, n> transpose() const;
	template<int x>
	Matrix<T, n, x> transpose_times(const Matrix<T, m, x>&) const;
	template<int x>
	Matrix<T, m, x> times_transpose(const Matrix<T, x, n>&) const;

	void factor_lu(Matrix<T, m, m>&, Vector<int, m>&);
	Vector<T, m> forward_substitution(Vector<T, m>&);
	Vector<T, m> backward_substitution(Vector<T, m>&);
	T determinant_triangular() const;
	T determinant_cholesky() const;
	Matrix<T, m, m> inverse_lu() const;
	Matrix<T, m, m> inverse_cholesky() const;
	void factor_cholesky(Matrix<T, m, m>& ) const;
};

template<class T, int m, int n>
Matrix<T, m, n>::Matrix(const Matrix& cpy) {
	for (int row = 0; row < m; row++) {
		for (int col = 0; col < n; col++) {
			arr[row][col] = (T)cpy.arr[row][col];
		}
	}
}

template<class T, int m, int n>
const T& Matrix<T, m, n>::operator() (int row, int col) const {
	return arr[row][col];
}

template<class T, int m, int n>
T& Matrix<T, m, n>::operator() (int row, int col) {
	return arr[row][col];
}

template<class T, int m, int n>
Matrix<T, m, n> Matrix<T, m, n>::operator+(const Matrix& rhs) const {
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
	*this = *this + rhs;
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
	*this = *this - rhs;
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
Matrix<T, m, n>& Matrix<T, m, n>::operator=(const Matrix& rhs) {
	for (int row = 0; row < m; row++) {
		for (int col = 0; col < n; col++) {
			arr[row][col] = rhs.arr[row][col];
		}
	}
	return *this;
}

template<class T, int m, int n>
Matrix<T, m, n> Matrix<T, m, n>::operator+ (const T rhs) const{
	Matrix<T, m, n> id(*this);
	for (int row = 0; row < m; row++) {
		id.arr[row][row] += rhs;
	}
	return id;
}

template<class T, int m, int n>
Matrix<T, m, n>& Matrix<T, m, n>::operator+= (const T rhs) {
	for (int row = 0; row < m; row++) {
		arr[row][row] += rhs;
	}
	return *this;
}

template<class T, int m, int n>
Matrix<T, m, n> Matrix<T, m, n>::operator- (const T rhs) const{
	Matrix<T, m, n> id(*this);
	for (int row = 0; row < m; row++) {
		id.arr[row][row] -= rhs;
	}
	return id;
}

template<class T, int m, int n>
Matrix<T, m, n>& Matrix<T, m, n>::operator-= (const T rhs) {
	for (int row = 0; row < m; row++) {
		arr[row][row] -= rhs;
	}
	return *this;
}

template<class T, int m, int n>
Matrix<T, n, m> Matrix<T, m, n>::t() const{
	Matrix<T, n, m> modMatrix;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			modMatrix(j,i) = arr[i][j];
		}
	} 
	return modMatrix;
}

template<class T, int m, int n>
Matrix<T, m, m> Matrix<T, m, n>::mk_id() {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			if (i == j) {
				arr[i][j] = 1;
			}
			else {
				arr[i][j] = 0;
			}
		}
	}
	return *this;
}

template<class T, int m, int n>
void Matrix<T, m, n>::mk_z() {
	//std::fill(&arr[0][0], &arr[0][0] + sizeof(arr), 0);
	memset(arr, 0, sizeof(arr[0][0]) * m * n);
}

template<class T, int m, int n>
template<int x> //You need two template declarations because each construct works on a different template argument
Matrix<T, m, x> Matrix<T, m, n>::operator*(const Matrix<T, n, x>& rhs) const {
	Matrix<T, m, x> modMatrix;
	T sum = 0;

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < x; j++) {
			for (int k = 0; k < n; k++) {
				sum += arr[i][k] * rhs.arr[k][j];
			}
			modMatrix(i, j) = sum;
			sum = 0;
		}
	}
	return modMatrix;
}

template<class T, int m, int n>
Vector<T, m> operator* (const Matrix<T, m, n>& lhs, const Vector<T, n>& rhs) {
	Vector<T, m> mod;
	T sum = 0;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			sum += lhs(i, j) * rhs[j];
		}
		mod[i] = sum;
		sum = 0;
	}
	return mod;
}

template<class T, int m, int n>
template<int x>
Matrix<T, n, x> Matrix<T, m, n>::transpose_times(const Matrix<T, m, x>& rhs) const{
	Matrix<T, n, x> modMatrix;
	T sum = 0;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < x; j++) {
			for (int k = 0; k < m; k++) {
				sum += arr[k][i] * rhs(k, j);
			}
			
			modMatrix(i, j) = sum;
			sum = 0;
		}	
	}
	return modMatrix;
}

template<class T, int m, int n>
template<int x>
Matrix<T, m, x> Matrix<T, m, n>::times_transpose(const Matrix<T, x, n>& rhs) const {
	Matrix<T, m, x> modMatrix;
	T sum = 0;

	for (int i = 0; i < x; i++) {
		for (int j = 0; j < m; j++) {
			for (int k = 0; k < n; k++) {
				sum += arr[j][k] * rhs(i, k);
			}
			modMatrix(j, i) = sum;
			sum = 0;
		}
	}
	return modMatrix;
}

template<class T, int m, int n>
Vector<T, n> Matrix<T, m, n>::transpose_times(const Vector<T, m>& rhs) const{
	Vector<T, n> modVector;
	T sum = 0;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			sum += arr[j][i] * rhs[j];
		}
		modVector[i] = sum;
		sum = 0;
	}
	return modVector;
}

template<class T, int m, int n>
Matrix<T, m, n> Matrix<T, m, n>::transpose() const {
	Matrix<T, n, m> transposeM;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			transposeM.arr[j][i] = arr[i][j];
		}
	}
	return transposeM;
}

template<class T, int m, int n>
Matrix<T, m, n> operator* (const T left, const Matrix<T, m, n>& rhs) {
	return rhs * left;
}

template<class T, int m, int n>
ostream& operator<< (ostream& out, const Matrix<T, m, n>& rhs) {
	for (int row = 0; row < m; row++) {
		out << "  ";
		for (int col = 0; col < n - 1; col++) {
			out << rhs.arr[row][col] << "  ";
		}
		out << rhs.arr[row][n - 1];
		if (row != m - 1) {
			out << endl;
		}
	}
	return out;
}

template<class T, int m, int n>
Matrix<T, m, n> outer_product(const Vector<T, m>& u, const Vector<T, n>& v) {
	Matrix<T, m, n> modMatrix;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			modMatrix(i, j) = u[i] * v[j];
		}
	}
	return modMatrix;
}

template<class T, int m>
Vector<T, m> solve_lu(const Matrix<T, m, m>& M, const Vector<T, m>& v) {
	Matrix<T, m, m> matrixU(M);
	Vector<T, m> vectorB(v);
	Vector<int, m> P;

	Matrix<T, m, m> matrixL;
	matrixU.factor_lu(matrixL, P);

	for (int i = 0; i < m; i++) {
		if (i != P[i]) {
			vectorB.arr[i] = v.arr[P[i]];
		}
	}

	Vector<T, m> vectorY = matrixL.forward_substitution(vectorB);
	return matrixU.backward_substitution(vectorY);
}

template<class T, int m, int n>
void Matrix<T, m, n>::factor_lu(Matrix<T, m, m>& matrixL, Vector<int,m>& P) {
	//Equation = PA = LU

	for (int i = 0; i < m; i++) {
		P[i] = i;
	}

	matrixL.mk_id(); 
	
	for (int i = 0; i < m; i++) {
		int largestIndex = i;
		for (int j = i + 1; j < m; j++) {
			//Find the largest index for partial pivot
			if (abs(arr[j][i]) > abs(arr[largestIndex][i])) {
				largestIndex = j;
			}
		}
		if (largestIndex != i) {
			//Swap Rows
			swap(P[i],P[largestIndex]);
			for (int j = 0; j < m; j++) {
				swap(arr[i][j], arr[largestIndex][j]);
			}
		}
		//Eliminate
		for (int k = i + 1; k < m; k++) {
			T factor = arr[k][i] / arr[i][i];
			for (int l = 0; l < m; l++) {
				arr[k][l] -= factor * arr[i][l];
			}
			matrixL.arr[k][i] = /*-1 **/ factor;
		}
	}
}

template<class T, int m, int n>
Vector<T, m> Matrix<T, m, n>::forward_substitution(Vector<T, m>& modVector) {
	Vector<T, m> finalVector;
	
	for (int i = 0; i < m; i++) {
		int k = 0;
		double sum = 0;
		for (int j = 0; j < i; j++) {
			 sum += arr[i][j] * finalVector.arr[k];
			++k;
		}
		finalVector.arr[i] = (modVector[i] - sum) / arr[i][i];
	}
	return finalVector;
}

template<class T, int m, int n>
Vector<T, m> Matrix<T, m, n>::backward_substitution(Vector<T, m>& modVector) {
	Vector<T, m> finalVector;
	
	for (int i = m - 1; i >= 0; i--) {
		int k = m - 1;
		double sum = 0;
		for (int j = m - 1; j >= i; j--) {
			sum += arr[i][j] * finalVector.arr[k];
			--k;
		}
		finalVector.arr[i] = (modVector[i] - sum) / arr[i][i];
	}
	return finalVector;
}

template<class T, int m, int n>
T Matrix<T, m, n>::determinant_triangular() const {
	Matrix<T, m, m> matrixU(*this);
	Matrix<T, m, m> matrixL;
	Vector<int, m> P;
	
	matrixU.factor_lu(matrixL, P);

	int rowSwap = 0;
	for (int i = 0; i < m; i++) {
		if (P.arr[i] != i) {
			rowSwap++;
		}
	}

	T product = 1;
	for (int i = 0; i < m; i++) {
		product *= matrixU.arr[i][i];
	}

	if (rowSwap > 0 && rowSwap%2 == 0) {
		return -1 * product;
	}
	return product;
}

template<class T, int m, int n>
T Matrix<T, m, n>::determinant_cholesky() const {
	Matrix<T, m, m> matrixM(*this);
	Matrix<T, m, m> matrixL;

	matrixM.factor_cholesky(matrixL);

	T product = 1;
	for (int i = 0; i < m; i++) {
		product *= matrixL.arr[i][i];
	}
	return product * product;
}

template<class T, int m, int n>
Matrix<T, m, m> Matrix<T, m, n>::inverse_lu() const{
	//matrixL.inverse_lu(matrixU, rowSwap, P)
	Matrix<T, m, m> U(*this);
	Matrix<T, m, m> L;
	Vector<int, m> P;
	U.factor_lu(L, P);

	Matrix<T, m, m> matrixU;
	Matrix<T, m, m> matrixL;
	matrixU.mk_id();
	matrixL.mk_id();

	for (int i = m - 1; i >= 0; i--) {
		for (int j = m - 1; j >= i; j--) {
			if (j == i) {
				for (int k = 0; k < m; k++) {
					matrixU.arr[i][k] /= U.arr[i][j];
				}
			}
			else {
				for (int k = 0; k < m; k++) {
					matrixU.arr[i][k] -= (matrixU.arr[j][k] * U.arr[i][j]);
				}
			}
		}
	}

	for (int i = 1; i < m; i++) {
		for (int j = 0; j < i; j++) {
			for (int k = 0; k < m; k++) {
				matrixL.arr[i][k] -= (matrixL.arr[j][k] * L.arr[i][j]);
			}

		}
	}

	int rowSwap = 0;
	for (int i = 0; i < m; i++) {
		if (P[i] != i) {
			rowSwap++;
		}
	}

	if (rowSwap > 0) {
		return matrixU * matrixL;
	}
	else {
		Matrix<T, m, m> matrixO(matrixU * matrixL);
		Matrix<T, m, m> matrixI;

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				matrixI.arr[j][i] = matrixO.arr[j][P[i]];
			}	
		}
		return matrixI;
	}
}

template<class T, int m, int n>
Matrix<T, m, m> Matrix<T, m, n>::inverse_cholesky() const{
	Matrix<T, m, m> matrixM(*this);
	Matrix<T, m, m> matrixL;
	Matrix<T, m, m> finalM;

	finalM.mk_id();
	matrixM.factor_cholesky(matrixL);
	
	for (int i = 0; i < m; i++) {
		for (int j = 0; j <= i; j++) {
			if (j == i) {
				for (int k = 0; k < m; k++) {
					finalM.arr[i][k] /= matrixL.arr[i][j];
				}
			}
			else {
				for (int k = 0; k < m; k++) {
					finalM.arr[i][k] -= (finalM.arr[j][k] * matrixL.arr[i][j]);
				}
			}
		}
	}
	return finalM.transpose_times(finalM);
}

template<class T, int m>
Vector<T, m> solve_cholesky(const Matrix<T, m, m>& M, const Vector<T, m>& v) {
	Vector<T, m> B(v);
	Matrix<T, m, m> matrixM(M);
	Matrix<T, m, m> matrixL;

	matrixM.factor_cholesky(matrixL);

	Vector<T, m> Y = matrixL.forward_substitution(B);
	return matrixL.transpose().backward_substitution(Y);
}

template<class T, int m, int n>
void Matrix<T, m, n>::factor_cholesky(Matrix<T, m, m>& matrixL) const{
	matrixL.mk_id();

	for (int j = 0; j < m; j++) {
		for (int i = j; i < m; i++) {
			T sum = 0;
			for (int k = j; k > 0; k--) {
				sum += matrixL.arr[i][k - 1] * matrixL.arr[j][k - 1];
			}
			if (i == j) {
				matrixL.arr[i][j] = sqrt(arr[i][j] - sum);
			}
			else {
				matrixL.arr[i][j] = (arr[i][j] - sum) / matrixL.arr[j][j];
			}
		}
	}
}

template<class T, int m>
Vector<T, m> solve(const Matrix<T, m, m>& M, const Vector<T, m>& u) {
	Matrix<T, m, m> modMatrix(M);
	Vector<T, m> modVector(u);
	Vector<T, m> finalVector;
	double factor = 0.0;
	double sum = 0.0;


	for (int i = 0; i < m; i++) {
		for (int k = i + 1; k < m; k++) {
			factor = (double)modMatrix.arr[k][i] / modMatrix.arr[i][i];
			if (factor != 0) {
				for (int l = 0; l < m; l++) {
					modMatrix.arr[k][l] -= factor * modMatrix.arr[i][l];
				}
				modVector[k] -= factor * modVector[i];
			}
		}
	}

	int k = m - 1; 
	for (int i = m - 1; i >= 0; i--) {
		for (int j = m - 1; j >= i; j--) {
			sum += modMatrix.arr[i][j] * finalVector.arr[k];
			--k;
		}
		finalVector.arr[i] = (modVector[i] - sum) / modMatrix.arr[i][i];
		sum = 0;
		k = m - 1;
	}
	return finalVector;
}

#endif
