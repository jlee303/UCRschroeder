#include "Matrix.h"
#include <vector>

using namespace std;

int main() {

	cout << "1. " << endl;
	Matrix<double, 3, 3> m(25.0, 15.0, -5.0, 15.0, 18.0, 0.0, -5.0, 0.0, 11.0);
	Vector<double, 3> v(15.0, 15.0, 17.0);
	cout << "Determinant: " << m.determinant_cholesky() << endl;
	cout << "Inverse: " << endl;
	cout << m.inverse_cholesky() << endl;
	cout << solve_cholesky(m, v) << endl;
	cout << endl;
	cout << "Inverse Check: " << endl;
	cout << m.inverse_cholesky() * m - 1 << endl;
	cout << "Cholesky Check: " << m * solve_cholesky(m, v) - v << endl;
	cout << endl;

	cout << "2. " << endl;
	Matrix<double, 3, 3> m1(4.0, 2.0, 6.0, 2.0, 2.0, 5.0, 6.0, 5.0, 22.0);
	Vector<double, 3> v1(7.0, 2.0, 3.0);
	cout << "Determinant: " << m1.determinant_cholesky() << endl;
	cout << "Inverse: " << endl;
	cout << m1.inverse_cholesky() << endl;
	cout << solve_cholesky(m1, v1) << endl;
	cout << endl;
	cout << "Inverse Check: " << endl;
	cout << m1.inverse_cholesky()* m1 - 1 << endl;
	cout << "Cholesky Check: " << m1 * solve_cholesky(m1, v1) - v1 << endl;
	cout << endl;
	
	cout << "3. " << endl;
	Matrix<double, 4, 4> te(0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 4.0, 12.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 4.0, 8.0);
	Vector<double, 4> teV(-0.2, 0.8, 1.5, 1.2);
	cout << "Determinant: " << te.determinant_triangular() << endl;
	cout << "Inverse: " << endl;
	cout << te.inverse_lu() << endl;
	cout << solve_lu(te, teV) << endl << endl; //-0.8, 6, -4.9, 1.2 


	cout << "4. " << endl;
	Matrix<double, 3, 3> t(1.0, 1.0, -1.0, 1.0, -2.0, 3.0, 2.0, 3.0, 1.0);
	Vector<double, 3> tV(4.0, -6.0, 7.0);
	cout << "Determinant: " << t.determinant_triangular() << endl;
	cout << "Inverse: " << endl;
	cout << t.inverse_lu() << endl;
	cout << solve_lu(t, tV) << endl << endl; // 1, 2, -1


	cout << "5. " << endl;
	Matrix<double, 3, 3> m2(4.0, -3.0, 1.0, -2.0, 1.0, -3.0, 1.0, -1.0, 2.0);
	Vector<double, 3> v2(-8.0, -4.0, 3.0);
	cout << "Determinant: " << m2.determinant_triangular() << endl;
	cout << "Inverse: " << endl;
	cout << m2.inverse_lu() << endl;
	cout << solve_lu(m2, v2) << endl << endl; //-2, 1, 3


	cout << "6. " << endl;
	Matrix<double, 2, 2> m3(0.0, 1.0, 1.0, 0.0);
	Vector<double, 2> v3(7.0, 4.0);
	cout << "Determinant: " << m3.determinant_triangular() << endl;
	cout << "Inverse: " << endl;
	cout << m3.inverse_lu() << endl;
	cout << solve_lu(m3, v3) << endl << endl; // 4, 7
}
	



	/*
	cout << "4. " << endl;
	Matrix<double, 3, 3> m1(1.0, 3.0, 2.0, 3.0, 7.0, 5.0, 2.0, 2.0, 4.0);
	Vector<double, 3> v1(11.0, 28.0, 14.0);
	cout << "Sol: " << solve(m1, v1) << endl;
	cout << "Sol(lu): " << solve_lu(m1, v1) <<  "   Incorrect => <5/3, 28, 70/3> != <11,28,14>" << endl;
	//cout << m1 * solve_lu(m1, v1) << endl;

	cout << endl;
	cout << "5. " << endl;
	Matrix<double, 3, 3> m2(4.0, -3.0, 1.0, -2.0, 1.0, - 3.0, 1.0, -1.0, 2.0);
	Vector<double, 3> v2(-8.0, -4.0, 3.0);
	cout << "Sol: " << solve(m2, v2) << endl;
	cout << "Sol(lu): " << solve_lu(m2, v2) << endl;

	cout << endl;
	cout << "6. " << endl;
	Matrix<double, 2, 2> m3(2.0, 3.0, 3.0, -2.0);
	Vector<double, 2> v3(7.0, 4.0);
	cout << "Sol: " << solve(m3, v3) << endl;
	cout << "Sol(lu): " << solve_lu(m3, v3) << endl;

	cout << endl;
	cout << "7. " << endl;
	Matrix<double, 3, 3> m4(4.0, -2.0, 2.0, -2.0, 1.0, 3.0, 2.0, -2.0, 2.0);
	Vector<double, 3> v4(2.0, 1.0, -4.0);
	cout << "Sol: " << solve(m4, v4) << endl;
	cout << "Sol(lu): " << solve_lu(m4, v4) << "   Incorrect => <2,-1,-2> != <2, 1, -4>, should be 3, 11/2, 1/2" << endl; 
	//cout << m4 * solve_lu(m4, v4) << endl; Line 410 in Matrix.h should have -1 * factor?
	
	
	
	
	
	cout << "Testing Outer Product: " << endl;
	Vector<int, 4> v4(7, 2, 3, 1);
	Vector<int, 3> v5(3, 2, 1);
	cout << v4 << " outer product " << v5 << endl;
	cout << outer_product(v4, v5) << endl;
	cout << endl;

	cout << "Testing Scalar addition/subtraction: " << endl;
	Matrix<int, 3, 3> m4(m1);
	cout << "Original Matrix: " << endl;
	cout << m4 << endl;
	cout << endl;

	cout << "Adding 2: " << endl;
	cout << m4 + 2 << endl;
	cout << "Subtracting 2: " << endl;
	cout << m4 - 2 << endl;
	cout << endl;

	cout << "Testing Transpose Times/Times Transpose: " << endl;
	Matrix<int, 2, 3> m5(1, 2, 1, 2, 5, 2);
	Matrix<int, 2, 3> m6(1, 0, 1, 1, 1, 0);
	cout << "Original Matrix A: " << endl;
	cout << m5 << endl;
	cout << "Original Matrix B: " << endl;
	cout << m6 << endl;
	cout << endl;

	cout << "A^T*B: " << endl;
	cout << m5.transpose_times(m6) << endl;
	cout << "B*A^T: " << endl;
	cout << m6.times_transpose(m5) << endl;
	cout << endl;

	cout << "Testing Tranpose Times Vector: " << endl;
	Matrix<int, 3, 2> m7(1, 2, 3, 4, 5, 6);
	Vector<int, 3> v6(7, 9, 11);
	cout << "Original Matrix A: " << endl;
	cout << m7 << endl;
	cout << "Original Vector B: " << endl;
	cout << v6 << endl;
	cout << endl;
	cout << "A^T*B: " << endl;
	cout << m7.transpose_times(v6) << endl;
	
	Matrix<int, 2, 2> x(1, 2, 3, 4);
	Matrix<int, 2, 2> y(1, 2, 3, 4);
	Vector<int, 2> z(1, 2);
	cout << "Matrix x: " << endl;
	cout << x << endl;
	cout << "Matrix y: " << endl;
	cout << y << endl;
	cout << "Vector z" << endl;
	cout << z << endl;
	cout << endl;

	cout << "Testing add/subtract" << endl;
	cout << "Add: " << endl;
	cout << x + y << endl;
	x += y;
	cout << "Subtract: " << endl;
	cout << x - y << endl;
	x -= y;
	cout << "Testing matrix multiplication" << endl;
	cout << "x * y:" << endl;
	cout << x * y << endl;
	cout << "Testing scalar multiplication + division" << endl;
	cout << "x * 2:" << endl;
	cout << 2 * x << endl;
	x *= 2;
	cout << "x / 2:" << endl;
	cout << x / 2 << endl;
	x /= 2;
	cout << "Testing transpose" << endl;
	cout << "x:" << endl;
	cout << x.t() << endl;
	cout << "Testing make identity" << endl;
	cout << "x:" << endl;
	x.mk_id();
	cout << x << endl;
	cout << "Testing make zero" << endl;
	cout << "x:" << endl;
	x.mk_z();
	cout << x << endl;
	cout << "Testing scalar/identity matrix" << endl;
	cout << "x + 2: " << endl;
	cout << x + 2 << endl;
	x += 2;
	cout << "x - 1:" << endl;
	cout << x - 1 << endl;
	x -= 1;
	cout << "Testing Matrix * Vector" << endl;
	cout << "y * z: " << endl;
	cout << y * z << endl;
	cout << "Testing Transpose Times Vector" << endl;
	cout << "y^t * z:" << endl;
	cout << y.transpose_times(z) << endl;
	cout << "Testing Transpose Times Matrix" << endl;
	cout << "y^t * x:" << endl;
	cout << y.transpose_times(x) << endl;
	cout << "Testing Times Tranpose Matrix" << endl;
	cout << "y * x^t:" << endl;
	cout << x.times_transpose(y) << endl;
	cout << endl;

	Vector<double, 3> a(2.2, 3.7, 4.5);
	Vector<double, 3> b(5.2, 6.1, 7.4);
	cout << "Vector a: " << endl;
	cout << a << endl;
	cout << "Vector b: " << endl;
	cout << b << endl;
	cout << endl;
	cout << "Testing Cross Product" << endl;
	cout << "a x b: " << endl;
	cout << a.cross(b) << endl;
	cout << "Testing Outer product" << endl;
	cout << "a outerpruduct b" << endl;
	cout << outer_product(a, b) << endl;
	*/