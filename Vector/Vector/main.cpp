#include "Vector.h"
#include "Matrix.h"
using namespace std;

int main() {
	cout << "VECTOR: " << endl;
	Vector<double, 4> u, v;
	Vector<int, 7> t(2, 1, 2, 2, 2, 2, 2);
	//Vector<int, 1> r(7, 6);
	//Vector<float, 3> r(7, 6);
	//Vector<int, 1> foo() { return 7; }

	for (int i = 0; i < 4; i++) u[i] = v[i] = i + 1;
	cout << "Original V: " << v << endl;
	cout << "Original U: " << u << endl;
	cout << endl;
	cout << u.dot(v) * u + v * 2.1 - 5.2 * u << endl;
	cout << u.magnitude() + v.magnitude_squared() << endl;
	double r = u[0] - u[2];
	cout << "u[0]: " << u[0] << ", u[2]: " << u[2] << endl;
	cout << "Subtracting u[0] & u[2]: " << r << endl;

	u *= 4.2;
	v /= 3;
	cout << "U multiplied by 4.2: " << u << endl;
	cout << "V divided by 3: " << v << endl;
	u -= v / 7.4;
	cout << "U subtract V (divided by 7.4): " << u << endl;
	cout << "3 times U: " << 3 * u << endl;
	v += u * 3;
	cout << "V plus 3 times U: " << v << endl;

	Vector<Vector<int, 2>, 3> z;
	std::cout << "Z: " << z << std::endl;

	cout << endl;
	cout << "MATRIX: " << endl;
	Matrix<int, 1, 3> al(3, 4, 2);
	cout << "al: " << endl;
	cout << al << endl;
	Matrix<int, 3, 4> all(13, 9, 7, 15, 8, 7, 4, 6, 6, 4, 0, 3);
	Matrix<double, 2, 2> k({{ 1.5, 2.3 }, { 3.2, 4.3 }});
	cout << "all: " << endl;
	cout << all << endl;
	cout << "al x all: " << endl;
	//cout << matrixM(al, all) << endl;
	cout << "k: " << endl;
	cout << k << endl;

	Matrix<double, 2, 2> i(1.0, 2.5, 3, 4);
	Matrix<double, 2, 2> j(5, 6, 7.4, 8);
	cout << "j: " << endl;
	cout << j << endl;
	cout << "i: " << endl;
	cout << i << endl;
	cout << "j - i:" << endl;
	cout << j - i << endl;
	cout << "j + i:" << endl;
	cout << j + i << endl;
	cout << "i * 3:" << endl;
	cout << i * 3 << endl;
	cout << "i / 3:" << endl;
	cout << i / 3 << endl;

	/*j -= i;
	cout << "j(-=): " << endl;
	cout << j << endl;
	j += i;
	cout << "j(+=): " << endl;
	cout << j << endl;
	j *= 4;
	cout << "j(*=): " << endl;
	cout << j << endl;
	i /= 3;
	cout << "i(/=): " << endl;
	cout << i << endl;*/

	/*
	Vector<int, 3> x, y;
	x[1] = 2;
	x[0] = 3;
	x[2] = 5;
	y[0] = y[1] = y[2] = 2;
	cout << x << " " << y << endl;
	cout << x - y << endl;
	cout << y - x << endl;
	cout << x + y << endl;
	cout << y + x << endl;
	cout << 3 * y << endl;
	cout << y / 4 << endl;
	*/


}




