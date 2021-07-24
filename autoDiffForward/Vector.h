#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <type_traits>

using namespace std;

template<class T, int d>
class Vector {
public:
	T arr[d];

	Vector() : arr{ 0 } {};

	template<class... A, class = typename std::enable_if<sizeof... (A) == d, void>::type>
	explicit Vector(const A&... args) : arr{ args... } {};

	Vector(const T(&cpy)[d]);

	T magnitude_squared();
	T magnitude();

	T dot(const Vector<T, d>&) const;
	const T& operator[](int) const;
	T& operator[] (int);

	Vector<T, d> operator/ (const T);
	Vector<T, d>& operator/= (const T);
	Vector<T, d> operator- (const Vector&) const;
	Vector<T, d>& operator-= (const Vector&);
	Vector<T, d> operator+ (const Vector&) const;
	Vector<T, d>& operator+= (const Vector&);
	Vector<T, d> operator* (const T) const;
	Vector<T, d>& operator*=(const T);

	void mk_z();
	Vector<T, 3> cross(const Vector<T, 3>&);
};

template<class T, int d>
Vector<T, d>::Vector(const T(&cpy)[d]) {
	for (int i = 0; i < d; i++) {
		arr[i] = cpy[i];
	}
}

template<class T, int d>
T Vector<T, d>::dot(const Vector<T, d>& rhs) const {
	T sum = 0;
	for (int i = 0; i < d; i++) {
		sum += (rhs[i] * arr[i]);
	}
	return sum;
}

template<class T, int d>
T Vector<T, d>::magnitude_squared() {
	return dot(*this);
}

template<class T, int d>
T Vector<T, d>::magnitude() {
	return sqrt(magnitude_squared());
}

template<class T, int d>
const T& Vector<T, d>::operator[](int index) const {
	return arr[index];
}

template<class T, int d>
T& Vector<T, d>::operator[] (int index) {
	return arr[index];
}

template<class T, int d>
Vector<T, d> Vector<T, d>::operator/ (const T rhs) {
	Vector<T, d> newArr;
	for (int i = 0; i < d; i++) {
		newArr[i] = arr[i] / rhs;
	}
	return newArr;
}

template<class T, int d>
Vector<T, d>& Vector<T, d>::operator/= (const T rhs) {
	for (int i = 0; i < d; i++) {
		arr[i] /= rhs;
	}
	return *this;
}

template<class T, int d>
Vector<T, d> Vector<T, d>::operator- (const Vector& v) const {
	Vector<T, d> newArr;
	for (int i = 0; i < d; i++) {
		newArr[i] = arr[i] - v[i];
	}
	return newArr;
}

template<class T, int d>
Vector<T, d>& Vector<T, d>::operator-= (const Vector& v) {
	for (int i = 0; i < d; i++) {
		arr[i] = arr[i] - v.arr[i];
	}
	return *this;
}

template<class T, int d>
Vector<T, d> Vector<T, d>::operator+ (const Vector& v) const {
	Vector<T, d> newArr;
	for (int i = 0; i < d; i++) {
		newArr[i] = arr[i] + v[i];
	}
	return newArr;
}

template<class T, int d>
Vector<T, d>& Vector<T, d>::operator+= (const Vector& v) {
	for (int i = 0; i < d; i++) {
		arr[i] = arr[i] + v.arr[i];
	}
	return *this;
}


template<class T, int d>
Vector<T, d>& Vector<T, d>::operator*= (const T rhs) {
	for (int i = 0; i < d; i++) {
		arr[i] = arr[i] * rhs;

	}
	return *this;
}

template<class T, int d>
Vector<T, d> Vector<T, d>::operator* (const T right) const {
	Vector<T, d> newArr;
	for (int i = 0; i < d; i++) {
		newArr[i] = arr[i] * right;
	}
	return newArr;
}

template<class T, int d>
void Vector<T, d>::mk_z() {
	*this *= 0;
}

template<class T, int d>
Vector<T, 3> Vector<T, d>::cross(const Vector<T, 3>& v) {
	Vector<T, 3> modVector(arr[1] * v[2] - arr[2] * v[1], arr[2] * v[0] - arr[0] * v[2], arr[0] * v[1] - arr[1] * v[0]);
	return modVector;
}

template<class T, int d> //You do not need two template declarations here because this is not a member function
Vector<T, d> operator* (const T left, const Vector<T, d>& right) {
	return right * left;
}

template<class T, int d>
ostream& operator<< (ostream& out, const Vector<T, d>& v) {
	out << "<";
	for (unsigned int i = 0; i < d - 1; i++) {
		out << v[i] << ", ";
	}
	out << v[d - 1] << ">";
	return out;
}

//Vector<float, 1> bar(float x) {return x;}

#endif
