#include <iostream>
#include "cstdarg"
#include <initializer_list>

#ifndef VECTOR_H
#define VECTOR_H

using namespace std;

template<class T, int d>
class Vector {
	public:
		T arr[d];

		template<class... A>
		Vector(const A&... args) : arr{ args... } {}; // Why does it need const

		T magnitude_squared();
		T magnitude();

		T dot(const Vector<T, d>&) const;
		const T& operator[](int ) const;
		T& operator[] (int );

		Vector<T, d> operator/ (const T );
		Vector<T, d>& operator/= (const T);
		Vector<T, d> operator- (const Vector&) const;
		Vector<T, d>& operator-= (const Vector&);
		Vector<T, d> operator+ (const Vector&) const;
		Vector<T, d>& operator+= (const Vector&);
		Vector<T, d> operator* (const T) const;
		Vector<T, d>& operator*=(const T);
};
/*template<class T, int d, class... A>
Vector<T, d>::Vector(const A&... args) : arr{ args... } {}*/

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
Vector<T, d> Vector<T, d>::operator- (const Vector& v) const{
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
Vector<T, d> Vector<T, d>::operator* (const T right) const{
	Vector<T, d> newArr;
	for (int i = 0; i < d; i++) {
		newArr[i] = arr[i] * right;
	}
	return newArr;
}



template<class T, int d, class W> //You do not need two template declarations here because this is not a member function
Vector<T, d> operator* (const W left, const Vector<T,d>& right) {
	return right * left;
}

template<class T, int d>
ostream& operator<< (ostream& out, const Vector<T, d>& v) {
	out << "<";
	for (int i = 0; i < d - 1; i++) {
		out << v[i] << ", ";
	}
	out << v[d - 1] << ">";
	return out;
}

#endif



