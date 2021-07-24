#include "Vector.h"
#include <math.h>

template<class T, int n>
struct DVS {
	Vector<T, n> u;
	Vector<T, n> du;

	//dot product => u dot v = du dot v + u dot dv
	//magnitude => w = w/ mag(w) * dw
	//+- Operator
};

template <class T>
struct DS {
	T s;
	T ds;

	DS(T scalar, T deriveS) : s(scalar), ds(deriveS) {};

	DS<T> operator*(const DS<T>& right) {
		DS<T> result(s * right.s, (ds * right.s) + (s * right.ds));
		return result;
	}
	DS<T> operator-(const DS<T>& right) {
		DS<T> result(s - right.s, ds - right.ds);
		return result;
	}
	DS<T> operator+(const DS<T>& right) {
		DS<T> result(s + right.s, ds + right.ds);
		return result;
	}
	DS<T> operator/(const DS<T>& right) {
		DS<T> result(s / right.s, (ds * right.s - s * right.ds) / (right.s * right.s));
		return result;
	}
};

template <class T>
DS<T> sin(const DS<T>& curr){
	DS<T> result(sin(curr.s), curr.ds* cos(curr.s));
	return result;
}

template <class T>
DS<T> cos(const DS<T>& curr){
	DS<T> result(cos(curr.s), -1 * curr.ds * sin(curr.s));
	return result;
}

template <class T>
DS<T> exp(const DS<T>& curr) {
	DS<T> result(exp(curr.s), curr.ds * exp(curr.s));
	return result;
}

template <class T>
DS<T> log(const DS<T>& curr) {
	DS<T> result(log(curr.s), curr.ds / curr.s);
	return result;
}