#include <iostream>
#include "Derive.h"

using namespace std;

int main() {
	double x = 5;
	double y = 2 * x;
	DS<double> test0(y, 2); //In terms of x
	DS<double> test1(x, 0);
	cout << (test0 / test1).s << " " << (test0 / test1).ds << endl;

}