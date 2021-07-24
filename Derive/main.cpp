#include <iostream>
#include <cmath>
#include <random>  
#include <iomanip> 

using namespace std;

void compute(double, double&, double&, double&);
void print(const double, const double&, double&, double&);
double check(const double, const double, const int, const int);
double precise(const int);

int main() {
    double f = 0.0;
    double df = 0.0;
    double ddf = 0.0;

    std::random_device rd; //random device engine, usually based on /dev/random on UNIX-like systems
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> di(0.0, 1.0);
    double x = di(generator);

    cout << "\tTesting cos(x)^(arctan(x)/x)" << endl;
    cout << "Type 1 => (f(x+h) - f(x))/h" << "\tType 2 => (f(x+h) - f(x-h))/(2*h)" << endl;
    for (int i = 0; i < 2; i++) {
        cout << "------------------------------------------------" << endl;
        compute(x, f, df, ddf);
        print(x, f, df, ddf);
        x = di(generator);
    }
}


void compute(double x, double& f, double& df, double& ddf) {

    double checker = 0;
    double dChecker = 0;
    double h = 0;

    double a = cos(x);
    double b = atan(x);
    double c = x;
    double d = c * c;
    double e = d + 1;

    double one = b / c;
    double two = log(2 + a);

    double da = sin(x);
    double db = 1 / e;
    double dd = 2 * c;

    double dOne = (c - b * e) / (d * e);
    double dTwo = -da / (2 + a);

    double ddOne = ((-dd * b) * (d * e) - (2 * c + 4 * c * d) * (c - b * e)) / (e * e * d * d);
    double ddTwo = (-2 * a - 1) / ((2 + a) * (2 + a));

    f = pow(2 + a, b / c); // Change to e^(b*lna) 
    df = f * (dOne * two + dTwo * one); // e^udu where u is b*lna
    ddf = f * (ddOne * two + ddTwo * one + 2 * dTwo * dOne) + df * (dOne * two + dTwo * one);

}

double check(const double x, const double org, const int i, const int type) {
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> di(0.0, 1.0);
    double h = 0.0;
    h = di(generator) / 1000000;

    double expf, expf1 = 0;
    double expDf, expDf1 = 0;
    double expDdf, expDdf1 = 0;

    if (i == 1 && type == 1) {
        compute(x + h, expf, expDf, expDdf);
        return (expf - org) / h;
    }
    else if (i == 1 && type == 2) {
        compute(x + h, expf, expDf, expDdf);
        compute(x - h, expf1, expDf1, expDdf1);
        return (expf - expf1) / (2 * h);
    }
    else if (i == 2 && type == 1) {
        compute(x + h, expf, expDf, expDdf);
        return (expDf - org) / h;
    }
    else if (i == 2 && type == 2) {
        compute(x + h, expf, expDf, expDdf);
        compute(x - h, expf1, expDf1, expDdf1);
        return (expDf - expDf1) / (2 * h);
    }
}

void print(const double x, const double& f, double& df, double& ddf) {
    std::cout << setprecision(4);

    double expDf = check(x, f, 1, 1);
    double expDf1 = check(x, f, 1, 2);

    double expDdf = check(x, df, 2, 1);
    double expDdf1 = check(x, df, 2, 2);

    double check = precise(df);
    double checkTwo = precise(ddf);
    double dif = 0.0;

    std::cout << "f(" << x << ") = " << f << std::endl;

    std::cout << "\tActual df(" << x << ") = " << df << std::endl;
    std::cout << "\t   Expected[Type 1] df(" << x << ") = " << expDf << std::endl;
    dif = abs(df - expDf);
    if (dif < check) {
        std::cout << "\t   PASS Difference : " << dif << std::endl;
    }
    else {
        std::cout << "\t   FAIL Difference : " << dif << std::endl;
    }
    std::cout << "\t   Expected[Type 2] df(" << x << ") = " << expDf1 << std::endl;
    dif = abs(df - expDf1);
    if (dif < check) {
        std::cout << "\t   PASS Difference : " << dif << std::endl;
    }
    else {
        std::cout << "\t   FAIL Difference : " << dif << std::endl;
    }



    std::cout << "\tActual ddf(" << x << ") = " << ddf << std::endl;
    std::cout << "\t   Expected[Type 1] ddf(" << x << ") = " << expDdf << std::endl;

    dif = abs(ddf - expDdf);
    if (dif < checkTwo) {
        std::cout << "\t   PASS Difference : " << dif << std::endl;
    }
    else {
        std::cout << "\t   FAIL Difference : " << dif << std::endl;
    }
    std::cout << "\t   Expected[Type 2] ddf(" << x << ") = " << expDdf1 << std::endl;
    dif = abs(ddf - expDdf1);
    if (dif < checkTwo) {
        std::cout << "\t   PASS Difference : " << dif << std::endl;
    }
    else {
        std::cout << "\t   FAIL Difference : " << dif << std::endl;
    }
}

double precise(const int c) {
    int size = to_string(c).length();

    if (size == 6) {
        return 10;
    }
    else if (size == 5) {
        return 1;
    }
    else if (size == 4) {
        return 0.1;
    }
    else if (size == 3) {
        return 0.01;
    }
    else {
        return 0.001;
    }
}