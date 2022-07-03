#include<iostream>
#include<cmath>
#include<gsa.h>

using namespace std;

double f1(double x) {
    return -4.0 * x + 1.0;
}

double f2(double x) {
    return 5.0 * x * x + 3.0 * x - 1.0;
}

double f3(double x) {
    return x * sin(x);
}

double f4(double x) {
    if (x != 0) {
        return x * sin(1 / x);
    } else {
        return 0.0;
    }
}

int main() {
    double a = 3.0, b = 4.0;
    double eps = 0.0001;
    double r = 2.0; // > 1

    double x_min;
    double x_opt;
    int n;

    // линейная функция
    gsa_method gsa(&f1, a, b, eps, r);
    x_min = gsa.solve(n);
    x_opt = 4.0;

    cout << "f1(x) = -4.0 * x + 1.0\n";
    cout << "[a; b] = [" << a << "; " << b << "]"<< endl;
    cout << "x_min = " << x_min << endl;
    cout << "Number of trials = " << n << endl;
    cout << "|f1(x_min) - f1(x_opt)| = " << abs(f1(x_min) - f1(x_opt)) << endl;
    cout << endl;

    // квадратичная функция
    a = -2.0; b = 2.0;
    x_opt = -0.3;
    gsa.setFunc(&f2);
    gsa.setA(a);
    gsa.setB(b);
    x_min = gsa.solve(n);

    cout << "f2(x) = 5.0 * x * x + 3.0 * x - 1.0\n";
    cout << "[a; b] = [" << a << "; " << b << "]"<< endl;
    cout << "x_min = " << x_min << endl;
    cout << "Number of trials = " << n << endl;
    cout << "|f2(x_min) - f2(x_opt)| = " << abs(f2(x_min) - f2(x_opt)) << endl;
    cout << endl;

    // 1 функция с sin
    a = 0.0; b = 20.0;
    x_opt = 17.336;
    gsa.setFunc(&f3);
    gsa.setA(a);
    gsa.setB(b);
    x_min = gsa.solve(n);

    cout << "f3(x) = x * sin(x)\n";
    cout << "[a; b] = [" << a << "; " << b << "]"<< endl;
    cout << "x_min = " << x_min << endl;
    cout << "Number of trials = " << n << endl;
    cout << "|f3(x_min) - f3(x_opt)| = " << abs(f3(x_min) - f3(x_opt)) << endl;
    cout << endl;

    // 2 функция с sin
    a = -0.4; b = 0.4;
    x_opt = -0.2225;
    gsa.setFunc(&f4);
    gsa.setA(a);
    gsa.setB(b);
    x_min = gsa.solve(n);

    cout << "f4(x) = x * sin(1 / x)\n";
    cout << "[a; b] = [" << a << "; " << b << "]"<< endl;
    cout << "x_min = " << x_min << endl;
    cout << "Number of trials = " << n << endl;
    cout << "|f4(x_min) - f4(x_opt)| = " << abs(f4(x_min) - f4(x_opt)) << endl;
    cout << endl;

    #if defined( _MSC_VER )
        cin.get();
    #endif
    return 0;
}