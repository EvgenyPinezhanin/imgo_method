#include<iostream>
#include<cmath>
#include<imgo.h>

using namespace std;

double f1(double x, int j) {
    switch(j) {
        case 1: return sin(x);
        case 2: return -2.0 * x + 3.0;
    }
    return -1.0;
}

double f2(double x, int j) {
    switch (j) {
        case 1: return x * x - 0.05;
        case 2: return -x + 0.1;
        case 3: return 5.0 * x * x + 3.0 * x - 1.0;
    }
    return -1.0;
}

/* double f3(double x) {
    return x * sin(x);
}

double f4(double x) {
    if (x != 0) {
        return x * sin(1 / x);
    } else {
        return 0.0;
    }
} */

int main() {
    double a = 2.0, b = 8.0;
    double eps = 0.001;
    double r = 3;
    int m;

    double x_min;
    double x_min_true;
    int n;

    m = 1;
    imgo_method imgo(&f1, m, a, b, eps, r);
    x_min = imgo.solve(n);
    x_min_true = 2.0 * M_PI;

    cout << "f1(x) = -2.0 * x + 3.0\n";
    cout << "g1(x) = sin(x)\n";
    cout << "[a; b] = [" << a << "; " << b << "]"<< endl;
    cout << "X_min_true = " << x_min_true << endl;
    cout << "X_min = " << x_min << endl;
    cout << "Number of trials = " << n << endl;
    cout << "|Error rate| = " << abs(x_min - x_min_true) << endl;
    cout << endl;

    m = 2;
    a = -2.0; b = 2.0;
    imgo.setFunc(&f2);
    imgo.setA(a);
    imgo.setB(b);
    imgo.setM(m);
    x_min = imgo.solve(n);
    x_min_true = 0.1;

    cout << "f2(x) = 5.0 * x * x + 3.0 * x - 1.0\n";
    cout << "g1(x) = x * x - 0.05\n";
    cout << "g2(x) = -x + 0.1\n";
    cout << "[a; b] = [" << a << "; " << b << "]"<< endl;
    cout << "X_min = " << x_min << endl;
    cout << "Number of trials = " << n << endl;
    cout << "|Error rate| = " << abs(x_min - x_min_true) << endl;
    cout << endl;
/*
    // 1 функция с sin
    a = 0.0; b = 20.0;
    x_opt = 17.336;
    gsa.setFunc(&f3);
    gsa.setA(a);
    gsa.setB(b);
    x_min = gsa.solve(n);

    cout << "f3(x) = x * sin(x)\n";
    cout << "[a; b] = [" << a << "; " << b << "]"<< endl;
    cout << "X_min = " << x_min << endl;
    cout << "Number of trials = " << n << endl;
    cout << "|Error rate| = " << abs(x_min - x_opt) << endl;
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
    cout << "X_min = " << x_min << endl;
    cout << "Number of trials = " << n << endl;
    cout << "|Error rate| = " << abs(x_min - x_opt) << endl;
    cout << endl; */
    return 0;
}