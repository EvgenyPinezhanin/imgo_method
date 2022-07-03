#include<iostream>
#include<iomanip>
#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include<math.h>
#else
    #include<cmath>
#endif
#include<gsa.h>

using namespace std;

struct task {
    double (*f)(double);
    double a, b;
    double x_min;

    task(double (*_f)(double), double _a, double _b, double _x_min)
        : f(_f), a(_a), b(_b), x_min(_x_min) {}
};

double f1(double x) {
    return -(-1.0 / 6.0 * pow(x, 6) + 52.0 / 25.0 * pow(x, 5) - 39.0 / 80.0 * pow(x, 4) - 71.0 / 10.0 * pow(x, 3) + 79.0 / 20.0 * x * x + x - 1.0 / 10.0);
}

double f2(double x) {
    return -(-sin(x) - sin(10.0 / 3.0 * x));
}

double f3(double x) {
    double sum = 0;
    for (int i = 1; i <= 5; i++) {
        sum += i * sin((i + 1.0) * x  + i);
    }
    return -sum;
}

double f4(double x) {
    return -(16.0 * x * x - 24.0 * x + 5.0) * exp(-x);
}

double f5(double x) {
    return -(-3.0 * x + 1.4) * sin(18.0 * x);
}

double f6(double x) {
    return -((x + sin(x)) * exp(-x * x));
}

double f7(double x) {
    return -(-sin(x) -sin(10.0 / 3.0 * x) - log(x) + 0.84 * x - 3);
}

double f8(double x) {
    double sum = 0;
    for (int i = 1; i <= 5; i++) {
        sum += i * cos((i + 1.0) * x  + i);
    }
    return -sum;
}

double f9(double x) {
    return -(-sin(x) - sin(2.0 / 3.0 * x));
}

double f10(double x) {
    return -(x * sin(x));
}

double f11(double x) {
    return -(-2.0 * cos(x) - cos(2.0 * x));
}

double f12(double x) {
    return -(-pow(sin(x), 3) - pow(cos(x), 3));
}

double f13(double x) {
    if (x * x - 1 < 0) {
        return -(pow(x, 2.0 / 3.0) + pow(-(x * x - 1), 1.0 / 3.0));
    }
    return -(pow(x, 2.0 / 3.0) - pow(x * x - 1, 1.0 / 3.0));
}

double f14(double x) {
    return -(exp(-x) * sin(2 * M_PI * x));
}

double f15(double x) {
    return -((-x * x + 5.0 * x -6.0) / (x * x + 1));
}

double f16(double x) {
    return -(-2.0 * (x - 3) * (x - 3) - exp(- x * x / 2));
}

double f17(double x) {
    return -(-pow(x, 6) + 15.0 * pow(x, 4) - 27.0 * x * x - 250.0);
}

double f18(double x) {
    if (x <= 3.0) {
        return (x - 2.0) * (x - 2.0);
    } 
    return -(-2.0 * log(x - 2.0) - 1);
}

double f19(double x) {
    return -(x - sin(3.0 * x) + 1);
}

double f20(double x) {
    return -(x - sin(x)) * exp(- x * x);
}

int main() {
    double eps = 0.001;
    double r = 2.0; // > 1

    double x_min;
    int n;

    vector<task> task_arr = { {f1, -1.5, 11.0, 10.0},
                              {f2, 2.7, 7.5, 5.145735},
                              {f3, -10.0, 10.0, 5.791785},
                              {f4, 1.9, 3.9, 2.868},
                              {f5, 0.0, 1.2, 0.96609},
                              {f6, -10.0, 10.0, 0.67956},
                              {f7, 2.7, 7.5, 5.19978},
                              {f8, -10.0, 10.0, -7.0835},
                              {f9, 3.1, 20.4, 17.039},
                              {f10, 0.0, 10.0, 7.9787},
                              {f11, -1.57, 6.28, 2.094},
                              {f12, 0.0, 6.28, 3.142},
                              {f13, 0.001, 0.99, 0.7071},
                              {f14, 0.0, 4.0, 0.224885},
                              {f15, -5.0, 5.0, 2.4142},
                              {f16, -3.0, 3.0, 3.0},
                              {f17, -4.0, 4.0, -3.0},
                              {f18, 0.0, 6.0, 2.0},
                              {f19, 0.0, 6.5, 5.87287},
                              {f20, -10.0, 10.0, 1.195137} };

    gsa_method gsa(&f1, 0.0, 0.0, eps, r);

    for (int i = 0; i < task_arr.size(); i++) {
        gsa.setFunc(task_arr[i].f);
        gsa.setA(task_arr[i].a);
        gsa.setB(task_arr[i].b);
        x_min = gsa.solve(n);

        cout << "f" << i + 1 << "(x)\n";
        cout << "[a; b] = [" << task_arr[i].a << "; " << task_arr[i].b << "]" << endl;
        cout << "x_min = " << setprecision(6) << x_min << endl;
        cout << "Number of trials = " << n << endl;
        cout << "|f" << i + 1 << "(x_min) - f" << i + 1 << "(x_opt)| = " << 
                abs(task_arr[i].f(x_min) - task_arr[i].f(task_arr[i].x_min)) << endl;
        cout << endl;
    }

    #if defined( _MSC_VER )
        cin.get();
    #endif
    return 0;
}