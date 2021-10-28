#include<iostream>
#include<vector>
#include<cmath>
#include<imgo.h>

using namespace std;

struct task {
    double (*f)(double, int);
    double a, b;
    double x_min;
    double m;

    task(double (*_f)(double, int), double _a, double _b, double _x_min, double _m)
        : f(_f), a(_a), b(_b), x_min(_x_min), m(_m) {}
};

double f1(double x, int j) {
    switch (j) {
        case 1: return exp(-sin(3.0 * x)) - 1.0 / 10.0 * pow(x - 1.0 / 2.0, 2.0) - 1.0;
        case 2: return -13.0 / 6.0 * x + sin(13.0 / 4.0 * (2.0 * x + 5.0)) - 53.0 / 12.0;
    }
    return -1.0;
}

double f4(double x, int j) {
    double sum = 0.0;
    switch (j) {
        case 1:
            for (int i = 1; i <= 5; i++) {
                sum += cos(5.0 / 4.0 * (i + 1) * x + i);
            }
            return 6.0 / 25.0 - sum;
        case 2: return 9.0 / 50.0 - 9.0 / 2.0 * exp(-(x - 1.0 / 10.0)) * 
                sin(2.0 * M_PI * (x - 1.0 / 10.0));
        case 3: return 4.0 * sin(M_PI / 4.0 * x + 1.0 / 20.0) * 
                pow(pow(sin(M_PI / 2.0 * x + 1.0 / 10.0), 3.0) + pow(cos(M_PI / 2.0 * x + 1.0 / 10.0), 3.0), 2.0);
    }
    return -1.0;
}

double f8(double x, int j) {
    double sum = 0.0;
    switch (j) {
        case 1:
            return exp(-sin(4.0 * x)) - 1.0 / 10.0 * pow(x - 1.0 / 2.0, 2) - 1.0;
        case 2: 
            for (int i = 1; i <= 5; i++) {
                sum += cos(5.0 * (i + 1) * (x + 1.0 / 2.0));
            }
            return 3.0 / 10.0 - sum;
        case 3: return (-21.0 / 20.0 * x - 13.0 / 8.0) * sin(63.0 / 10.0 * x + 63.0 / 4.0) + 1.0 / 5.0;
        case 4: return cos(7.0 / 4.0 * x + 241.0 / 40.0) - sin(35.0 / 4.0 * x + 241.0 / 8.0) - 5.0;
    }
    return -1.0;
}

int main() {
    double eps = 0.001;
    double r = 3; // > 1

    double x_min;
    int n;

    vector<task> task_arr = {{f1, -2.5, 1.5, 1.05738, 1},
                             {f4, 0.0, 4.0, 2.45956, 2},
                             {f8, -2.5, 1.5, -1.12724, 3}};

    imgo_method imgo(&f1, 0.0, 0.0, 0.0, eps, r);

    for (int i = 0; i < task_arr.size(); i++) {
        imgo.setFunc(task_arr[i].f);
        imgo.setA(task_arr[i].a);
        imgo.setB(task_arr[i].b);
        imgo.setM(task_arr[i].m);

        x_min = imgo.solve(n);

        cout << "f" << i + 1 << "(x)\n";
        cout << "[a; b] = [" << task_arr[i].a << "; " << task_arr[i].b << "]"<< endl;
        cout << "X_min = " << x_min << endl;
        cout << "X_min_true = " << task_arr[i].x_min << endl;
        cout << "Number of trials = " << n << endl;
        cout << "|Error rate| = " << abs(x_min - task_arr[i].x_min) << endl;
        cout << endl;
    }
    return 0;
}