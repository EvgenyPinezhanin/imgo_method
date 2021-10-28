#include<iostream>
#include<vector>
#include<cmath>
#include<imgo.h>

using namespace std;

struct task {
    double (*f)(double, int);
    double a, b;
    double x_min;

    task(double (*_f)(double, int), double _a, double _b, double _x_min)
        : f(_f), a(_a), b(_b), x_min(_x_min) {}
};

double f1(double x) {
    return ;
}

int main() {
    double eps = 0.001;
    double r = 3; // > 1

    double x_min;
    int n;

    int numFunc = 20;
    double (*f_arr[])(double) = {f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14,
                                 f15, f16, f17, f18, f19, f20};
    double intervals[numFunc][2] = {{-1.5, 11}, {2.7, 7.5}, {-10, 10}, {1.9, 3.9}, {0, 1.2}, 
                                    {-10.0, 10.0}, {2.7, 7.5}, {-10.0, 10.0}, {3.1, 20.4},
                                    {0, 10}, {-1.57, 6.28}, {0, 6.28}, {0.001, 0.99}, {0, 4.0}, {-5.0, 5.0},
                                    {-3.0, 3.0}, {-4.0, 4.0}, {0, 6.0}, {0, 6.5}, {-10.0, 10.0}};
    double optimum[] = {10, 5.145735, -6.775, 2.868, 0.96609, 0.67956, 5.19978, -7.084, 17.039,
                        7.9787, 2.094, 3.142, 0.7071, 0.224885, 2.4142, 1.5907, -3.0, 2.0, 5.87287, 1.195137};

    imgo_method gsa(&f1, 0.0, 0.0, eps, r);

    for (int i = 0; i < numFunc; i++) {
        gsa.setFunc(*f_arr[i]);
        gsa.setA(intervals[i][0]);
        gsa.setB(intervals[i][1]);
        x_min = gsa.solve(n);

        cout << "f" << i + 1 << "(x)\n";
        cout << "[a; b] = [" << intervals[i][0] << "; " << intervals[i][1] << "]"<< endl;
        cout << "X_min = " << x_min << endl;
        cout << "Number of trials = " << n << endl;
        cout << "|Error rate| = " << abs(x_min - optimum[i]) << endl;
        cout << endl;
    }
    return 0;
}