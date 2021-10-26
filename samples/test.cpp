#include<iostream>
#include<cmath>
#include<imgo.h>

using namespace std;

double f1(double x) {
    return -(-1.0 / 6.0 * pow(x, 6) + 52.0 / 25.0 * pow(x, 5) - 39.0 / 80.0 * pow(x, 4) - 71.0 / 10.0 * pow(x, 3) + 79.0 / 20.0 * x * x + x - 1.0 / 10.0);
}

double f2(double x) {
    return -(-sin(x) - sin(10.0 / 3.0 * x));
}

double f3(double x) {
    double sum = 0;
    for (int i = 1; i <= 5; i++) {
        sum += i * sin((i + 1) * x  + i);
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
        sum += i * cos((i + 1) * x  + i);
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
    //double eps = 0.5;
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
    // optimum f3: -0.49139, 5.791785
    // optimum f8: -0.8003
    // optimum f11: 4.189
    // optimum f12: 4.712
    // optimum f17: 3.0

    gsa_method gsa(&f1, 0.0, 0.0, eps, r);

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