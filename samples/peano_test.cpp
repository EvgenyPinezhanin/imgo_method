#include <iostream>
#include <cmath>
#include <Grishagin/GrishaginProblemFamily.hpp>
#include <imgo.h>

double f_test(vector<double> x, int j) {
    switch (j) {
        case 1: return 0.01 * (pow((x[0] - 2.2), 2.0) + pow((x[1] - 1.2), 2.0) - 2.25);
        case 2: return 100.0 * (1.0 - pow((x[0] - 2.0), 2.0) / 1.44 - pow(0.5 * x[1], 2.0));
        case 3: return 10.0 * (x[1] - 1.5 - 1.5 * sin(6.283 * (x[0] - 1.75)));
        case 4: return -1.5 * x[0] * x[0] * exp(1.0 - x[0] * x[0] - 20.25 * pow((x[0] - x[1]), 2.0)) - 
                       pow(0.5 * (x[0] - 1.0) * (x[1] - 1.0), 4.0) * exp(2.0 - pow(0.5 * (x[0] - 1.0), 4.0) - 
                       pow(x[1] - 1.0, 4.0));
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f_SPPR1(vector<double> x, int j) {
    switch (j) {
        case 1: return 1.0 - 2.0 * x[0] - x[0] * x[1];
        case 2: return 2.0 * x[0] * x[0] + x[1];
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f_SPPR2(vector<double> x, int j) {
    switch (j) {
        case 1: return sin(x[0]) - x[1];
        case 2: return cos(x[0]) + x[1];
        default: return numeric_limits<double>::quiet_NaN();
    }
}

TGrishaginProblemFamily grishaginProblem;
int current_func;
double f_grishagin(vector<double> x, int j) {
    switch (j) {
        case 1: return grishaginProblem[current_func]->ComputeFunction(x);
        default: return numeric_limits<double>::quiet_NaN();
    }
}

int main() {
    vector<double> A(2), B(2), X_opt(2), X(2);
    int m, number_trials;

    int n = 2, den = 10, key = 1;
    double eps = 0.001, r = 2.0, d = 0.0;

    A[0] = 0, B[0] = 4.0;
    A[1] = -1.0, B[1] = 3.0;
    X_opt[0] = 0.942, X_opt[1] = 0.944;
    m = 3;

    imgo_method imgo(f_test, n, m, A, B, r, d, eps, 0, den, key);
    imgo.solve(number_trials, X);

    std::cout << "Test_function" << std::endl;
    std::cout << "Parameters for constructing the Peano curve:" << std::endl;
    std::cout << "n = " << n << " m = " << den << " key = " << key << std::endl;
    std::cout << "Trials result:" << std::endl;
    std::cout << "Number of trials = " << number_trials << std::endl;
    std::cout << "x*_min = " << X_opt[0] << " y*_min = " << X_opt[1] << std::endl;
    std::cout << "x_min = " << X[0] << " y_min = " << X[1] << std::endl;
    std::cout << "||X* - X|| = " << sqrt((X_opt[0] - X[0]) * (X_opt[0] - X[0]) + (X_opt[1] - X[1]) * (X_opt[1] - X[1])) << std::endl;
    std::cout << std::endl;

/*     A[0] = 0.5, B[0] = 3.0;
    A[1] = -2.0, B[1] = 2.0;
    X_opt[0] = 0.629960524947, X_opt[1] = -0.412598948032;
    m = 1;

    imgo.setA(A);
    imgo.setB(B);
    imgo.setM(m);
    imgo.setFunc(f_SPPR1);

    imgo.solve(number_trials, X);

    std::cout << "Test_function(SPPR-1)" << std::endl;
    std::cout << "Parameters for constructing the Peano curve:" << std::endl;
    std::cout << "n = " << n << " m = " << den << " key = " << key << std::endl;
    std::cout << "Trials result:" << std::endl;
    std::cout << "Number of trials = " << number_trials << std::endl;
    std::cout << "x*_min = " << X_opt[0] << " y*_min = " << X_opt[1] << std::endl;
    std::cout << "x_min = " << X[0] << " y_min = " << X[1] << std::endl;
    std::cout << "||X* - X|| = " << sqrt((X_opt[0] - X[0]) * (X_opt[0] - X[0]) + (X_opt[1] - X[1]) * (X_opt[1] - X[1])) << std::endl;
    std::cout << std::endl;

    A[0] = 0.0, B[0] = M_PI;
    A[1] = 0.0, B[1] = 2.0;
    X_opt[0] = M_PI, X_opt[1] = 0.0;
    m = 1;

    imgo.setA(A);
    imgo.setB(B);
    imgo.setM(m);
    imgo.setFunc(f_SPPR2);

    imgo.solve(number_trials, X);

    std::cout << "Test_function(SPPR-2)" << std::endl;
    std::cout << "Parameters for constructing the Peano curve:" << std::endl;
    std::cout << "n = " << n << " m = " << den << " key = " << key << std::endl;
    std::cout << "Trials result:" << std::endl;
    std::cout << "Number of trials = " << number_trials << std::endl;
    std::cout << "x*_min = " << X_opt[0] << " y*_min = " << X_opt[1] << std::endl;
    std::cout << "x_min = " << X[0] << " y_min = " << X[1] << std::endl;
    std::cout << "||X* - X|| = " << sqrt((X_opt[0] - X[0]) * (X_opt[0] - X[0]) + (X_opt[1] - X[1]) * (X_opt[1] - X[1])) << std::endl;
    std::cout << std::endl; */

    return 0;
}
