#include <iostream>
#include <cmath>
#include <Grishagin/GrishaginProblemFamily.hpp>
#include <Grishagin/GrishaginConstrainedProblemFamily.hpp>
#include <imgo.h>

double f(vector<double> x, int j) {
    switch (j) {
        case 1: return 0.01 * (pow((x[0] - 2.2), 2.0) + pow((x[1] - 1.2), 2.0) - 2.25);
        case 2: return 100.0 * (1.0 - pow((x[0] - 2.0), 2.0) / 1.44 - pow(0.5 * x[1], 2.0));
        case 3: return 10.0 * (x[1] - 1.5 - 1.5 * sin(6.283 * (x[0] - 1.75)));
        case 4: return -1.5 * x[0] * x[0] * exp(1.0 - x[0] * x[0] - 20.25 * pow((x[0] - x[1]), 2.0)) - 
                       pow(0.5 * (x[0] - 1.0) * (x[1] - 1.0), 4.0) * exp(2.0 - pow(0.5 * (x[0] - 1.0), 4.0) - 
                       pow(x[1] - 1.0, 4.0));
    }
    return -1;
}

TGrishaginProblemFamily grishaginProblem;
int current_func;
double f_grishagin(vector<double> x, int j) {
    return grishaginProblem[current_func]->ComputeFunction(x);
}

GrishaginConstrainedProblem grishaginConstrainedProblem(cptInFeasibleDomain, 0.3, 0, 1);
double f_grishaginConstrained(vector<double> x, int j) {
    if (j > grishaginConstrainedProblem.GetConstraintsNumber()) {
        return grishaginConstrainedProblem.ComputeFunction(x);
    } else {
        return grishaginConstrainedProblem.ComputeConstraint(j, x);
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

    imgo_method imgo(f, n, m, A, B, eps, r, d, den, key);
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

    imgo.setFunc(f_grishagin);
    n = grishaginProblem[current_func]->GetDimension();
    imgo.setN(n);
    grishaginProblem[current_func]->GetBounds(A, B);
    imgo.setA(A);
    imgo.setB(B);
    n = grishaginProblem[current_func]->GetDimension();
    imgo.setN(n);
    m = 0;
    imgo.setM(m);
    r = 2.9;
    imgo.setR(r);
    for (current_func = 0; current_func < 100; current_func++) {
        imgo.solve(number_trials, X);
        X_opt = grishaginProblem[current_func]->GetOptimumPoint();

        std::cout << "grishagin_function " << current_func + 1 << std::endl;
        std::cout << "Parameters for constructing the Peano curve:" << std::endl;
        std::cout << "n = " << n << " m = " << den << " key = " << key << std::endl;
        std::cout << "Trials result:" << std::endl;
        std::cout << "Number of trials = " << number_trials << std::endl;
        std::cout << "x*_min = " << X_opt[0] << " y*_min = " << X_opt[1] << std::endl;
        std::cout << "x_min = " << X[0] << " y_min = " << X[1] << std::endl;
        std::cout << "||X* - X|| = " << sqrt((X_opt[0] - X[0]) * (X_opt[0] - X[0]) + (X_opt[1] - X[1]) * (X_opt[1] - X[1])) << std::endl;
        std::cout << std::endl;
    }

   //grishaginConstrainedProblem.GetBounds(A, B);
   //imgo.setA(A);
   //imgo.setB(B);
   //X_opt = grishaginConstrainedProblem.GetOptimumPoint();
   //n = grishaginConstrainedProblem.GetDimension();
   //imgo.setN(n);
   //m = grishaginConstrainedProblem.GetConstraintsNumber();
   //imgo.setM(m);
   //imgo.setFunc(f_grishaginConstrained);

   //imgo.solve(number_trials, X);

   //std::cout << "grishagin_constrained_function" << std::endl;
   //std::cout << "Parameters for constructing the Peano curve:" << std::endl;
   //std::cout << "n = " << n << " m = " << den << " key = " << key << std::endl;
   //std::cout << "Trials result:" << std::endl;
   //std::cout << "Number of trials = " << number_trials << std::endl;
   //std::cout << "x*_min = " << X_opt[0] << " y*_min = " << X_opt[1] << std::endl;
   //std::cout << "x_min = " << X[0] << " y_min = " << X[1] << std::endl;
   //std::cout << "||X* - X|| = " << sqrt((X_opt[0] - X[0]) * (X_opt[0] - X[0]) + (X_opt[1] - X[1]) * (X_opt[1] - X[1])) << std::endl;

    return 0;
}