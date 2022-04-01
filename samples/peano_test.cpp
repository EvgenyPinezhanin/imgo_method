#include <iostream>
#include <cmath>
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

int main() {
    double a1, b1, a2, b2, x_opt, y_opt;
    int m, number_trials;
    vector<double> X(2);

    int n = 2, den = 10, key = 1;
    double eps = 0.0001, r = 2.0, d = 0.0;

    a1 = 0, b1 = 4.0;
    a2 = -1.0, b2 = 3.0;
    x_opt = 0.942, y_opt = 0.944;
    m = 3;

    imgo_method imgo(f, n, m, vector<double>{a1, a2}, vector<double>{b1, b2}, eps, r, d, den, key);
    imgo.solve(number_trials, X);

    std::cout << "Parameters for constructing the Peano curve:" << std::endl;
    std::cout << "n = " << n << " m = " << den << " key = " << key << std::endl;
    std::cout << "Trials result:" << std::endl;
    std::cout << "Number of trials = " << number_trials << std::endl;
    std::cout << "x*_min = " << x_opt << " y*_min = " << y_opt << std::endl;
    std::cout << "x_min = " << X[0] << " y_min = " << X[1] << std::endl;
    std::cout << "||X* - X|| = " << sqrt((x_opt - X[0]) * (x_opt - X[0]) + (y_opt - X[1]) * (y_opt - X[1])) << std::endl;

    return 0;
}