#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <imgo.h>

//#define FUNC_1
#define FUNC_2
//#define FUNC_3
//#define FUNC_4

#if defined(FUNC_1)
double f(vector<double> x, int j) {
    switch (j) {
        case 1: return (x[0] - 1.0) * (x[0] - 1.0) / 5.0 + (x[1] - 1.0) * (x[1] - 1.0) / 5.0;
    }
    return -1;
}
#elif defined(FUNC_2)
const double k = 0.3;
double f(vector<double> x, int j) {
    switch (j) {
        case 1: return k - x[0] - x[1];
        case 2: return x[0] * x[0] / 5.0 + x[1] * x[1] / 5.0;
    }
    return -1;
}
#elif defined(FUNC_3)
double f(vector<double> x, int j) {
    switch (j) {
        case 1: return 1.0 - x[0] - x[1];
    }
    return -1;
}
#endif

int main() {
    std::ofstream ofstr("trial_points.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";

#if defined(FUNC_1)
    double a1 = -4.0, b1 = 4.0;
    double a2 = -4.0, b2 = 4.0;
    double x1_opt = 1.0, x2_opt = 1.0;
    int m = 0;
#elif defined(FUNC_2)
    double a1 = -0.5, b1 = 0.5;
    double a2 = -0.5, b2 = 0.5;
    double x1_opt = k / 2.0, x2_opt = k / 2.0;
    int m = 1;
#elif defined(FUNC_3)
    double a1 = -0.5, b1 = 0.5;
    double a2 = -0.5, b2 = 0.5;
    double x1_opt = 4.0, x2_opt = 4.0;
    int m = 0;
#endif

    vector<double> X(2);
    std::vector<vector<double>> trial_vec;
    int number_trials;
    int n = 2, den = 8, key = 2;
    double eps = 0.0001, r = 2.5, d = 0.01;

    imgo_method imgo(f, n, m, vector<double>{a1, a2}, vector<double>{b1, b2}, eps, r, d, den, key);
    imgo.solve(number_trials, X);

    std::cout << "Parameters for constructing the Peano curve:" << std::endl;
    std::cout << "n = " << n << " m = " << den << " key = " << key << std::endl;
    std::cout << "Trials result:" << std::endl;
    std::cout << "Number of trials = " << number_trials << std::endl;
    std::cout << "x*_min = " << x1_opt << " y*_min = " << x2_opt << std::endl;
    std::cout << "x_min = " << X[0] << " y_min = " << X[1] << std::endl;

    // Построение графика
    ofstr << X[0] << " " << X[1] << " " << f(X, m + 1) << endl;
    ofstr << endl << endl;
    ofstr << x1_opt << " " << x2_opt << " " << f(vector<double>{x1_opt, x2_opt}, m + 1) << endl;
    ofstr << endl << endl;
    imgo.getPoints(trial_vec);
    for (int i = 0; i < trial_vec.size(); i++) {
        ofstr << trial_vec[i][0] << " " << trial_vec[i][1] << " " << f(trial_vec[i], m + 1) << std::endl;
    }
    ofstr.close();

#if defined(__linux__)
    int error;
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x chart.gp");
    if (error != 0) {
        std::cerr << "Error chmod" << std::endl;
    }

#if defined(FUNC_1)
    error = system("gnuplot -p -c chart.gp func_1");
#elif defined(FUNC_2)
    char str[100];
    sprintf(str, "gnuplot -p -c chart.gp func_2 %g", k);
    error = system(str);
#elif defined(FUNC_3)
    error = system("gnuplot -p -c chart.gp func_3");
#endif
#endif

    if (error != 0) {
        std::cerr << "Error gnuplot" << std::endl;
    }
    return 0;
}