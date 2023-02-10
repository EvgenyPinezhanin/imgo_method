#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <mggsa.h>
#include <map.h>

using namespace std;

double f(vector<double> x, int j) {
    switch (j) {
        case 1: return x[0] * x[0] + x[1] * x[1] - cos(18.0 * x[0]) - cos(18.0 * x[1]);
        default: return numeric_limits<double>::quiet_NaN();
    }
}

int main() {
    ofstream ofstr("output_data/evolvent_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_points("output_data/evolvent_test_points.txt");
    if (!ofstr_points.is_open()) cerr << "File opening error\n";

    vector<double> A{-1.0 / 2.0, -1.0 / 2.0}, B{1.0, 1.0}, X_opt{0.0, 0.0};
    double eps = 0.01, r = 2.0, d = 0.0;
    int constr = 0, count_trials, Nmax = 1000;
    int m = 10, n = 2, key = 2, incr = 10;
    vector<double> X(n); 
    Stop stop = Stop::ACCURACY;

    double k = (key != 3) ? 1.0 / (pow(2.0, n * m) - 1.0) :
                            1.0 / (pow(2.0, m * n) * (pow(2.0, n) - 1.0)) + 0.0000000001;
    if (key == 3) m++;
    for (double i = 0.0; i <= 1.0; i += k) {
        mapd(i, m, X.data(), n, key);
        ofstr << X[0] * (B[0] - A[0]) + (A[0] + B[0]) / 2.0 << " " 
              << X[1] * (B[1] - A[1]) + (A[1] + B[1]) / 2.0 << endl;
    }
    mapd(1.0, m, X.data(), n, key);
    ofstr << X[0] * (B[0] - A[0]) + (A[0] + B[0]) / 2.0 << " "
          << X[1] * (B[1] - A[1]) + (A[1] + B[1]) / 2.0 << endl;
    if (key == 3) m--;
    ofstr.close();

    mggsa_method mggsa(f, n, constr, A, B, r, d, m, key, eps, Nmax, incr);

    vector<double> mu;
    vector<vector<double>> points;

    mggsa.solve(count_trials, X, stop);
    mggsa.getMu(mu);

    cout << "Function: " << "x^2 + y^2 - cos(18.0 * x) - cos(18.0 * y)" << endl;
    cout << "Dimension = " << n << endl;
    cout << "Number of constrained = " << constr << endl;
    cout << "[A; B] = [(" << A[0] << ", " << A[1] << "); (" << 
                             B[0] << ", " << B[1] << ")]"<< endl;
    cout << "X* = (" << X_opt[0] << ", " << X_opt[1] << ")" << endl;
    cout << "f(X*) = " << f(X_opt, constr + 1) << endl;
    cout << "Parameters for method:" << endl;
    cout << "eps = " << eps << " r = " << r << " d = " << d << endl;
    cout << "Parameters for constructing the Peano curve:" << endl;
    cout << "m = " << m << " key = " << key << " incr = " << incr << endl;
    cout << "Trials result:" << endl;
    cout << "Number of trials = " << count_trials << endl;
    cout << "Number of points = " << mggsa.getCountPoints() << endl;
    cout << "Estimation of the Lipschitz constant = " << mu[0] << endl;
    cout << "X = (" << X[0] << ", " << X[1] << ")" << endl;
    cout << "f(X) = " << f(X, constr + 1) << endl;
    cout << "||X* - X|| = " << sqrt((X_opt[0] - X[0]) * (X_opt[0] - X[0]) + 
                                    (X_opt[1] - X[1]) * (X_opt[1] - X[1])) << endl;
    cout << "|f(X*) - f(X)| = " << abs(f(X_opt, constr + 1) - f(X, constr + 1)) << endl;
    cout << endl;

    ofstr_points << X_opt[0] << " " << X_opt[1] << " " << f(X_opt, constr + 1) << endl;
    ofstr_points << endl << endl;
    ofstr_points << X[0] << " " << X[1] << " " << f(X, constr + 1) << endl;
    ofstr_points << endl << endl;
    mggsa.getPoints(points);
    for (int j = 0; j < points.size(); j++) {
        ofstr_points << points[j][0] << " " << points[j][1] << " " << f(points[j], constr + 1) << endl;
    }
    ofstr_points.close();

    // Plotting the function(works with gnuplot)
    int error;
#if defined(__linux__)
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/map_test.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
#endif

    char str[100];
    sprintf(str, "gnuplot -c scripts/map_test.gp %d", key);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}
