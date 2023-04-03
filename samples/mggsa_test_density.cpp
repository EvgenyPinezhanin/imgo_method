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
    ofstream ofstr("output_data/mggsa_test_density.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("output_data/mggsa_test_density_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    vector<double> A{-0.5, -0.5}, B{1.0, 1.0}, X_opt{0.0, 0.0};
    double eps = 0.01, r = 2.0, d = 0.0;
    int constr = 0;
    int countIters, countTrials, countEvals;
    int maxIters = 100000, maxEvals = 100000;
    int n = 2, key = 3, incr = 10;
    vector<int> den{ 8, 10, 12, 14, 16 };
    vector<double> X;

    mggsa_method mggsa(f, n, constr, A, B, r, d, -1, key, eps, maxIters, maxEvals, incr);

    vector<double> mu;
    vector<vector<double>> points;

    cout << "Function: " << "x^2 + y^2 - cos(18.0 * x) - cos(18.0 * y)" << endl;
    cout << "Dimension = " << n << endl;
    cout << "Number of constrained = " << constr << endl;
    cout << "[A; B] = [(" << A[0] << ", " << A[1] << "); (" << 
                             B[0] << ", " << B[1] << ")]"<< endl;
    cout << "X* = (" << X_opt[0] << ", " << X_opt[1] << ")" << endl;
    cout << "f(X*) = " << f(X_opt, constr + 1) << endl;
    cout << "Parameters for method:" << endl;
    cout << "eps = " << eps << " r = " << r << " d = " << d << endl;

    ofstr << X_opt[0] << " " << X_opt[1] << " " << f(X_opt, constr + 1) << endl;
    ofstr << endl << endl;

    for (int i = 0; i < den.size(); i++) {
        mggsa.setDen(den[i]);

        if (i == 0) {
            mggsa.solve(countIters, countTrials, countEvals, X);
        } else {
            mggsa.solve(countIters, countTrials, countEvals, X, TypeSolve::RESOLVE);
        }
        mggsa.getLambda(mu);

        cout << "Parameters for constructing the Peano curve:" << endl;
        cout << "m = " << den[i] << " key = " << key << " incr = " << incr << endl;
        cout << "Trials result:" << endl;
        cout << "Number of iters = " << countIters << endl;
        cout << "Number of trials = " << countTrials << endl;
        cout << "Number of evals = " << countEvals << endl;
        cout << "Estimation of the Lipschitz constant = " << mu[0] << endl;
        cout << "X = (" << X[0] << ", " << X[1] << ")" << endl;
        cout << "f(X) = " << f(X, constr + 1) << endl;
        cout << "||X* - X|| = " << sqrt((X_opt[0] - X[0]) * (X_opt[0] - X[0]) + 
                                        (X_opt[1] - X[1]) * (X_opt[1] - X[1])) << endl;
        cout << "|f(X*) - f(X)| = " << abs(f(X_opt, constr + 1) - f(X, constr + 1)) << endl;
        cout << endl;

        ofstr << X[0] << " " << X[1] << " " << f(X, constr + 1) << endl;
        ofstr << endl << endl;
        mggsa.getPoints(points);
        for (int j = 0; j < points.size(); j++) {
            ofstr << points[j][0] << " " << points[j][1] << " " << f(points[j], constr + 1) << endl;
        }
        ofstr << endl << endl;
    }
    ofstr.close();

    ofstr_opt << "array Den[" << den.size() << "]" << endl;
    for (int i = 0; i < den.size(); i++) {
        ofstr_opt << "Den[" << i + 1 << "]=\"" << den[i] << "\"" << endl;
    }
    ofstr_opt.close();

    // Plotting the function(works with gnuplot)
    int error;
#if defined(__linux__)
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/mggsa_test_density.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
#endif

    char str[100];
    sprintf(str, "gnuplot -c scripts/mggsa_test_density.gp %d", (int)den.size());
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}
