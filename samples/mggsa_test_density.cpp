#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <mggsa.h>
#include <map.h>
#include <output_results.h>

using namespace std;

const int displayType = 1; // 0 - application, 1 - gif

double f(vector<double> x, int j) {
    switch (j) {
        case 1: return x[0] * x[0] + x[1] * x[1] - cos(18.0 * x[0]) - cos(18.0 * x[1]);
        default: return numeric_limits<double>::quiet_NaN();
    }
}

int main() {
    ofstream ofstr("output_data/mggsa_test_density.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstrOpt("output_data/mggsa_test_density_opt.txt");
    if (!ofstrOpt.is_open()) cerr << "File opening error\n";

    vector<double> A{ -0.5, -0.5 }, B{ 1.0, 1.0 };
    double eps = 0.01, r = 2.0, d = 0.0;
    int numberConstraints = 0;
    int maxTrials = 100000, maxFevals = 100000;
    int n = 2, key = 3, incr = 0;

    MggsaMethod mggsa(f, n, numberConstraints, A, B, r, d, -1, key, eps, maxTrials, maxFevals, incr);

    vector<double> XOpt{ 0.0, 0.0 }, X;
    vector<int> den{ 4, 6, 8, 10, 12 };
    int numberTrials, numberFevals;
    vector<double> L;
    vector<vector<double>> points;
    vector<TrialConstrained> trials;

    addPointGnuplot(ofstr, XOpt, f(XOpt, numberConstraints + 1));

    for (int i = 0; i < den.size(); i++) {
        mggsa.setDen(den[i]);

        if (i == 0) {
            mggsa.solve(numberTrials, numberFevals, X);
        } else {
            mggsa.solve(numberTrials, numberFevals, X, TypeSolve::RESOLVE);
        }
        mggsa.getL(L);

        printResultMggsa("x^2 + y^2 - cos(18.0 * x) - cos(18.0 * y)", n, numberConstraints, A, B, vector<double>(), XOpt,
                         f(XOpt, numberConstraints + 1), maxTrials, maxFevals, eps, r, d, den[i], key, -1, numberTrials,
                         numberFevals, L, X, f(X, numberConstraints + 1));

        addPointGnuplot(ofstr, X, f(X, numberConstraints + 1));

        mggsa.getPoints(points);
        addPointsGnuplot(ofstr, points);
    }
    ofstr.close();

    size_t sizeDen = den.size();
    initArrayGnuplot(ofstrOpt, "den", sizeDen);
    for (int i = 0; i < den.size(); i++) {
        setValueInArrayGnuplot(ofstrOpt, "den", i + 1, den[i]);
    }
    ofstrOpt.close();

    vector<int> args{ displayType, (int)den.size() };
    drawGraphGnuplot("scripts/mggsa_test_density.gp", args);

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}
