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

const int key = 3;

double f(vector<double> x, int j) {
    switch (j) {
        case 1: return x[0] * x[0] + x[1] * x[1] - cos(18.0 * x[0]) - cos(18.0 * x[1]);
        default: return numeric_limits<double>::quiet_NaN();
    }
}

int main() {
    ofstream ofstr("output_data/mggsa_test_evolvent.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstrPoints("output_data/mggsa_test_evolvent_points.txt");
    if (!ofstrPoints.is_open()) cerr << "File opening error\n";

    vector<double> A{ -0.5, -0.5 }, B{ 1.0, 1.0 };
    double eps = 0.01, r = 2.0, d = 0.0;
    int numberConstraints = 0;
    int maxTrials = 100000, maxFevals = 100000;
    int den = 8, n = 2, incr = 0;

    MggsaMethod mggsa(f, n, numberConstraints, A, B, r, d, den, key, eps, maxTrials, maxFevals, incr);

    vector<double> XMap(n);
    double k = (key != 3) ? 1.0 / (pow(2.0, n * den) - 1.0) :
                            1.0 / (pow(2.0, den * n) * (pow(2.0, n) - 1.0)) + 0.0000000001;
    if (key == 3) den++;
    for (double i = 0.0; i <= 1.0; i += k) {
        mapd(i, den, XMap.data(), n, key);
        ofstr << XMap[0] * (B[0] - A[0]) + (A[0] + B[0]) / 2.0 << " " 
              << XMap[1] * (B[1] - A[1]) + (A[1] + B[1]) / 2.0 << endl;
    }
    mapd(1.0, den, XMap.data(), n, key);
    ofstr << XMap[0] * (B[0] - A[0]) + (A[0] + B[0]) / 2.0 << " "
          << XMap[1] * (B[1] - A[1]) + (A[1] + B[1]) / 2.0 << endl;
    if (key == 3) den--;
    ofstr.close();

    vector<double> XOpt{ 0.0, 0.0 }, X, L;
    int numberTrials, numberFevals;
    vector<vector<double>> points;
    vector<TrialConstrained> trials;

    mggsa.solve(numberTrials, numberFevals, X);
    mggsa.getL(L);

    printResultMggsa("x^2 + y^2 - cos(18.0 * x) - cos(18.0 * y)", n, numberConstraints, A, B, vector<double>(), XOpt,
                     f(XOpt, numberConstraints + 1), maxTrials, maxFevals, eps, r, d, den, key, incr, numberTrials,
                     numberFevals, L, X, f(X, numberConstraints + 1));

    addPointGnuplot(ofstrPoints, XOpt, f(XOpt, numberConstraints + 1));
    addPointGnuplot(ofstrPoints, X, f(X, numberConstraints + 1));

    mggsa.getPoints(points);
    mggsa.getTrialPoints(trials);
    addPointsGnuplot(ofstrPoints, points, trials);

    ofstrPoints.close();

    drawGraphGnuplot("scripts/mggsa_test_evolvent.gp", key);

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}
