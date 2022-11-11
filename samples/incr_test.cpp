#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

#include <mggsa.h>
#include <map.h>

using namespace std;

// #define CALC

double f(vector<double> x, int j) {
    switch (j) {
        case 1: return x[0] * x[0] + x[1] * x[1] - cos(18.0 * x[0]) - cos(18.0 * x[1]);
        default: return numeric_limits<double>::quiet_NaN();
    }
}

const int type = 0; // 0 - count trials, 1 - accuracy
const int incr_min = 1, incr_max = 40;
const int m_min = 5, m_max = 12;

int main() {
#if defined(CALC)
    ofstream ofstr_incr("incr_test.txt");

    vector<double> X_opt_num(2), A{-1.0 / 2.0, -1.0 / 2.0}, B{1.0, 1.0}, X_opt{0.0, 0.0};
    double eps = 0.01, r = 2.1, d = 0.0;
    int constr = 0, count, Nmax = 1000;
    Stop stop = Stop::ACCURNUMBER;
    int n = 2, key = 3;

    vector<vector<double>> accuracy(m_max - m_min + 1);
    vector<vector<int>> count_trials(m_max - m_min + 1);
    for (int i = m_min; i <= m_max; i++) {
        accuracy[i - m_min].resize(incr_max - incr_min + 1);
        count_trials[i - m_min].resize(incr_max - incr_min + 1);
    }

    mggsa_method mggsa(f, n, constr, A, B, r, d, -1, key, eps, Nmax, -1);

    cout << "Function: " << "x^2 + y^2 - cos(18.0 * x) - cos(18.0 * y)" << endl;
    cout << "Number of constrained = " << constr << endl;
    cout << "Dimension = " << n << endl;
    cout << "Parameters for method:" << endl;
    cout << "eps = " << 0.01 << " r = " << 2.0 << " d = " << 0.0 << endl;
    cout << endl;

    double accur;
    for (int i = m_min; i <= m_max; i++) {
        for (int j = incr_min; j <= incr_max; j++) {
            mggsa.setDen(i);
            mggsa.setIncr(j);

            mggsa.solve(count, X_opt_num, stop);

            accur = sqrt((X_opt[0] - X_opt_num[0]) * (X_opt[0] - X_opt_num[0]) + 
                                            (X_opt[1] - X_opt_num[1]) * (X_opt[1] - X_opt_num[1]));
            accuracy[i - m_min][j - incr_min] = accur;
            count_trials[i - m_min][j - incr_min] = count;

            cout << "Parameters for constructing the Peano curve:" << endl;
            cout << "m = " << i << " key = " << key << " incr = " << j << endl;
            cout << "Trials result:" << endl;
            cout << "Number of trials = " << count << endl;
            cout << "X* = " << X_opt[0] << " Y* = " << X_opt[1] << endl;
            cout << "X = " << X_opt_num[0] << " Y = " << X_opt_num[1] << endl;
            cout << "||X* - X|| = " << accur << endl;
            cout << "f(X*) = " << f(X_opt, constr + 1) << endl;
            cout << "f(X) = " << f(X_opt_num, constr + 1) << endl;
            cout << "|f(X*) - f(X)| = " << abs(f(X_opt, constr + 1) - f(X_opt_num, constr + 1)) << endl;
            cout << endl;
        }
    }

    for (int i = m_min; i <= m_max; i++) {
        for (int j = incr_min; j <= incr_max; j++) {
            ofstr_incr << j << " " << count_trials[i - m_min][j - incr_min]
                            << " " << accuracy[i - m_min][j - incr_min] << endl;
        }
        ofstr_incr << endl << endl;
    }
#endif

    int error;
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/incr_test.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
    char str[100];
    sprintf(str, "gnuplot -c scripts/incr_test.gp %d %d %d %d % d", type, incr_min, incr_max, m_min, m_max);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }
    return 0;
}

    // ofstr_count << incr_max << " ";
    // ofstr_accur << incr_max << " ";
    // for (int j = 0; j < incr_max; j++) {
	// 	ofstr_count << j + 1 << " ";
    //     ofstr_accur << j + 1 << " ";
	// }
    // ofstr_count << endl;
    // ofstr_accur << endl;
    // for (int i = 0; i < m_max; i++) {
    //     ofstr_count << i + 1 << " ";
    //     ofstr_accur << i + 1 << " ";
    //     for (int j = 0; j < incr_max; j++) {
    //         ofstr_count << count_trials[i][j] << " ";
    //         ofstr_accur << accuracy[i][j] << " ";
    //     }
    //     ofstr_count << endl;
    //     ofstr_accur << endl;
    // }
	// ofstr_count.close();
    // ofstr_accur.close();