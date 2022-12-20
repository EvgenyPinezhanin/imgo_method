#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

#include <omp.h>
#include <mggsa.h>
#include <map.h>

using namespace std;

#define CALC
// #define OUTPUT_INFO

double euclidean_distance(vector<double> val1, vector<double> val2) {
    double res = 0.0;
    size_t size = val1.size();
    for (int i = 0; i < size; i++) {
        res += (val1[i] - val2[i]) * (val1[i] - val2[i]);
    }
    return sqrt(res);
}

int N;
double f_rastrigin(vector<double> x, int j) {
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += x[i] * x[i] - 10.0 * cos(2.0 * M_PI * x[i]);
    }
    switch (j) {
        case 1: return 10.0 * N + sum;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

const int type = 3; // 0 - count trials, 1 - count points, 
                    // 2 - accuracy, 3 - c_points / c_trials
const int n = 2; // n_min ... n_max
// const int n_min = 2, n_max = 5;
const int n_min = 2, n_max = 3;
const int incr_min = 1, incr_max = 30;
const int m_min = 8, m_max = 12;
const int chunk = 5;

int main() {
#if defined(CALC)
    ofstream ofstr_incr("output_data/incr_test.txt");

    vector<double> X_opt_num, X_opt{0.0}, A{-1.0 / 2.0}, B{1.0};
    double eps = 0.01, r = 2.1, d = 0.0;
    int constr = 0, count, Nmax = 2000;
    Stop stop = Stop::ACCURNUMBER;
    int key = 3;

    vector<vector<vector<double>>> accuracy(n_max - n_min + 1);
    vector<vector<vector<int>>> count_trials(n_max - n_min + 1);
    vector<vector<vector<double>>> count_points(n_max - n_min + 1);
    for (int i = n_min; i <= n_max; i++) {
        accuracy[i - n_min].resize(m_max - m_min + 1);
        count_trials[i - n_min].resize(m_max - m_min + 1);
        count_points[i - n_min].resize(m_max - m_min + 1);
        for (int j = m_min; j <= m_max; j++) {
            accuracy[i - n_min][j - m_min].resize(incr_max - incr_min + 1);
            count_trials[i - n_min][j - m_min].resize(incr_max - incr_min + 1);
            count_points[i - n_min][j - m_min].resize(incr_max - incr_min + 1);
        }
    }

    mggsa_method mggsa(f_rastrigin, -1, constr, A, B, r, d, -1, key, eps, Nmax, -1);

    double t1, t2, dt;

    // X_opt.push_back(0.0);
    // A.push_back(-1.0 / 2.0); 
    // B.push_back(1.0);

    t1 = omp_get_wtime ();
    double accur;
    int count_pnts;
    for (int i = n_min; i <= n_max; i++) {
        cout << "Rastrigin function" << endl;
        cout << "Number of constrained = " << constr << endl;
        cout << "Dimension = " << i << endl;
        cout << "Parameters for method:" << endl;
        cout << "eps = " << eps << " r = " << r << " d = " << d << endl;
        cout << endl;

        N = i;
        mggsa.setN(i);
        X_opt.push_back(0.0);
        A.push_back(-1.0 / 2.0); 
        B.push_back(1.0);
        mggsa.setAB(A, B);

    #pragma omp parallel for schedule(dynamic, chunk) num_threads(omp_get_num_procs()) collapse(2) \
            shared(accuracy, count_trials, count_points, stop) firstprivate(mggsa) private(accur, count_pnts, count, X_opt_num)
        for (int j = m_min; j <= m_max; j++) {
            for (int k = incr_min; k <= incr_max; k++) {
                mggsa.setDen(j);
                mggsa.setIncr(k);
 
                mggsa.solve(count, X_opt_num, stop);
 
                accur = euclidean_distance(X_opt, X_opt_num);
                count_pnts = mggsa.getCountPoints();
                accuracy[i - n_min][j - m_min][k - incr_min] = accur;
                count_trials[i - n_min][j - m_min][k - incr_min] = count;
                count_points[i - n_min][j - m_min][k - incr_min] = count_pnts;

            #if defined(OUTPUT_INFO)
                cout << "Parameters for constructing the Peano curve:" << endl;
                cout << "m = " << j << " key = " << key << " incr = " << k << endl;
                cout << "Trials result:" << endl;
                cout << "Number of trials = " << count << endl;
                cout << "Number of points = " << count_pnts << endl;
                cout << "X* = (";
                for (int l = 0; l < X_opt.size() - 1; l++) {
                    cout << X_opt[l] << ", ";
                }
                cout << X_opt[X_opt.size() - 1] << ")" << endl;
                cout << "X = (";
                for (int l = 0; l < X_opt_num.size() - 1; l++) {
                    cout << X_opt_num[l] << ", ";
                }
                cout << X_opt_num[X_opt_num.size() - 1] << ")" << endl;
                cout << "||X* - X|| = " << accur << endl;
                cout << "f(X*) = " << f_rastrigin(X_opt, constr + 1) << endl;
                cout << "f(X) = " << f_rastrigin(X_opt_num, constr + 1) << endl;
                cout << "|f(X*) - f(X)| = " << abs(f_rastrigin(X_opt, constr + 1) - f_rastrigin(X_opt_num, constr + 1)) << endl;
                cout << endl;
            #endif
            #if not defined(OUTPUT_INFO)
                string str = "n = " + to_string(i) + " m = " + to_string(j) + " incr = " + to_string(k) +
                             " t_num = " + to_string(omp_get_thread_num()) + "\n";
                cout << str;
            #endif
            }
        }
    }

    t2 = omp_get_wtime ();
    dt = t2 - t1;
    cout << "Time: " << dt << endl;

    for (int i = n_min; i <= n_max; i++) {
        for (int j = m_min; j <= m_max; j++) {
            for (int k = incr_min; k <= incr_max; k++) {
                ofstr_incr << k << " " << count_trials[i - n_min][j - m_min][k - incr_min]
                                << " " << count_points[i - n_min][j - m_min][k - incr_min] 
                                << " " << accuracy[i - n_min][j - m_min][k - incr_min] << endl;
            }
            ofstr_incr << endl << endl;
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
    sprintf(str, "gnuplot -c scripts/incr_test.gp %d %d %d %d %d %d %d %d", type, n, n_min, n_max, m_min, m_max, incr_min, incr_max);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }
    return 0;
}
