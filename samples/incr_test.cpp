#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

#include <omp.h>
#include <mggsa.h>
#include <map.h>

using namespace std;

// #define CALC
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
const int n = 3; // n_min ... n_max
const int n_min = 2, n_max = 3;
const int n_count = 2;
const int incr_array[n_count][2] = { {1, 60},
                                     {1, 100} };
const int m_min = 8, m_max = 10;
const int chunk = 2;

int main() {
    ofstream ofstr_opt("output_data/incr_test_opt.txt");
#if defined(CALC)
    ofstream ofstr_incr("output_data/incr_test.txt");

    vector<double> X_opt_num, X_opt, A, B;
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
            accuracy[i - n_min][j - m_min].resize(incr_array[i - n_min][1] - incr_array[i - n_min][0] + 1);
            count_trials[i - n_min][j - m_min].resize(incr_array[i - n_min][1] - incr_array[i - n_min][0] + 1);
            count_points[i - n_min][j - m_min].resize(incr_array[i - n_min][1] - incr_array[i - n_min][0] + 1);
        }
    }

    mggsa_method mggsa(f_rastrigin, -1, constr, A, B, r, d, -1, key, eps, Nmax, -1);

    for (int i = 0; i < n_min - 1; i++) {
        X_opt.push_back(0.0);
        A.push_back(-1.0 / 2.0); 
        B.push_back(1.0);
    }

    double t1, t2, dt;
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

    #pragma omp parallel for schedule(dynamic, chunk) proc_bind(spread) num_threads(omp_get_num_procs()) collapse(2) \
            shared(incr_array, accuracy, count_trials, count_points, stop) firstprivate(mggsa) private(accur, count_pnts, count, X_opt_num)
        for (int j = m_min; j <= m_max; j++) {
            for (int k = incr_array[i - n_min][0]; k <= incr_array[i - n_min][1]; k++) {
                mggsa.setDen(j);
                mggsa.setIncr(k);
 
                mggsa.solve(count, X_opt_num, stop);
 
                accur = euclidean_distance(X_opt, X_opt_num);
                count_pnts = mggsa.getCountPoints();
                accuracy[i - n_min][j - m_min][k - incr_array[i - n_min][0]] = accur;
                count_trials[i - n_min][j - m_min][k - incr_array[i - n_min][0]] = count;
                count_points[i - n_min][j - m_min][k - incr_array[i - n_min][0]] = count_pnts;

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
            for (int k = incr_array[i - n_min][0]; k <= incr_array[i - n_min][1]; k++) {
                ofstr_incr << k << " " << count_trials[i - n_min][j - m_min][k - incr_array[i - n_min][0]]
                                << " " << count_points[i - n_min][j - m_min][k - incr_array[i - n_min][0]] 
                                << " " << accuracy[i - n_min][j - m_min][k - incr_array[i - n_min][0]] << endl;
            }
            ofstr_incr << endl << endl;
        }
        ofstr_incr << endl << endl;
    }
    ofstr_incr.close();

#endif

    ofstr_opt << "array I_MIN[" << n_count << "]" << endl;
    ofstr_opt << "array I_MAX[" << n_count << "]" << endl;
    for (int i = 0; i < n_count; i++) {
        ofstr_opt << "I_MIN[" << i + 1 << "] = " << incr_array[i][0] << endl;
        ofstr_opt << "I_MAX[" << i + 1 << "] = " << incr_array[i][1] << endl;
    }
    ofstr_opt.close();

    int error;
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/incr_test.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
    char str[100];
    sprintf(str, "gnuplot -c scripts/incr_test.gp %d %d %d %d %d %d", type, n, n_min, n_max, m_min, m_max);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }
    return 0;
}
