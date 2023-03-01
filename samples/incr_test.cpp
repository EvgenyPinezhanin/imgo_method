#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS
#endif

#include <fstream>
#include <iostream>
#include <vector>
#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include <math.h>
#else
    #include <cmath>
#endif

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
const int n_type = 2; // n_min ... n_max

const int n_min = 2, n_max = 3;
const int n_count = n_max - n_min + 1;
const int incr_array[n_count][2] = { {1, 60},
                                     {1, 100} };
const int m_min = 8, m_max = 10;

const int chunk = 2;

int main() {
    ofstream ofstr_opt("output_data/incr_test_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";
#if defined(CALC)
    ofstream ofstr("output_data/incr_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";

    vector<double> X, X_opt, A, B;
    double eps = 0.01, r = 2.1, d = 0.0;
    int constr = 0, Nmax = 10000, key = 3;
    Stop stop = Stop::ACCURNUMBER;

    vector<vector<vector<double>>> accuracy_vec(n_count), count_points_vec(n_count);
    vector<vector<vector<int>>> count_trials_vec(n_count);
    for (int i = n_min; i <= n_max; i++) {
        accuracy_vec[i - n_min].resize(m_max - m_min + 1);
        count_trials_vec[i - n_min].resize(m_max - m_min + 1);
        count_points_vec[i - n_min].resize(m_max - m_min + 1);
        for (int j = m_min; j <= m_max; j++) {
            accuracy_vec[i - n_min][j - m_min].resize(incr_array[i - n_min][1] - incr_array[i - n_min][0] + 1);
            count_trials_vec[i - n_min][j - m_min].resize(incr_array[i - n_min][1] - incr_array[i - n_min][0] + 1);
            count_points_vec[i - n_min][j - m_min].resize(incr_array[i - n_min][1] - incr_array[i - n_min][0] + 1);
        }
    }

    mggsa_method mggsa(f_rastrigin, -1, constr, A, B, r, d, -1, key, eps, Nmax, -1);

    for (int i = 0; i < n_min - 1; i++) {
        X_opt.push_back(0.0);
        A.push_back(-1.0 / 2.0);
        B.push_back(1.0);
    }

    double accuracy;
    int count_points, count_trials;

    int total_start_time = clock();
    for (int i = n_min; i <= n_max; i++) {
        N = i;
        mggsa.setN(N);
        X_opt.push_back(0.0);
        A.push_back(-1.0 / 2.0);
        B.push_back(1.0);
        mggsa.setAB(A, B);

    #if defined(OUTPUT_INFO)
        cout << "Rastrigin function" << endl;
        cout << "Dimension = " << N << endl;
        cout << "Number of constrained = " << constr << endl;
        cout << "[A; B] = [(";
        for (int j = 0; j < A.size() - 1; j++) {
            cout << A[j] << ", "; 
        }
        cout << A[A.size() - 1] << "); (" << endl;
        for (int j = 0; j < B.size() - 1; j++) {
            cout << B[j] << ", "; 
        }
        cout << B[B.size() - 1] << ")]" << endl;
        cout << "X* = (";
        for (int l = 0; l < X_opt.size() - 1; l++) {
            cout << X_opt[l] << ", ";
        }
        cout << X_opt[X_opt.size() - 1] << ")" << endl;
        cout << "f(X*) = " << f_rastrigin(X_opt, constr + 1) << endl;
        cout << "Parameters for method:" << endl;
        cout << "eps = " << eps << " r = " << r << " d = " << d << endl;
        cout << "Parameters for constructing the Peano curve:" << endl;
        cout << "key = " << key << endl;
        cout << endl;
    #endif

        vector<double> mu;
    #pragma omp parallel for schedule(dynamic, chunk) proc_bind(spread) num_threads(omp_get_num_procs()) collapse(2) \
            shared(incr_array, accuracy_vec, count_trials_vec, count_points_vec, stop) \
            firstprivate(mggsa) private(accuracy, count_points, count_trials, X, mu)
        for (int j = m_min; j <= m_max; j++) {
            for (int k = incr_array[i - n_min][0]; k <= incr_array[i - n_min][1]; k++) {
                double start_time = omp_get_wtime();

                mggsa.setDen(j);
                mggsa.setIncr(k);
 
                mggsa.solve(count_trials, X, stop);
                mggsa.getLambda(mu);
 
                accuracy = euclidean_distance(X_opt, X);
                count_points = mggsa.getCountPoints();
                accuracy_vec[i - n_min][j - m_min][k - incr_array[i - n_min][0]] = accuracy;
                count_trials_vec[i - n_min][j - m_min][k - incr_array[i - n_min][0]] = count_trials;
                count_points_vec[i - n_min][j - m_min][k - incr_array[i - n_min][0]] = count_points;

                double end_time = omp_get_wtime();
                double work_time = end_time - start_time;

            #if defined(OUTPUT_INFO)
                cout << "Parameters for constructing the Peano curve:" << endl;
                cout << "m = " << j << " incr = " << k << endl;
                cout << "Trials result:" << endl;
                cout << "Number of trials = " << count_trials << endl;
                cout << "Number of points = " << count_points << endl;
                cout << "Estimation of the Lipschitz constant = " << mu[0] << endl;
                cout << "X = (";
                for (int l = 0; l < X.size() - 1; l++) {
                    cout << X[l] << ", ";
                }
                cout << X[X.size() - 1] << ")" << endl;
                cout << "f(X) = " << f_rastrigin(X, constr + 1) << endl;
                cout << "||X* - X|| = " << accuracy << endl;

                cout << "|f(X*) - f(X)| = " << abs(f_rastrigin(X_opt, constr + 1) - f_rastrigin(X, constr + 1)) << endl;
                cout << endl;
            #else
                string str = "Rastrigin: n = " + to_string(i) + " m = " + to_string(j) + " incr = " + to_string(k) +
                             " time: " + to_string(work_time) + " t_num = " + to_string(omp_get_thread_num()) + "\n";
                cout << str;
            #endif
            }
        }
    }
    int total_end_time = clock();
    double total_work_time = (double)(total_end_time - total_start_time) / CLOCKS_PER_SEC;
    cout << "Total time: " << total_work_time << endl;

    for (int i = n_min; i <= n_max; i++) {
        for (int j = m_min; j <= m_max; j++) {
            for (int k = incr_array[i - n_min][0]; k <= incr_array[i - n_min][1]; k++) {
                ofstr << k << " " << count_trials_vec[i - n_min][j - m_min][k - incr_array[i - n_min][0]]
                                << " " << count_points_vec[i - n_min][j - m_min][k - incr_array[i - n_min][0]] 
                                << " " << accuracy_vec[i - n_min][j - m_min][k - incr_array[i - n_min][0]] << endl;
            }
            ofstr << endl << endl;
        }
        ofstr << endl << endl;
    }
    ofstr.close();
#endif

    ofstr_opt << "array I_MIN[" << n_count << "]" << endl;
    ofstr_opt << "array I_MAX[" << n_count << "]" << endl;
    for (int i = 0; i < n_count; i++) {
        ofstr_opt << "I_MIN[" << i + 1 << "] = " << incr_array[i][0] << endl;
        ofstr_opt << "I_MAX[" << i + 1 << "] = " << incr_array[i][1] << endl;
    }
    ofstr_opt.close();

    int error;
#if defined(__linux__)
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/incr_test.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
#endif

    char str[100];
    sprintf(str, "gnuplot -c scripts/incr_test.gp %d %d %d %d %d %d", type, n_type, n_min, n_max, m_min, m_max);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

    return 0;
}
