#include <iostream>
#include <fstream>
#include <vector>

#include <omp.h>
#include <Grishagin/grishagin_function.hpp>
#include <Grishagin/GrishaginConstrainedProblem.hpp>
#include <GKLS/GKLSProblem.hpp>
#include <GKLS/GKLSConstrainedProblem.hpp>
#include <mggsa.h>
#include <map.h>
#include <task.h>

using namespace std;

using vector_4d = vector<vector<vector<vector<double>>>>;

// #define CALC
// #define OUTPUT_INFO

void calculation(mggsa_method &mggsa, vector_4d &lipschitz_const, problem_single problem, int num_func,
                                                                     double r, int key, int m, int incr);

const int type = 0; // 0 - grishagin, 1 - GKLS
                    // 2 - constrained grisagin, 3 - constrained GKLS

const int incr_min = 1, incr_max = 40;
const int m_min = 8, m_max = 12;
const int key_min = 1, key_max = 3;

const int chunk = 4;

int main() {
    ofstream ofstr_opt("output_data/lipschitz_test_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";
    
    const int count_func = 4;

    TGrishaginProblem grishaginProblem;
    GrishaginConstrainedProblem grishaginConstrainedProblem;
    TGKLSProblem gklsProblem;
    TGKLSConstrainedProblem gklsConstrainedProblem;

    vector<problem_single> problems{ problem_single("Grishagin", &grishaginProblem, type_constraned::NONCONSTR),
                                     problem_single("GKLS", &gklsProblem, type_constraned::NONCONSTR),
                                     problem_single("GrishaginConstrained", &grishaginConstrainedProblem, type_constraned::CONSTR),
                                     problem_single("GKLSConstrained", &gklsConstrainedProblem, type_constraned::CONSTR) };

#if defined(CALC)
    ofstream ofstr("output_data/lipschitz_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";

    double eps = 0.01, d = 0.0, r = 2.0;
    int maxIters = 100000, maxEvals = 100000;
    vector<double> r_vec{2.8, 3.5, 2.8, 2.8};

    mggsa_method mggsa(nullptr, -1, -1, vector<double>(), vector<double>(), r, d, -1, -1, eps, maxIters, maxEvals, -1);

    vector_4d lipschitz_const(count_func);
    for (int i = 0; i < count_func; i++) {
        lipschitz_const[i].resize(key_max - key_min + 1);
        for (int j = 0; j < key_max - key_min + 1; j++) {
            lipschitz_const[i][j].resize(m_max - m_min + 1);
            for (int k = 0; k < m_max - m_min + 1; k++) {
                lipschitz_const[i][j][k].resize(incr_max - incr_min + 1);
            }
        }
    }

    cout << "Parameters for method:" << endl;
    cout << "eps = " << eps << " r = " << r << " d = " << d << endl;

    int total_start_time = clock();
#pragma omp parallel for schedule(dynamic, chunk) proc_bind(spread) num_threads(omp_get_num_procs()) collapse(4) \
        shared(count_func, key_min, key_max, m_min, m_max, incr_min, incr_max, problems, r_vec, lipschitz_const) \
        firstprivate(mggsa)
    for (int i = 0; i < count_func; i++) {
        for (int j = key_min; j <= key_max; j++) {
            for (int k = m_min; k <= m_max; k++) {
                for (int l = incr_min; l <= incr_max; l++) {
                    if (j != key_max) {
                        if (l == incr_min) {
                            calculation(mggsa, lipschitz_const, problems[i], i, r_vec[i], j, k, l);
                            for (int m = incr_min + 1; m <= incr_max; m++) {
                                lipschitz_const[i][j - key_min][k - m_min][m - incr_min] = 
                                    lipschitz_const[i][j - key_min][k - m_min][0];
                            }
                        }
                    } else {
                        calculation(mggsa, lipschitz_const, problems[i], i, r_vec[i], j, k, l);
                    }
                }
            }
        }
    }
    int total_end_time = clock();
    double total_work_time = (double)(total_end_time - total_start_time) / CLOCKS_PER_SEC;
    cout << "Total time: " << total_work_time << endl;

    for (int i = 0; i < count_func; i++) {
        for (int j = key_min; j <= key_max; j++) {
            for (int k = m_min; k <= m_max; k++) {
                for (int l = incr_min; l <= incr_max; l++) {
                    ofstr << l << " " << lipschitz_const[i][j - key_min][k - m_min][l - incr_min] << endl;
                }
                ofstr << endl << endl;
            }
        }
    }
    ofstr.close();
#endif

    ofstr_opt << "array Name[" << count_func << "]" << endl;
    for (int i = 0; i < count_func; i++) {
        ofstr_opt << "Name[" << i + 1 << "] = \"" << problems[i].name << "\"" << endl;
    }
    ofstr_opt.close();

    int error;
#if defined(__linux__)
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/lipschitz_test.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
#endif

    char str[100];
    sprintf(str, "gnuplot -c scripts/lipschitz_test.gp %d %d %d %d %d %d %d", type, key_min, key_max, m_min, m_max, incr_min, incr_max);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

    return 0;
}

void calculation(mggsa_method &mggsa, vector_4d &lipschitz_const, problem_single problem, int num_func,
                                                                     double r, int key, int m, int incr) {
    double accuracy, f_X_opt, f_X;                                                            
    int constr, count_points, n;
    int countIters, countTrials, countEvals;
    vector<double> A, B, X_opt, X, mu;
    functor_single func;
    functor_single_constr func_constr;

    double start_time = omp_get_wtime();

    if (problem.type == type_constraned::CONSTR) {
        func_constr.constr_opt_problem = static_cast<IConstrainedOptProblem*>(problem.optProblem);
        func_constr.constr_opt_problem->GetBounds(A, B);
        n = func_constr.constr_opt_problem->GetDimension();
        constr = func_constr.constr_opt_problem->GetConstraintsNumber();
        X_opt = func_constr.constr_opt_problem->GetOptimumPoint();
        mggsa.setF(func_constr);
        f_X_opt = func_constr(X_opt, constr + 1);
    } else {
        func.opt_problem = static_cast<IOptProblem*>(problem.optProblem);
        func.opt_problem->GetBounds(A, B);
        n = func.opt_problem->GetDimension();
        constr = 0;
        X_opt = func.opt_problem->GetOptimumPoint();
        mggsa.setF(func);
        f_X_opt = func(X_opt, 1);
    }

    mggsa.setM(constr);
    mggsa.setN(n);
    mggsa.setAB(A, B);
    mggsa.setKey(key);
    mggsa.setDen(m);
    mggsa.setIncr(incr);
    mggsa.setR(r);

    mggsa.solve(countIters, countTrials, countEvals, X);
    mggsa.getLambda(mu);

    if (problem.type == type_constraned::CONSTR) {
        f_X = func_constr(X, constr + 1);
    } else {
        f_X = func(X, 1);
    }
    lipschitz_const[num_func][(size_t)key - key_min][(size_t)m - m_min][(size_t)incr - incr_min] = mu[0];
    accuracy = sqrt((X_opt[0] - X[0]) * (X_opt[0] - X[0]) + 
                    (X_opt[1] - X[1]) * (X_opt[1] - X[1]));

    double end_time = omp_get_wtime();
    double work_time = end_time - start_time;

#if defined(OUTPUT_INFO)
    cout << "Function: " << problem.name << endl;
    cout << "Dimension = " << n << endl;
    cout << "Number of constrained = " << constr << endl;
    cout << "[A; B] = [(" << A[0] << ", " << A[1] << "); (" <<
                             B[0] << ", " << B[1] << ")]"<< endl;
    cout << "X* = (" << X_opt[0] << ", " << X_opt[1] << ")" << endl;
    cout << "f(X*) = " << f_X_opt << endl;
    cout << "Parameters for constructing the Peano curve:" << endl;
    cout << "m = " << m << " key = " << key << " incr = " << incr << endl;
    cout << "Trials result:" << endl;
    cout << "Number of iters = " << countIters << endl;
    cout << "Number of trials = " << countTrials << endl;
    cout << "Number of evals = " << countEvals << endl;
    cout << "Estimation of the Lipschitz constant:" << endl;
    cout << "L(f(x)) = " << mu[constr] << endl;
    for (int j = 0; j < constr; j++) {
        cout << "L(g" << j + 1 << ") = " << mu[j] << endl;
    }
    cout << "X = (" << X[0] << ", " << X[1] << ")" << endl;
    cout << "f(X) = " << f_X << endl;
    cout << "||X* - X|| = " << accuracy << endl;
    cout << "|f(X*) - f(X)| = " << abs(f_X_opt - f_X) << endl;
    cout << endl;
#else
    string str_output = problem.name + " key = " + to_string(key) + " m = " + to_string(m) + " incr = " + to_string(incr) +
                        " time: " + to_string(work_time) + " t_num = " + to_string(omp_get_thread_num()) + "\n";
    cout << str_output;
#endif
}
