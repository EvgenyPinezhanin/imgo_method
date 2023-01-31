#include <fstream>
#include <iostream>
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

#define CALC
// #define OUTPUT_INFO

TGrishaginProblem grishaginProblem;
double f_grishagin(vector<double> x, int j) {
    switch (j) {
        case 1: return grishaginProblem.ComputeFunction(x);
        default: return numeric_limits<double>::quiet_NaN();
    }
};

GrishaginConstrainedProblem grishaginConstrainedProblem;
double f_constr_grishagin(vector<double> x, int j) {
    int constr = grishaginConstrainedProblem.GetConstraintsNumber();
    if (j >= 1 && j <= constr) {
        return grishaginConstrainedProblem.ComputeConstraint(j - 1, x);
    } else if (j - 1 == constr) {
        return grishaginConstrainedProblem.ComputeFunction(x);
    } else {
        return numeric_limits<double>::quiet_NaN();
    }
};

TGKLSProblem GKLSProblem;
double f_gkls(vector<double> x, int j) {
    switch (j) {
        case 1: return GKLSProblem.ComputeFunction({ x });
        default: return numeric_limits<double>::quiet_NaN();
    }
};

TGKLSConstrainedProblem GKLSConstrainedProblem;
double f_constr_gkls(vector<double> x, int j) {
    int constr = GKLSConstrainedProblem.GetConstraintsNumber();
    if (j >= 1 && j <= constr) {
        return GKLSConstrainedProblem.ComputeConstraint(j - 1, x);
    } else if (j - 1 == constr) {
        return GKLSConstrainedProblem.ComputeFunction(x);
    } else {
        return numeric_limits<double>::quiet_NaN();
    }
};

void calculation(mggsa_method &mggsa, vector_4d &lipschitz_const, class_problems_om problem, int num_func,
                                                            double r, int key, int m, int incr, Stop stop);

const int type = 0; // 0 - grishagin, 1 - GKLS
                    // 2 - constrained grisagin, 3 - constrained GKLS

const int incr_min = 1, incr_max = 40;
const int m_min = 8, m_max = 12;
const int key_min = 1, key_max = 3;

const int chunk = 4;

int main() {
#if defined(CALC)
    ofstream ofstr_lipschitz("output_data/lipschitz_test.txt");

    int count_func = 4;
    vector<class_problems_om> problems{ class_problems_om("Grishagin", &grishaginProblem, type_constraned::NONCONSTR, f_grishagin),
                                        class_problems_om("GKLS", &GKLSProblem, type_constraned::NONCONSTR, f_gkls),
                                        class_problems_om("GrishaginConstrained", &grishaginConstrainedProblem, 
                                                                                      type_constraned::CONSTR, f_constr_grishagin),
                                        class_problems_om("GKLSConstrained", &GKLSConstrainedProblem, 
                                                                                           type_constraned::CONSTR, f_constr_gkls) };
    double eps = 0.01, d = 0.0, r = 2.0;
    int Nmax = 2000;
    Stop stop = Stop::ACCURNUMBER;
    vector<double> r_vec{2.8, 3.5, 2.8, 2.8};

    mggsa_method mggsa(nullptr, -1, -1, vector<double>(), vector<double>(), r, d, -1, -1, eps, Nmax, -1);

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

#pragma omp parallel for schedule(dynamic, chunk) proc_bind(spread) num_threads(omp_get_num_procs()) collapse(4) \
    shared(count_func, key_min, key_max, m_min, m_max, incr_min, incr_max, problems, r_vec, lipschitz_const, stop) \
    firstprivate(mggsa)
    for (int i = 0; i < count_func; i++) {
        for (int j = key_min; j <= key_max; j++) {
            for (int k = m_min; k <= m_max; k++) {
                for (int l = incr_min; l <= incr_max; l++) {
                    if (j != key_max) {
                        if (l == incr_min) {
                            calculation(mggsa, lipschitz_const, problems[i], i, r_vec[i], j, k, l, stop);
                            for (int m = incr_min + 1; m <= incr_max; m++) {
                                lipschitz_const[i][j - key_min][k - m_min][m - incr_min] = 
                                    lipschitz_const[i][j - key_min][k - m_min][0];
                            }
                        }
                    } else {
                        calculation(mggsa, lipschitz_const, problems[i], i, r_vec[i], j, k, l, stop);
                    }
                }
            }
        }
    }

    for (int i = 0; i < count_func; i++) {
        for (int j = key_min; j <= key_max; j++) {
            for (int k = m_min; k <= m_max; k++) {
                for (int l = incr_min; l <= incr_max; l++) {
                    ofstr_lipschitz << l << " " << lipschitz_const[i][j - key_min][k - m_min][l - incr_min] << endl;
                }
                ofstr_lipschitz << endl << endl;
            }
        }
    }
    ofstr_lipschitz.close();

#endif

    int error;
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/lipschitz_test.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }

    char str[100];
    sprintf(str, "gnuplot -c scripts/lipschitz_test.gp %d %d %d %d %d %d %d", type, key_min, key_max, m_min, m_max, incr_min, incr_max);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

    return 0;
}

void calculation(mggsa_method &mggsa, vector_4d &lipschitz_const, class_problems_om problem, int num_func,
                                                            double r, int key, int m, int incr, Stop stop) {
    double accuracy;                                                            
    int constr, count_trials, count_points, n;
    vector<double> A, B, X_opt, X, mu;
    IConstrainedOptProblem *constr_opt_problem;
    IOptProblem *opt_problem;

    if (problem.type == type_constraned::CONSTR) {
        constr_opt_problem = dynamic_cast<IConstrainedOptProblem*>(problem.problem);
        n = constr_opt_problem->GetDimension();
        constr_opt_problem->GetBounds(A, B);
        constr = constr_opt_problem->GetConstraintsNumber();
        X_opt = constr_opt_problem->GetOptimumPoint();
    } else {
        opt_problem = dynamic_cast<IOptProblem*>(problem.problem);
        n = opt_problem->GetDimension();
        opt_problem->GetBounds(A, B);
        constr = 0;
        X_opt = opt_problem->GetOptimumPoint();
    }

    mggsa.setM(constr);
    mggsa.setN(n);
    mggsa.setAB(A, B);
    mggsa.setKey(key);
    mggsa.setDen(m);
    mggsa.setIncr(incr);
    mggsa.setR(r);
    mggsa.setF(problem.f);

    mggsa.solve(count_trials, X, stop);
    mggsa.getMu(mu);

    lipschitz_const[num_func][key - key_min][m - m_min][incr - incr_min] = mu[0];
    accuracy = sqrt((X_opt[0] - X[0]) * (X_opt[0] - X[0]) + 
                    (X_opt[1] - X[1]) * (X_opt[1] - X[1]));
    count_points = mggsa.getCountPoints();

#if defined(OUTPUT_INFO)
    cout << "Function: " << problem.name << endl;
    cout << "Dimension = " << n << endl;
    cout << "Number of constrained = " << constr << endl;
    cout << "[A; B] = [(" << A[0] << ", " << A[1] << "); (" << 
                             B[0] << ", " << B[1] << ")]"<< endl;
    cout << "X* = (" << X_opt[0] << ", " << X_opt[1] << ")" << endl;
    cout << "f(X*) = " << problem.f(X_opt, constr + 1) << endl;
    cout << "Parameters for constructing the Peano curve:" << endl;
    cout << "m = " << m << " key = " << key << " incr = " << incr << endl;
    cout << "Trials result:" << endl;
    cout << "Number of trials = " << count_trials << endl;
    cout << "Number of points = " << count_points << endl;
    cout << "X = (" << X[0] << ", " << X[1] << ")" << endl;
    cout << "f(X) = " << problem.f(X, constr + 1) << endl;
    cout << "||X* - X|| = " << accuracy << endl;
    cout << "|f(X*) - f(X)| = " << abs(problem.f(X_opt, constr + 1) - problem.f(X, constr + 1)) << endl;
    cout << endl;
#else
    string str = problem.name + " key = " + to_string(key) + " m = " + to_string(m) + " incr = " + to_string(incr) +
                                " t_num = " + to_string(omp_get_thread_num()) + "\n";
    cout << str;
#endif
}