#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS    
#endif

#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>

#include <omp.h>
#include <Grishagin/GrishaginProblemFamily.hpp>
#include <Grishagin/GrishaginConstrainedProblemFamily.hpp>
#include <GKLS/GKLSProblemFamily.hpp>
#include <GKLS/GKLSConstrainedProblemFamily.hpp>
#include <mggsa.h>
#include <task.h>

using namespace std;

#define CALC

const int family_number = 3; // 0 - Grishagin, 1 - GKLS,
                             // 2 - Grishagin(constrained), 3 - GKLS(constrained),
const int number_family = 4;
int current_func;

TGrishaginProblemFamily grishaginProblems;
double f_grishagin(vector<double> x, int j) {
    switch (j) {
        case 1: return grishaginProblems[current_func]->ComputeFunction(x);
        default: return numeric_limits<double>::quiet_NaN();
    }
};

TGrishaginConstrainedProblemFamily grishaginConstrainedProblems;
double f_constr_grishagin(vector<double> x, int j) {
    int constr = grishaginConstrainedProblems[current_func]->GetConstraintsNumber();
    if (j >= 1 && j <= constr) {
        return grishaginConstrainedProblems[current_func]->ComputeConstraint(j - 1, x);
    } else if (j - 1 == constr) {
        return grishaginConstrainedProblems[current_func]->ComputeFunction(x);
    } else {
        return numeric_limits<double>::quiet_NaN();
    }
};

TGKLSProblemFamily GKLSProblems;
double f_gkls(vector<double> x, int j) {
    switch (j) {
        case 1: return GKLSProblems[current_func]->ComputeFunction({ x });
        default: return numeric_limits<double>::quiet_NaN();
    }
};

TGKLSConstrainedProblemFamily GKLSConstrainedProblems;
double f_constr_gkls(vector<double> x, int j) {
    int constr = GKLSConstrainedProblems[current_func]->GetConstraintsNumber();
    if (j >= 1 && j <= constr) {
        return GKLSConstrainedProblems[current_func]->ComputeConstraint(j - 1, x);
    } else if (j - 1 == constr) {
        return GKLSConstrainedProblems[current_func]->ComputeFunction(x);
    } else {
        return numeric_limits<double>::quiet_NaN();
    }
};

int main() {
#if defined(CALC)
    ofstream ofstr("output_data/mggsa_operational_characteristics_r_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("output_data/mggsa_operational_characteristics_r_test_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    int count_func;

    vector<vector<int>> K{ {0, 700, 25},
                           {0, 1200, 25},
                           {0, 2200, 25},
                           {0, 4500, 25} };

    int count_successful = 0;
    int count_trials;

    int den = 10, m = 0;
    double eps = 0.01;
    int Nmax = 5000, incr = 30, key = 3;
    int chunk = 2;

    vector<double> A, B, X_opt;
    vector<double> r_min{1.0, 1.0, 2.2, 3.0};
    vector<double> r_max{3.0, 4.3, 3.0, 4.5};
    vector<double> step_array{0.1, 0.1, 0.05, 0.05};
    vector<double> d_array{0.0, 0.0, 0.01, 0.01};
    vector<vector<double>> r_array(number_family);
    for (int i = 0; i < number_family; i++) {
        for (double j = r_min[i]; j <= r_max[i]; j += step_array[i]) {
            r_array[i].push_back(j);
        }
    }
    vector<vector<double>> P_vector(number_family);
    for (int i = 0; i < number_family; i++) {
        P_vector[i].resize(r_array[i].size());
    }

    vector<int> count_trials_vec;

    vector<class_problems_fm> problems{ class_problems_fm("GrishaginProblemFamily", &grishaginProblems, type_constraned::NONCONSTR,
                                                          f_grishagin, "Grishagin"),
                                        class_problems_fm("GKLSProblemFamily", &GKLSProblems, type_constraned::NONCONSTR, 
                                                          f_gkls, "GKLS"),
                                        class_problems_fm("GrishaginProblemConstrainedFamily", &grishaginConstrainedProblems, 
                                                          type_constraned::CONSTR, f_constr_grishagin, "GrishaginConstrained"),
                                        class_problems_fm("GKLSProblemConstrainedFamily", &GKLSConstrainedProblems, 
                                                          type_constraned::CONSTR, f_constr_gkls, "GKLSConstrained") };

    mggsa_method mggsa(nullptr, -1, -1, A, B, -1.0, -1.0, den, key, eps, Nmax, incr);

    IOptProblemFamily *opt_problem_family;
    IConstrainedOptProblemFamily *constr_opt_problem_family;
    for (int i = 0; i < number_family; i++) {
        if (problems[i].type == type_constraned::CONSTR) {
            constr_opt_problem_family = static_cast<IConstrainedOptProblemFamily*>(problems[i].problem);
            mggsa.setN((*constr_opt_problem_family)[0]->GetDimension());
            mggsa.setM((*constr_opt_problem_family)[0]->GetConstraintsNumber());
            (*constr_opt_problem_family)[0]->GetBounds(A, B);
            count_trials_vec.resize(constr_opt_problem_family->GetFamilySize());
        } else {
            opt_problem_family = static_cast<IOptProblemFamily*>(problems[i].problem);
            mggsa.setN((*opt_problem_family)[0]->GetDimension());
            mggsa.setM(0);
            (*opt_problem_family)[0]->GetBounds(A, B);
            count_trials_vec.resize(opt_problem_family->GetFamilySize());
        }
        mggsa.setAB(A, B);
        mggsa.setF(problems[i].f);
        mggsa.setD(d_array[i]);
        count_func = problems[i].problem->GetFamilySize();

        cout << problems[i].name << endl;
    #pragma omp parallel for schedule(dynamic, chunk) proc_bind(spread) num_threads(omp_get_num_procs()) \
            shared(P_vector, r_array, count_func, current_func, K, problems, constr_opt_problem_family, opt_problem_family) \
            firstprivate(count_trials_vec, mggsa) private(count_trials, X_opt, count_successful)
        for (int j = 0; j < r_array[i].size(); j++) {
            mggsa.setR(r_array[i][j]);
            for (int k = 0; k < count_func; k++) {
                current_func = k;
                count_trials = K[i][1];
                if (problems[i].type == type_constraned::CONSTR) {
                    X_opt = (*constr_opt_problem_family)[k]->GetOptimumPoint();
                } else {
                    X_opt = (*opt_problem_family)[k]->GetOptimumPoint();
                }
                if (mggsa.solve_test(X_opt, count_trials, Stop::ACCURNUMBER)) {
                    count_trials_vec[k] = count_trials;
                } else {
                    count_trials_vec[k] = count_trials + 1;
                }
            }

            int k = K[i][1];
            count_successful = (int)count_if(count_trials_vec.begin(), count_trials_vec.end(), [k](double elem){ return elem <= k; });
            P_vector[i][j] = (double)count_successful / count_func;

            string str = "r = " + to_string(r_array[i][j]) + " t_num = " + to_string(omp_get_thread_num()) + "\n";
            cout << str;
        }
    }

    for (int i = 0; i < number_family; i++) {
        for (int j = 0; j < P_vector[i].size(); j++) {
            ofstr << r_array[i][j] << " " << P_vector[i][j] << endl;
        }
        ofstr << endl << endl;
    }
    ofstr.close();

    for (int i = 0; i < number_family; i++) {
        ofstr_opt << "Name[" << i + 1 << "]=\"" << problems[i].short_name << "\"" << endl;
    }
    ofstr_opt.close();
#endif

    // Plotting operational characteristics(works with gnuplot)
    int error;
#if defined(__linux__)
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/mggsa_operational_characteristics_r_test.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
#endif

    char str[100];
    sprintf(str, "gnuplot -c scripts/mggsa_operational_characteristics_r_test.gp %d", family_number);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

#if defined(_MSC_VER)
    cin.get();
#endif
	return 0;
}
