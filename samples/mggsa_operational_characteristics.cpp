#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS    
#endif

#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>

#include <Grishagin/GrishaginProblemFamily.hpp>
#include <Grishagin/GrishaginConstrainedProblemFamily.hpp>
#include <GKLS/GKLSProblemFamily.hpp>
#include <GKLS/GKLSConstrainedProblemFamily.hpp>
#include <mggsa.h>
#include <task.h>

using namespace std;

#define CALC
const int family_number = 4; // 0 - Grishagin, 1 - GKLS,
                             // 2 - Grishagin(constrained), 3 - GKLS(constrained),
                             // 4 - comparison Grishagin and GKLS, 5 - comparison Grishagin and GKLS (constrained)
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
    ofstream ofstr("output_data/mggsa_operational_characteristics.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("output_data/mggsa_operational_characteristics_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    int count_func;

    int start_time, end_time;
    double work_time;

    vector<vector<int>> K{ {0, 700, 25},
                           {0, 1500, 25},
                           {0, 3000, 25},
                           {0, 4500, 25} };

    int count_successful = 0;
    int count_trials;

    int den = 10, key = 1, m = 0;
    double eps = 0.01, r = 0.0, d = 0.0;
    int Nmax = 5000;

    vector<double> A, B, X_opt;
    vector<vector<double>> r_array{ {2.5, 3.0, 3.5},
                                    {3.5, 4.3, 5.0},
                                    {2.5, 3.0, 3.5},
                                    {4.0, 4.5, 4.9} };

    vector<vector<int>> count_trials_vec(number_family);
    count_trials_vec[0].resize(grishaginProblems.GetFamilySize(), 0);
    count_trials_vec[1].resize(GKLSProblems.GetFamilySize(), 0);
    count_trials_vec[2].resize(grishaginConstrainedProblems.GetFamilySize(), 0);
    count_trials_vec[3].resize(GKLSConstrainedProblems.GetFamilySize(), 0);

    vector<class_problems_fm> problems{ class_problems_fm("GrishaginProblemFamily", &grishaginProblems, type_constraned::NONCONSTR,
                                                          f_grishagin, "Grishagin"),
                                        class_problems_fm("GKLSProblemFamily", &GKLSProblems, type_constraned::NONCONSTR, 
                                                          f_gkls, "GKLS"),
                                        class_problems_fm("GrishaginProblemConstrainedFamily", &grishaginConstrainedProblems, 
                                                          type_constraned::CONSTR, f_constr_grishagin, "GrishaginConstrained"),
                                        class_problems_fm("GKLSProblemConstrainedFamily", &GKLSConstrainedProblems, 
                                                          type_constraned::CONSTR, f_constr_gkls, "GKLSConstrained") };

    mggsa_method mggsa(nullptr, -1, -1, A, B, -1.0, d, den, key, eps, Nmax);

    IOptProblemFamily *opt_problem_family;
    IConstrainedOptProblemFamily *constr_opt_problem_family;
    for (int i = 0; i < number_family; i++) {
        ofstr_opt << "Name[" << i + 1 << "]=\"" << problems[i].short_name << "\"" << endl;
        for (int j = 0; j < r_array[i].size(); j++) {
            ofstr_opt << "R[" << (i * 3) + j + 1 << "]=\"" << r_array[i][j] << "\"" << endl; 
        }

        if (problems[i].type == type_constraned::CONSTR) {
            constr_opt_problem_family = static_cast<IConstrainedOptProblemFamily*>(problems[i].problem);
            mggsa.setN((*constr_opt_problem_family)[0]->GetDimension());
            mggsa.setM((*constr_opt_problem_family)[0]->GetConstraintsNumber());
            (*constr_opt_problem_family)[0]->GetBounds(A, B);
        } else {
            opt_problem_family = static_cast<IOptProblemFamily*>(problems[i].problem);
            mggsa.setN((*opt_problem_family)[0]->GetDimension());
            mggsa.setM(0);
            (*opt_problem_family)[0]->GetBounds(A, B);
        }
        mggsa.setAB(A, B);
        mggsa.setF(problems[i].f);
        count_func = problems[i].problem->GetFamilySize();

        cout << problems[i].name << endl;
        ofstr << "# " << problems[i].name << endl;
        for (int j = 0; j < r_array[i].size(); j++) {
            cout << "r = " << r_array[i][j] << endl;
            mggsa.setR(r_array[i][j]);
            start_time = clock();
            for (int k = 0; k < count_func; k++) {
                current_func = k;
                count_trials = K[i][1];
                if (problems[i].type == type_constraned::CONSTR) {
                    X_opt = (*constr_opt_problem_family)[k]->GetOptimumPoint();
                } else {
                    X_opt = (*opt_problem_family)[k]->GetOptimumPoint();
                }
                if (mggsa.solve_test(X_opt, count_trials, Stop::ACCURNUMBER)) {
                    count_trials_vec[i][k] = count_trials;
                } else {
                    count_trials_vec[i][k] = count_trials + 1;
                }
            }
            for (int k = K[i][0]; k <= K[i][1]; k += K[i][2]) {
                count_successful = (int)count_if(count_trials_vec[i].begin(), count_trials_vec[i].end(), [k](double elem){ return elem <= k; });
                cout << "K = " << k << " success rate = " << (double)count_successful / count_func << endl;
                ofstr << k << " " << (double)count_successful / count_func << endl;
            }
            ofstr << endl << endl;
            end_time = clock();
            work_time = ((double)end_time - start_time) / CLOCKS_PER_SEC;
            cout << "time: " << work_time << endl;
        }
    }

    ofstr.close();
    ofstr_opt.close();
#endif

    // Plotting operational characteristics(works with gnuplot)
    int error;
#if defined(__linux__)
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/mggsa_operational_characteristics.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
#endif

    char str[100];
    sprintf(str, "gnuplot -c scripts/mggsa_operational_characteristics.gp %d", family_number);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

#if defined(_MSC_VER)
    cin.get();
#endif
	return 0;
}
