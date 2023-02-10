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

// #define CALC

const int family_number = 3; // 0 - Grishagin, 1 - GKLS,
                             // 2 - Grishagin(constrained), 3 - GKLS(constrained),

int main() {
#if defined(CALC)
    ofstream ofstr("output_data/mggsa_operational_characteristics_key_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("output_data/mggsa_operational_characteristics_key_test_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    int count_func, count_successful, count_trials;

    int start_time, end_time;
    double work_time;

    vector<vector<int>> K{ {0, 700, 25},
                           {0, 1200, 25},
                           {0, 2200, 25},
                           {0, 4500, 25} };

    int den = 10, m = 0, Nmax = 5000, incr = 30;
    double eps = 0.01;

    vector<double> A, B, X_opt;
    vector<vector<double>> r_array{ {3.0, 2.4, 1.6, 1.0},
                                    {4.2, 3.8, 2.0, 1.0},
                                    {3.5, 2.6, 1.8, 1.0},
                                    {4.7, 3.0, 2.0, 1.0} };
    vector<int> key_array{1, 3, 3, 3};
    vector<double> d_array{0.0, 0.0, 0.01, 0.01};

    const int number_family = 4;

    TGrishaginProblemFamily grishaginProblems;
    TGrishaginConstrainedProblemFamily grishaginConstrainedProblems;
    TGKLSProblemFamily GKLSProblems;
    TGKLSConstrainedProblemFamily GKLSConstrainedProblems;

    vector<vector<int>> count_trials_vec(number_family);
    count_trials_vec[0].resize(grishaginProblems.GetFamilySize(), 0);
    count_trials_vec[1].resize(GKLSProblems.GetFamilySize(), 0);
    count_trials_vec[2].resize(grishaginConstrainedProblems.GetFamilySize(), 0);
    count_trials_vec[3].resize(GKLSConstrainedProblems.GetFamilySize(), 0);

    vector<class_problems_f> problems{ class_problems_f("GrishaginProblemFamily", &grishaginProblems, type_constraned::NONCONSTR,
                                                        "Grishagin"),
                                       class_problems_f("GKLSProblemFamily", &GKLSProblems, type_constraned::NONCONSTR, "GKLS"),
                                       class_problems_f("GrishaginProblemConstrainedFamily", &grishaginConstrainedProblems, 
                                                          type_constraned::CONSTR, "GrishaginConstrained"),
                                       class_problems_f("GKLSProblemConstrainedFamily", &GKLSConstrainedProblems, 
                                                          type_constraned::CONSTR, "GKLSConstrained") };

    mggsa_method mggsa(nullptr, -1, -1, A, B, -1.0, -1.0, den, -1, eps, Nmax, incr);

    for (int i = 0; i < number_family; i++) {
        functor_non_constr func_non_constr;
        functor_constr func_constr;
        if (problems[i].type == type_constraned::CONSTR) {
            func_constr.constr_opt_problem_family = static_cast<IConstrainedOptProblemFamily*>(problems[i].problem);
            (*func_constr.constr_opt_problem_family)[0]->GetBounds(A, B);
            mggsa.setN((*func_constr.constr_opt_problem_family)[0]->GetDimension());
            mggsa.setM((*func_constr.constr_opt_problem_family)[0]->GetConstraintsNumber());
        } else {
            func_non_constr.opt_problem_family = static_cast<IOptProblemFamily*>(problems[i].problem);
            (*func_non_constr.opt_problem_family)[0]->GetBounds(A, B);
            mggsa.setN((*func_non_constr.opt_problem_family)[0]->GetDimension());
            mggsa.setM(0);
        }

        count_trials_vec.resize(problems[i].problem->GetFamilySize());
        count_func = problems[i].problem->GetFamilySize();
        mggsa.setAB(A, B);
        mggsa.setD(d_array[i]);

        cout << problems[i].name << endl;
        ofstr << "# " << problems[i].name << endl;
        for (int j = 0; j < r_array[i].size(); j++) {
            cout << "r = " << r_array[i][j] << endl;
            mggsa.setKey(key_array[j]);
            mggsa.setR(r_array[i][j]);
            start_time = clock();
            for (int k = 0; k < count_func; k++) {
                count_trials = K[i][1];
                if (problems[i].type == type_constraned::CONSTR) {
                    func_constr.current_func = k;
                    X_opt = (*func_constr.constr_opt_problem_family)[k]->GetOptimumPoint();
                    mggsa.setF(function<double(vector<double>, int)>(func_constr));
                } else {
                    func_non_constr.current_func = k;
                    X_opt = (*func_non_constr.opt_problem_family)[k]->GetOptimumPoint();
                    mggsa.setF(function<double(vector<double>, int)>(func_non_constr));
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

    int size = key_array.size();
    ofstr_opt << "count_key = " << size << endl;
    ofstr_opt << "array Name[" << number_family << "]" << endl;
    ofstr_opt << "array R[" << size * number_family << "]" << endl;
    ofstr_opt << "array Key[" << size << "]" << endl;
    for (int i = 0; i < number_family; i++) {
        ofstr_opt << "Name[" << i + 1 << "]=\"" << problems[i].short_name << "\"" << endl;
        ofstr_opt << "Key[" << i + 1 << "]=\"" << key_array[i] << "\"" << endl;
        for (int j = 0; j < r_array[i].size(); j++) {
            ofstr_opt << "R[" << (i * size) + j + 1 << "]=\"" << r_array[i][j] << "\"" << endl; 
        }
    }
    ofstr_opt.close();
#endif

    // Plotting operational characteristics(works with gnuplot)
    int error;
#if defined(__linux__)
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/mggsa_operational_characteristics_key_test.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
#endif

    char str[100];
    sprintf(str, "gnuplot -c scripts/mggsa_operational_characteristics_key_test.gp %d", family_number);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

#if defined(_MSC_VER)
    cin.get();
#endif

	return 0;
}
