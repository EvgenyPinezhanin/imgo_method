#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS    
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <algorithm>
#include <functional>

#include <omp.h>
#include <Grishagin/GrishaginProblemFamily.hpp>
#include <Grishagin/GrishaginConstrainedProblemFamily.hpp>
#include <GKLS/GKLSProblemFamily.hpp>
#include <GKLS/GKLSConstrainedProblemFamily.hpp>
#include <mggsa.h>
#include <task.h>

using namespace std;

// #define CALC

const int type = 0; // 0 - count trials, 1 - count points, 2 - c_points / c_trials
const int family_number = 2; // 0 - Grishagin, 1 - GKLS,
                             // 2 - Grishagin(constrained), 3 - GKLS(constrained)

int main() {
#if defined(CALC)
    ofstream ofstr("output_data/mggsa_operational_characteristics_incr_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("output_data/mggsa_operational_characteristics_incr_test_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    const int number_family = 4;

    int start_time, end_time;
    double work_time;

    int num_threads_family = min(number_family, omp_get_num_procs());
    int chunk_family = 1;

    int count_func, count_trials;

    vector<double> A, B, X_opt;
    vector<vector<double>> r_array{ {3.0, 2.4, 1.6, 1.0},
                                    {4.2, 3.8, 2.0, 1.0},
                                    {3.5, 2.6, 1.8, 1.0},
                                    {4.7, 3.0, 2.0, 1.0} };
    vector<int> incr_array{0, 20, 40, 60};
    vector<double> d_array{0.0, 0.0, 0.01, 0.01};
    int den = 10, Nmax = 30000, key = 3;
    double eps = 0.01;

    vector<int> count_trials_vec, count_points_vec;

    vector<vector<vector<int>>> max_count_trials(number_family), max_count_points(number_family);
    for (int i = 0; i < number_family; i++) {
        max_count_trials[i].resize(r_array[i].size());
        max_count_points[i].resize(r_array[i].size());
        for (int j = 0; j < r_array[i].size(); j++) {
            max_count_trials[i][j].resize(incr_array.size());
            max_count_points[i][j].resize(incr_array.size());
        }
    }

    TGrishaginProblemFamily grishaginProblems;
    TGKLSProblemFamily GKLSProblems;
    TGrishaginConstrainedProblemFamily grishaginConstrainedProblems;
    TGKLSConstrainedProblemFamily GKLSConstrainedProblems;

    vector<class_problems_f> problems{ class_problems_f("GrishaginProblemFamily", &grishaginProblems, type_constraned::NONCONSTR,
                                                        "Grishagin"),
                                       class_problems_f("GKLSProblemFamily", &GKLSProblems, type_constraned::NONCONSTR, "GKLS"),
                                       class_problems_f("GrishaginProblemConstrainedFamily", &grishaginConstrainedProblems, 
                                                        type_constraned::CONSTR, "GrishaginConstrained"),
                                       class_problems_f("GKLSProblemConstrainedFamily", &GKLSConstrainedProblems, 
                                                        type_constraned::CONSTR, "GKLSConstrained") };

    mggsa_method mggsa(nullptr, -1, -1, A, B, -1.0, -1.0, den, key, eps, Nmax, -1);

#pragma omp parallel for schedule(static, chunk_family) proc_bind(spread) num_threads(num_threads_family) \
        shared(number_family, problems, d_array, r_array, incr_array, max_count_trials, max_count_points), firstprivate(mggsa) \
        private(A, B, count_func, start_time, end_time, work_time, X_opt, count_trials, count_trials_vec, count_points_vec)
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
        count_points_vec.resize(problems[i].problem->GetFamilySize());
        count_func = problems[i].problem->GetFamilySize();
        mggsa.setAB(A, B);
        mggsa.setD(d_array[i]);

        string str_input;
        for (int j = 0; j < r_array[i].size(); j++) {
            mggsa.setR(r_array[i][j]);
            for (int k = 0; k < incr_array.size(); k++) {
                mggsa.setIncr(incr_array[k]);
                start_time = clock();
                for (int l = 0; l < count_func; l++) {
                    if (problems[i].type == type_constraned::CONSTR) {
                        func_constr.current_func = l;
                        X_opt = (*func_constr.constr_opt_problem_family)[l]->GetOptimumPoint();
                        mggsa.setF(function<double(vector<double>, int)>(func_constr));
                    } else {
                        func_non_constr.current_func = l;
                        X_opt = (*func_non_constr.opt_problem_family)[l]->GetOptimumPoint();
                        mggsa.setF(function<double(vector<double>, int)>(func_non_constr));
                    }
                    if (mggsa.solve_test(X_opt, count_trials, Stop::ACCURNUMBER)) {
                        count_trials_vec[l] = count_trials;
                    } else {
                        count_trials_vec[l] = count_trials + 1;
                    }
                    count_points_vec[l] = mggsa.getCountPoints();
                }
                end_time = clock();
                work_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

                auto iter = max_element(count_trials_vec.begin(), count_trials_vec.end());
                max_count_trials[i][j][k] = *iter;
                iter = max_element(count_points_vec.begin(), count_points_vec.end());
                max_count_points[i][j][k] = *iter;

                str_input = problems[i].name + " r = " + to_string(r_array[i][j]) + " incr = " + to_string(incr_array[k]) + 
                            " count trials = " + to_string(max_count_trials[i][j][k]) + " count points = " + 
                            to_string(max_count_points[i][j][k]) + " time: " + to_string(work_time) + " t_num: " + to_string(omp_get_thread_num()) + "\n";
                cout << str_input;
            }
        }
    }

    for (int i = 0; i < number_family; i++) {
        for (int j = 0; j < r_array[i].size(); j++) {
            for (int k = 0; k < incr_array.size(); k++) {
                ofstr << incr_array[k] << " " << max_count_trials[i][j][k] << " " << max_count_points[i][j][k] << endl;
            }
            ofstr << endl << endl;
        }
        ofstr << endl << endl;
    }
    ofstr.close();

    int size = incr_array.size();
    ofstr_opt << "count_key = " << size << endl;
    ofstr_opt << "array Name[" << number_family << "]" << endl;
    ofstr_opt << "array R[" << size * number_family << "]" << endl;
    for (int i = 0; i < number_family; i++) {
        ofstr_opt << "Name[" << i + 1 << "]=\"" << problems[i].short_name << "\"" << endl;
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
    error = system("chmod +x scripts/mggsa_operational_characteristics_incr_test.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
#endif

    char str[100];
    sprintf(str, "gnuplot -c scripts/mggsa_operational_characteristics_incr_test.gp %d %d", type, family_number);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

#if defined(_MSC_VER)
    cin.get();
#endif

    return 0;
}
