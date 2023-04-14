#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS
    #define PROC_BIND
#else
    #define PROC_BIND proc_bind(spread)
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

const int number_family = 4;

const int chunk = 4;

int main() {
    ofstream ofstr_opt("output_data/mggsa_operational_characteristics_incr_test_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    vector<vector<double>> r_array{ {3.0, 2.4, 1.6, 1.0},
                                    {4.2, 3.8, 2.0, 1.0},
                                    {3.5, 2.6, 1.8, 1.0},
                                    {4.7, 3.0, 2.0, 1.0} };
    vector<int> incr_array{0, 20, 40, 60};
    vector<double> d_array{0.0, 0.0, 0.01, 0.01};

    TGrishaginProblemFamily grishaginProblems;
    TGKLSProblemFamily GKLSProblems;
    TGrishaginConstrainedProblemFamily grishaginConstrainedProblems;
    TGKLSConstrainedProblemFamily GKLSConstrainedProblems;

    vector<problem_family> problems{ problem_family("GrishaginProblemFamily", &grishaginProblems, type_constraned::NONCONSTR,
                                                    "Grishagin"),
                                     problem_family("GKLSProblemFamily", &GKLSProblems, type_constraned::NONCONSTR, "GKLS"),
                                     problem_family("GrishaginProblemConstrainedFamily", &grishaginConstrainedProblems, 
                                                    type_constraned::CONSTR, "GrishaginConstrained"),
                                     problem_family("GKLSProblemConstrainedFamily", &GKLSConstrainedProblems, 
                                                    type_constraned::CONSTR, "GKLSConstrained") };

#if defined(CALC)
    ofstream ofstr("output_data/mggsa_operational_characteristics_incr_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";

    int den = 10, key = 3;
    double eps = 0.01;
    int maxIters = 30000, maxEvals = 100000;

    vector<vector<vector<int>>> maxNumberIters(number_family), maxNumberTrialPoints(number_family);
    for (int i = 0; i < number_family; i++) {
        maxNumberIters[i].resize(r_array[i].size());
        maxNumberTrialPoints[i].resize(r_array[i].size());
        for (int j = 0; j < r_array[i].size(); j++) {
            maxNumberIters[i][j].resize(incr_array.size());
            maxNumberTrialPoints[i][j].resize(incr_array.size());
        }
    }

    mggsa_method mggsa(nullptr, -1, -1, vector<double>{}, vector<double>{}, -1.0, -1.0, den, key, eps, maxIters, maxEvals, -1);

    int r_size = r_array[0].size();

    int total_start_time = clock();
#pragma omp parallel for schedule(static, chunk) PROC_BIND num_threads(omp_get_num_procs()) collapse(2) \
        shared(number_family, problems, d_array, r_array, r_size, incr_array, maxNumberIters, maxNumberTrialPoints) \
        firstprivate(mggsa)
    for (int i = 0; i < number_family; i++) {
        for (int j = 0; j < r_size; j++) {
            vector<double> A, B;
            functor_family func;
            functor_family_constr func_constr;

            if (problems[i].type == type_constraned::CONSTR) {
                func_constr.constr_opt_problem_family = static_cast<IConstrainedOptProblemFamily*>(problems[i].optProblemFamily);
                (*func_constr.constr_opt_problem_family)[0]->GetBounds(A, B);
                mggsa.setN((*func_constr.constr_opt_problem_family)[0]->GetDimension());
                mggsa.setNumberConstraints((*func_constr.constr_opt_problem_family)[0]->GetConstraintsNumber());
            } else {
                func.opt_problem_family = static_cast<IOptProblemFamily*>(problems[i].optProblemFamily);
                (*func.opt_problem_family)[0]->GetBounds(A, B);
                mggsa.setN((*func.opt_problem_family)[0]->GetDimension());
                mggsa.setNumberConstraints(0);
            }

            vector<int> count_iters_vec(problems[i].optProblemFamily->GetFamilySize());
            vector<int> count_trials_vec(problems[i].optProblemFamily->GetFamilySize());
            int count_func = problems[i].optProblemFamily->GetFamilySize();
            mggsa.setAB(A, B);
            mggsa.setD(d_array[i]);
            mggsa.setR(r_array[i][j]);

            vector<double> X_opt; 
            int index;
            int countIters, countEvals;
            string str_input;
            double start_time, end_time, work_time;

            for (int k = 0; k < incr_array.size(); k++) {
                mggsa.setIncr(incr_array[k]);

                start_time = omp_get_wtime();
                for (int l = 0; l < count_func; l++) {
                    if (problems[i].type == type_constraned::CONSTR) {
                        func_constr.current_func = l;
                        X_opt = (*func_constr.constr_opt_problem_family)[l]->GetOptimumPoint();
                        mggsa.setF(func_constr);
                    } else {
                        func.current_func = l;
                        X_opt = (*func.opt_problem_family)[l]->GetOptimumPoint();
                        mggsa.setF(func);
                    }
                    if (mggsa.solve_test(X_opt, countIters, countEvals)) {
                        count_iters_vec[l] = countIters;
                    } else {
                        count_iters_vec[l] = countIters + 1;
                    }
                    count_trials_vec[l] = mggsa.getNumberTrialPoints();
                }
                end_time = omp_get_wtime();
                work_time = end_time - start_time;

                auto iter = max_element(count_iters_vec.begin(), count_iters_vec.end());
                maxNumberIters[i][j][k] = *iter;
                for (int l = 0; l < count_iters_vec.size(); l++) {
                    if (*iter == count_iters_vec[l]) index = l;
                }
                maxNumberTrialPoints[i][j][k] = count_trials_vec[index];

                str_input = problems[i].name + " r = " + to_string(r_array[i][j]) + " incr = " + to_string(incr_array[k]) + 
                            " count iters = " + to_string(maxNumberIters[i][j][k]) + " count trials = " + 
                            to_string(maxNumberTrialPoints[i][j][k]) + " time: " + to_string(work_time) + " t_num: " +
                            to_string(omp_get_thread_num()) + "\n";
                cout << str_input;
            }
        }
    }
    int total_end_time = clock();
    double total_work_time = ((double)total_end_time - total_start_time) / CLOCKS_PER_SEC;
    cout << "Total time: " << total_work_time << endl;

    for (int i = 0; i < number_family; i++) {
        for (int j = 0; j < r_array[i].size(); j++) {
            for (int k = 0; k < incr_array.size(); k++) {
                ofstr << incr_array[k] << " " << maxNumberIters[i][j][k] << " " << maxNumberTrialPoints[i][j][k] << endl;
            }
            ofstr << endl << endl;
        }
        ofstr << endl << endl;
    }
    ofstr.close();
#endif

    size_t size = incr_array.size();
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
