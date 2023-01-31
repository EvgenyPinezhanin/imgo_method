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

#define CALC

const int type = 3; // 0 - P_max, 1 - count_trials, 2 - count points, 3 - c_points / c_trials
const int family_number = 0; // 0 - Grishagin, 1 - GKLS,
                             // 2 - Grishagin(constrained), 3 - GKLS(constrained)

int main() {
#if defined(CALC)
    ofstream ofstr("output_data/mggsa_operational_characteristics_r_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("output_data/mggsa_operational_characteristics_r_test_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    const int number_family = 4;

    int num_threads_family = min(number_family, omp_get_num_procs());
    int chunk_family = 1;
    int num_threads_r = max(1, (omp_get_num_procs() - num_threads_family) / num_threads_family + 1);
    int chunk_r = 10;

    vector<vector<int>> K{ {0, 700, 10},
                           {0, 1200, 10},
                           {0, 2200, 10},
                           {0, 4500, 10} };

    TGrishaginProblemFamily grishaginProblems;
    TGrishaginConstrainedProblemFamily grishaginConstrainedProblems;
    TGKLSProblemFamily GKLSProblems;
    TGKLSConstrainedProblemFamily GKLSConstrainedProblems;

    int den = 10, key = 3, Nmax = 5000, incr = 30;
    double eps = 0.01;
    vector<double> r_min{1.0, 1.0, 1.0, 1.0};
    vector<double> r_max{5.0, 5.0, 5.0, 5.0};
    vector<double> step_array{0.05, 0.05, 0.01, 0.01};
    vector<double> d_array{0.0, 0.0, 0.01, 0.01};
    vector<vector<double>> r_array(number_family);
    for (int i = 0; i < number_family; i++) {
        for (double j = r_min[i]; j <= r_max[i]; j += step_array[i]) {
            r_array[i].push_back(j);
        }
    }

    vector<vector<double>> P_vector(number_family);
    vector<vector<int>> max_count_trials(number_family), max_count_points(number_family);
    for (int i = 0; i < number_family; i++) {
        P_vector[i].resize(r_array[i].size());
        max_count_trials[i].resize(r_array[i].size());
        max_count_points[i].resize(r_array[i].size());
    }

    vector<class_problems_f> problems{ class_problems_f("GrishaginProblemFamily", &grishaginProblems, type_constraned::NONCONSTR,
                                                        "Grishagin"),
                                       class_problems_f("GKLSProblemFamily", &GKLSProblems, type_constraned::NONCONSTR, "GKLS"),
                                       class_problems_f("GrishaginProblemConstrainedFamily", &grishaginConstrainedProblems, 
                                                        type_constraned::CONSTR, "GrishaginConstrained"),
                                       class_problems_f("GKLSProblemConstrainedFamily", &GKLSConstrainedProblems, 
                                                        type_constraned::CONSTR, "GKLSConstrained") };

    mggsa_method mggsa(nullptr, -1, -1, vector<double>{}, vector<double>{}, -1.0, -1.0, den, key, eps, Nmax, incr);

    omp_set_nested(1);
#pragma omp parallel for schedule(static, chunk_family) proc_bind(spread) num_threads(num_threads_family) \
        shared(number_family, d_array, r_array, P_vector, max_count_trials, max_count_points) \
        firstprivate(K, problems, mggsa)
    for (int i = 0; i < number_family; i++) {
        vector<double> A, B;
        vector<int> count_trials_vec, count_points_vec;
        int count_func, thread_num;
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

        thread_num = omp_get_thread_num();

    #pragma omp parallel for schedule(dynamic, chunk_r) proc_bind(close) num_threads(num_threads_r) \
            shared(d_array, r_array, P_vector, max_count_trials, max_count_points) \
            firstprivate(K, problems, mggsa, count_trials_vec, count_points_vec, count_func, func_constr, func_non_constr, thread_num)
        for (int j = 0; j < r_array[i].size(); j++) {
            vector<double> X_opt;
            int count_successful, count_trials;
            int start_time, end_time;
            double work_time;
            string str_input;

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
                    count_trials_vec[k] = count_trials;
                } else {
                    count_trials_vec[k] = count_trials + 1;
                }
                count_points_vec[k] = mggsa.getCountPoints();
            }
            end_time = clock();
            work_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

            int k = K[i][1];
            count_successful = (int)count_if(count_trials_vec.begin(), count_trials_vec.end(), [k](double elem){ return elem <= k; });
            P_vector[i][j] = (double)count_successful / count_func;

            auto iter = max_element(count_trials_vec.begin(), count_trials_vec.end());
            max_count_trials[i][j] = *iter;
            int index;
            for (int k = 0; k < count_trials_vec.size(); k++) {
                if (*iter == count_trials_vec[k]) index = k;
            }
            max_count_points[i][j] = count_points_vec[index];

            str_input = problems[i].name + " r = " + to_string(r_array[i][j]) + " count trials = " + to_string(max_count_trials[i][j]) + 
                        " count points = " + to_string(max_count_points[i][j]) + " time: " + to_string(work_time) + " t_num_ext: " +
                        to_string(thread_num) + " t_num_int: " + to_string(omp_get_thread_num()) + "\n";
            cout << str_input;
        }
    }

    for (int i = 0; i < number_family; i++) {
        for (int j = 0; j < P_vector[i].size(); j++) {
            ofstr << r_array[i][j] << " " << P_vector[i][j] << " " << 
                     max_count_trials[i][j] << " " << max_count_points[i][j] << endl;
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
    sprintf(str, "gnuplot -c scripts/mggsa_operational_characteristics_r_test.gp %d %d", type, family_number);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

#if defined(_MSC_VER)
    cin.get();
#endif

	return 0;
}
