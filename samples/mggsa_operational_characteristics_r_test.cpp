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

const int type = 1; // 0 - P_max, 1 - count_trials, 2 - count points, 3 - c_points / c_trials
const int family_number = 2; // 0 - Grishagin, 1 - GKLS,
                             // 2 - Grishagin(constrained), 3 - GKLS(constrained)

int main() {
    ofstream ofstr("output_data/mggsa_operational_characteristics_r_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("output_data/mggsa_operational_characteristics_r_test_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    const int number_family = 4;

    int count_func;

    vector<int> K{ 700, 1200, 2500, 4500 };

    vector<double> r_min{1.0, 1.0, 1.0, 1.0};
    vector<double> r_max{5.0, 5.0, 5.0, 5.0};
    vector<double> step_array{0.05, 0.05, 0.01, 0.01};
    vector<vector<double>> r_array(number_family);
    for (int i = 0; i < number_family; i++) {
        for (double j = r_min[i]; j <= r_max[i]; j += step_array[i]) {
            r_array[i].push_back(j);
        }
    }

    vector<vector<vector<int>>> trials_data(number_family), points_data(number_family);
    for (int i = 0; i < number_family; i++) {
        trials_data[i].resize(r_array[i].size());
        points_data[i].resize(r_array[i].size());
    }

    vector<vector<double>> P_vector(number_family);
    vector<vector<int>> max_count_trials(number_family), max_count_points(number_family);
    for (int i = 0; i < number_family; i++) {
        P_vector[i].resize(r_array[i].size());
        max_count_trials[i].resize(r_array[i].size());
        max_count_points[i].resize(r_array[i].size());
    }

    TGrishaginProblemFamily grishaginProblems;
    TGrishaginConstrainedProblemFamily grishaginConstrainedProblems;
    TGKLSProblemFamily GKLSProblems;
    TGKLSConstrainedProblemFamily GKLSConstrainedProblems;

    vector<class_problems_f> problems{ class_problems_f("GrishaginProblemFamily", &grishaginProblems, type_constraned::NONCONSTR,
                                                        "Grishagin"),
                                       class_problems_f("GKLSProblemFamily", &GKLSProblems, type_constraned::NONCONSTR, "GKLS"),
                                       class_problems_f("GrishaginProblemConstrainedFamily", &grishaginConstrainedProblems,
                                                        type_constraned::CONSTR, "GrishaginConstrained"),
                                       class_problems_f("GKLSProblemConstrainedFamily", &GKLSConstrainedProblems,
                                                        type_constraned::CONSTR, "GKLSConstrained") };

#if defined(CALC)
    ofstream ofstr_data("output_data/mggsa_operational_characteristics_r_test_data.txt");
    if (!ofstr_data.is_open()) cerr << "File opening error\n";

    int num_threads_family = min(number_family, omp_get_num_procs());
    int chunk_family = 1;
    int num_threads_r = max(1, (omp_get_num_procs() - num_threads_family) / num_threads_family + 1);
    int chunk_r = 10;

    int den = 10, key = 3, Nmax = 30000, incr = 30;
    double eps = 0.01;
    vector<double> d_array{0.0, 0.0, 0.01, 0.01};

    mggsa_method mggsa(nullptr, -1, -1, vector<double>{}, vector<double>{}, -1.0, -1.0, den, key, eps, Nmax, incr);

    omp_set_nested(1);
#pragma omp parallel for schedule(static, chunk_family) proc_bind(spread) num_threads(num_threads_family) \
        shared(number_family, d_array, r_array, P_vector, trials_data, points_data) \
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

        count_func = problems[i].problem->GetFamilySize();
        size_t size = r_array[i].size();
        for (int j = 0; j < size; j++) {
            trials_data[i][j].resize(count_func);
            points_data[i][j].resize(count_func);
        }
        count_trials_vec.resize(count_func);
        count_points_vec.resize(count_func);
        mggsa.setAB(A, B);
        mggsa.setD(d_array[i]);

        thread_num = omp_get_thread_num();

    #pragma omp parallel for schedule(dynamic, chunk_r) proc_bind(close) num_threads(num_threads_r) \
            shared(d_array, r_array, P_vector, max_count_trials, max_count_points) \
            firstprivate(K, problems, mggsa, count_trials_vec, count_points_vec, count_func, func_constr, func_non_constr, thread_num)
        for (int j = 0; j < r_array[i].size(); j++) {
            vector<double> X_opt;
            int count_trials;
            int start_time, end_time;
            double work_time;
            string str_input;

            mggsa.setR(r_array[i][j]);
            start_time = clock();
            for (int k = 0; k < count_func; k++) {
                count_trials = K[i];
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
            for (int k = 0; k < count_func; k++) {
                trials_data[i][j][k] = count_trials_vec[k];
                points_data[i][j][k] = count_points_vec[k];
            }
            end_time = clock();
            work_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

            str_input = problems[i].name + " r = " + to_string(r_array[i][j]) + " time: " + to_string(work_time) +
                        " t_num_ext: " + to_string(thread_num) + " t_num_int: " + to_string(omp_get_thread_num()) + "\n";
            cout << str_input;
        }
    }

    for (int i = 0; i < number_family; i++) {
        for (int j = 0; j < r_array[i].size(); j++) {
            for (int k = 0; k < trials_data[i][j].size(); k++) {
                ofstr_data << trials_data[i][j][k] << " ";
            }
            ofstr_data << endl;
            for (int k = 0; k < points_data[i][j].size(); k++) {
                ofstr_data << points_data[i][j][k] << " ";
            }
            ofstr_data << endl;
        }
    }
    ofstr_data.close();
#else
    ifstream ifstr_data("output_data/mggsa_operational_characteristics_r_test_data.txt");
    if (!ifstr_data.is_open()) cerr << "File opening error\n";

    size_t size;
    for (int i = 0; i < number_family; i++) {
        count_func = problems[i].problem->GetFamilySize();
        size = r_array[i].size();
        for (int j = 0; j < size; j++) {
            trials_data[i][j].resize(count_func);
            points_data[i][j].resize(count_func);

            for (int k = 0; k < count_func; k++) {
                ifstr_data >> trials_data[i][j][k];
            }
            for (int k = 0; k < count_func; k++) {
                ifstr_data >> points_data[i][j][k];
            }
        }
    }
    ifstr_data.close();
#endif

    int count_successful, k, index;

    for (int i = 0; i < number_family; i++) {
        k = K[i];
        count_func = problems[i].problem->GetFamilySize();
        for (int j = 0; j < r_array[i].size(); j++) {
            count_successful = (int)count_if(trials_data[i][j].begin(), trials_data[i][j].end(), [k](double elem){ return elem <= k; });
            P_vector[i][j] = (double)count_successful / count_func;

            auto iter = max_element(trials_data[i][j].begin(), trials_data[i][j].end());
            max_count_trials[i][j] = *iter;
            for (int k = 0; k < trials_data[i][j].size(); k++) {
                if (*iter == trials_data[i][j][k]) index = k;
            }
            max_count_points[i][j] = points_data[i][j][index];
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

    ofstr_opt << "array Name[" << number_family << "]" << endl;
    ofstr_opt << "array K_max[" << number_family << "]" << endl;
    for (int i = 0; i < number_family; i++) {
        ofstr_opt << "Name[" << i + 1 << "]=\"" << problems[i].short_name << "\"" << endl;
        ofstr_opt << "K_max[" << i + 1 << "]= " << K[i] << endl;
    }
    ofstr_opt.close();

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
