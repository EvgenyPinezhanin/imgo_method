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
const int family_number = 3; // 0 - Grishagin, 1 - GKLS,
                             // 2 - Grishagin(constrained), 3 - GKLS(constrained)

int main() {
    ofstream ofstr("output_data/mggsa_operational_characteristics_r_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("output_data/mggsa_operational_characteristics_r_test_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    const int number_family = 4;

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

    vector<int> K{ 700, 1200, 2500, 4500 };

    vector<double> r_min{1.0, 1.0, 1.0, 1.0};
    vector<double> r_max{5.0, 5.0, 5.0, 5.0};
    vector<double> step_array{0.05, 0.05, 0.01, 0.01};
    vector<vector<double>> r_array(number_family);
    vector<int> r_size(number_family);
    for (int i = 0; i < number_family; i++) {
        for (double j = r_min[i]; j <= r_max[i]; j += step_array[i]) {
            r_array[i].push_back(j);
        }
        r_size[i] = r_array[i].size();
    }
    vector<int> key_array{1, 3};

    int number_functions;
    vector<vector<vector<vector<int>>>> trials_data(number_family), points_data(number_family);
    for (int i = 0; i < number_family; i++) {
        number_functions = problems[i].problem->GetFamilySize();
        trials_data[i].resize(key_array.size());
        points_data[i].resize(key_array.size());
        for (int j = 0; j < key_array.size(); j++) {
            trials_data[i][j].resize(r_size[i]);
            points_data[i][j].resize(r_size[i]);
            for (int k = 0; k < r_size[i]; k++) {
                trials_data[i][j][k].resize(number_functions);
                points_data[i][j][k].resize(number_functions);
            }
        }
    }

    vector<vector<vector<double>>> P_vector(number_family);
    vector<vector<vector<int>>> max_count_trials(number_family), max_count_points(number_family);
    for (int i = 0; i < number_family; i++) {
        P_vector[i].resize(key_array.size());
        max_count_trials[i].resize(key_array.size());
        max_count_points[i].resize(key_array.size());
        for (int j = 0; j < key_array.size(); j++) {
            P_vector[i][j].resize(r_array[i].size());
            max_count_trials[i][j].resize(r_array[i].size());
            max_count_points[i][j].resize(r_array[i].size());
        }
    }

#if defined(CALC)
    ofstream ofstr_data("output_data/mggsa_operational_characteristics_r_test_data.txt");
    if (!ofstr_data.is_open()) cerr << "File opening error\n";

    int num_threads = min(number_family * (int)key_array.size(), omp_get_num_procs());
    int chunk = 10;

    int den = 10, incr = 20;
    double eps = 0.01;
    vector<double> d_array{0.0, 0.0, 0.01, 0.01};
    vector<double> Nmax_array{10000, 15000, 25000, 30000};

    mggsa_method mggsa(nullptr, -1, -1, vector<double>{}, vector<double>{}, -1.0, -1.0, den, -1, eps, -1, incr);

    int count_r = 0;
    for (int i = 0; i < r_size.size(); i++) {
        count_r += r_size[i];
    }

#pragma omp parallel for schedule(static, chunk) proc_bind(spread) num_threads(omp_get_num_procs()) collapse(2) \
        shared(count_r, d_array, Nmax_array, r_array, P_vector, trials_data, points_data) \
        firstprivate(K, r_size, key_array, problems, mggsa)
    for (int t = 0; t < count_r; t++) {
        for (int j = 0; j < key_array.size(); j++) {
            int i, k, t_tmp = t;
            for (int l = 0; l < r_size.size(); l++) {
                if (t_tmp < r_size[l]) {
                    i = l; k = t_tmp;
                    break;
                } else {
                    t_tmp -= r_size[l];
                }
            }

            vector<double> A, B;
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

            int number_functions = problems[i].problem->GetFamilySize();
            vector<int> count_trials_vec(number_functions), count_points_vec(number_functions);
            mggsa.setAB(A, B);
            mggsa.setR(r_array[i][k]);
            mggsa.setD(d_array[i]);
            mggsa.setKey(key_array[j]);
            mggsa.setNmax(Nmax_array[i]);

            vector<double> X_opt;
            int count_trials;
            int start_time = clock();
            for (int l = 0; l < number_functions; l++) {
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
            trials_data[i][j][k] = count_trials_vec;
            points_data[i][j][k] = count_points_vec;
            int end_time = clock();
            double work_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

            string str_input = problems[i].name + " key = " + to_string(key_array[j]) + " r = " + to_string(r_array[i][k]) + 
                               " time: " + to_string(work_time) + " t_num: " + to_string(omp_get_thread_num()) + "\n";
            cout << str_input;
        }
    }

    for (int i = 0; i < number_family; i++) {
        for (int j = 0; j < key_array.size(); j++) {
            for (int k = 0; k < r_array[i].size(); k++) {
                for (int l = 0; l < trials_data[i][j][k].size(); l++) {
                    ofstr_data << trials_data[i][j][k][l] << " ";
                }
                ofstr_data << endl;
                for (int l = 0; l < points_data[i][j][k].size(); l++) {
                    ofstr_data << points_data[i][j][k][l] << " ";
                }
                ofstr_data << endl;
            }
        }
    }
    ofstr_data.close();
#else
    ifstream ifstr_data("output_data/mggsa_operational_characteristics_r_test_data.txt");
    if (!ifstr_data.is_open()) cerr << "File opening error\n";

    size_t size;
    for (int i = 0; i < number_family; i++) {
        for (int j = 0; j < key_array.size(); j++) {
            for (int k = 0; k < r_array[i].size(); k++) {
                size = trials_data[i][j][k].size();
                for (int l = 0; l < size; l++) {
                    ifstr_data >> trials_data[i][j][k][l];
                }
                for (int l = 0; l < size; l++) {
                    ifstr_data >> points_data[i][j][k][l];
                }
            }
        }
    }
    ifstr_data.close();
#endif

    int count_successful, h, index;
    for (int i = 0; i < number_family; i++) {
        number_functions = problems[i].problem->GetFamilySize();
        h = K[i];
        for (int j = 0; j < key_array.size(); j++) {
            for (int k = 0; k < r_array[i].size(); k++) {
                count_successful = (int)count_if(trials_data[i][j][k].begin(), trials_data[i][j][k].end(), [h](double elem){ return elem <= h; });
                P_vector[i][j][k] = (double)count_successful / number_functions;

                auto iter = max_element(trials_data[i][j][k].begin(), trials_data[i][j][k].end());
                max_count_trials[i][j][k] = *iter;
                for (int l = 0; l < trials_data[i][j][k].size(); l++) {
                    if (*iter == trials_data[i][j][k][l]) index = l;
                }
                max_count_points[i][j][k] = points_data[i][j][k][index];
            }
        }
    }

    for (int i = 0; i < number_family; i++) {
        for (int j = 0; j < key_array.size(); j++) {
            for (int k = 0; k < r_array[i].size(); k++) {
                ofstr << r_array[i][k] << " " << P_vector[i][j][k] << " " <<
                         max_count_trials[i][j][k] << " " << max_count_points[i][j][k] << endl;
            }
            ofstr << endl << endl;
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
