#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS
    #define PROC_BIND
#else
    #define PROC_BIND proc_bind(spread)
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
#include <direct_method.h>
#include <mggsa.h>
#include <task.h>

using namespace std;

// #define CALC_DIRECT
// #define CALC_MGGSA

const int family_number = 3; // 0 - Grishagin, 1 - GKLS,
                             // 2 - Grishagin(constrained), 3 - GKLS(constrained),

double f(int n, const double *X, int *undefined_flag, void *data) {
    data_direct_oper_character *f_data = static_cast<data_direct_oper_character*>(data);
    f_data->count_evals++;
    vector<double> point(n), opt_point(n);
    for (int i = 0; i < n; i++) {
        point[i] = X[i];
    }
    f_data->points.push_back(point);

    double f;
    if (f_data->type == type_constraned::CONSTR) {
        functor_family_constr *problem = static_cast<functor_family_constr*>(f_data->functor);

        opt_point = (*problem->constr_opt_problem_family)[problem->current_func]->GetOptimumPoint();

        bool is_constr = true;
        int constr = (*problem->constr_opt_problem_family)[problem->current_func]->GetConstraintsNumber();
        for (int i = 0; i < constr; i++) {
            if ((*problem)(point, i) > 0.0) is_constr = false;
        }
        if (is_constr) *undefined_flag = 1;

        f = (*problem)(point, constr + 1);
    } else {
        functor_family *problem = static_cast<functor_family*>(f_data->functor);

        opt_point = (*problem->opt_problem_family)[problem->current_func]->GetOptimumPoint();

        f = (*problem)(point, 1);
    }

    if (!f_data->converge) {
        double dist = 0.0;
        size_t size = point.size();
        for (int i = 0; i < size; i++) {
            dist += (point[i] - opt_point[i]) * (point[i] - opt_point[i]);
        }
        if (dist <= f_data->eps) {
            f_data->converge = true;
            f_data->min_count_evals = f_data->count_evals;
        }
    }

    return f;
}

int main() {
    ofstream ofstr_opt("output_data/direct_operational_characteristics_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";
    
    const int number_family = 4;

    TGrishaginProblemFamily grishaginProblems;
    TGrishaginConstrainedProblemFamily grishaginConstrainedProblems;
    TGKLSProblemFamily gklsProblems;
    TGKLSConstrainedProblemFamily gklsConstrainedProblems;

    vector<problem_family> problems{ problem_family("GrishaginProblemFamily", &grishaginProblems, type_constraned::NONCONSTR, "Grishagin"),
                                     problem_family("GKLSProblemFamily", &gklsProblems, type_constraned::NONCONSTR, "GKLS"),
                                     problem_family("GrishaginProblemConstrainedFamily", &grishaginConstrainedProblems, 
                                                    type_constraned::CONSTR, "GrishaginConstrained"),
                                     problem_family("GKLSProblemConstrainedFamily", &gklsConstrainedProblems, 
                                                    type_constraned::CONSTR, "GKLSConstrained") };

    vector<vector<int>> K{ {0, 700, 25},
                           {0, 1500, 25},
                           {0, 3000, 25},
                           {0, 4500, 25} };

    double eps = 0.01;

    // Parameters of DIRECT
    data_direct_oper_character f_data;
    int max_iter = 10000;
    double magic_eps = 1.0e-4;
    double volume_reltol = 0.0;
    double sigma_reltol  = 0.0;
    vector<direct_algorithm> algorithms{ DIRECT_ORIGINAL, DIRECT_GABLONSKY };

    direct_method direct(f, &f_data, -1, vector<double>{}, vector<double>{}, -1, max_iter, magic_eps,
                                                               volume_reltol, sigma_reltol, nullptr, DIRECT_ORIGINAL);
    f_data.eps = eps;

    // Parameters of mggsa
    int den = 10, incr = 30;
    int maxEvals = 100000;
    vector<vector<double>> r_array{ {3.0, 2.4, 1.6, 1.0},
                                    {4.2, 3.8, 2.0, 1.0},
                                    {3.5, 2.6, 1.8, 1.0},
                                    {4.7, 3.0, 2.0, 1.0} };
    vector<int> key_array{1, 3, 3, 3};
    vector<double> d_array{0.0, 0.0, 0.01, 0.01};

    mggsa_method mggsa(nullptr, -1, -1, vector<double>{}, vector<double>{}, -1.0, -1.0, den, -1, eps, -1, maxEvals, incr);

    int r_size = r_array[0].size();
    vector<vector<vector<double>>> success_rate(number_family);
    for (int i = 0; i < number_family; i++) {
        success_rate[i].resize(r_size + 2);
        for (int j = 0; j < r_size + 2; j++) {
            success_rate[i][j].resize((K[i][1] - K[i][0]) / K[i][2] + 1);
        }
    }

    int total_start_time = clock();
#if defined(CALC_DIRECT)
    ofstream ofstr_direct("output_data/direct_mggsa_operational_characteristics_direct.txt");
    if (!ofstr_direct.is_open()) cerr << "File opening error\n";

    for (int i = 0; i < number_family; i++) {
        for (int j = 0; j < algorithms.size(); j++) {
            vector<double> A, B;
            functor_family func;
            functor_family_constr func_constr;

            if (problems[i].type == type_constraned::CONSTR) {
                func_constr.constr_opt_problem_family = static_cast<IConstrainedOptProblemFamily*>(problems[i].optProblemFamily);
                (*func_constr.constr_opt_problem_family)[0]->GetBounds(A, B);
                direct.setN((*func_constr.constr_opt_problem_family)[0]->GetDimension());
                f_data.functor = &func_constr;
                f_data.type = type_constraned::CONSTR;
            } else {
                func.opt_problem_family = static_cast<IOptProblemFamily*>(problems[i].optProblemFamily);
                (*func.opt_problem_family)[0]->GetBounds(A, B);
                direct.setN((*func.opt_problem_family)[0]->GetDimension());
                f_data.functor = &func;
                f_data.type = type_constraned::NONCONSTR;
            }
            direct.setAB(A, B);
            direct.setMaxFeval(K[i][1]);
            direct.setAlghorithm(algorithms[j]);

            int count_func = problems[i].optProblemFamily->GetFamilySize();
            vector<int> count_evals(count_func);
            int count_successful;

            int start_time = clock();
            for (int k = 0; k < count_func; k++) {
                f_data.count_evals = 0;
                f_data.converge = false;

                if (problems[i].type == type_constraned::CONSTR) {
                    func_constr.current_func = k;
                } else {
                    func.current_func = k;
                }
                direct.solve_test();

                if (f_data.converge) {
                    count_evals[k] = f_data.min_count_evals;
                } else {
                    count_evals[k] = K[i][1] + 1;
                }
            }
            for (int k = K[i][0]; k <= K[i][1]; k += K[i][2]) {
                count_successful = (int)count_if(count_evals.begin(), count_evals.end(), [k](double elem){ return elem <= k; });
                success_rate[i][j][k / K[i][2]] = (double)count_successful / count_func;
            }
            int end_time = clock();
            double work_time = ((double)end_time - start_time) / CLOCKS_PER_SEC;

            string type_direct = (algorithms[j] == DIRECT_ORIGINAL) ? "ORIGINAL" : "GABLONSKY";
            string str_input = "DIRECT " + type_direct + " " + problems[i].name + " time: " + to_string(work_time) + "\n";
            cout << str_input;
        }
    }
    for (int i = 0; i < number_family; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = K[i][0]; k <= K[i][1]; k += K[i][2]) {
                ofstr_direct << k << " " << success_rate[i][j][k / K[i][2]] << endl;
            }
            ofstr_direct << endl << endl;
        }
    }
    ofstr_direct.close();
#endif

#if defined(CALC_MGGSA)
    ofstream ofstr_mggsa("output_data/direct_mggsa_operational_characteristics_mggsa.txt");
    if (!ofstr_mggsa.is_open()) cerr << "File opening error\n";

    const int chunk_mggsa = 2;

#pragma omp parallel for schedule(static, chunk_mggsa) PROC_BIND num_threads(omp_get_num_procs()) collapse(2) \
        shared(number_family, problems, r_array, r_size, success_rate, K) \
        firstprivate(mggsa, key_array, d_array)
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
            
            mggsa.setMaxIters(K[i][1]);
            mggsa.setAB(A, B);
            mggsa.setD(d_array[i]);
            mggsa.setKey(key_array[j]);
            mggsa.setR(r_array[i][j]);

            vector<int> count_iters(problems[i].optProblemFamily->GetFamilySize());
            int count_func = problems[i].optProblemFamily->GetFamilySize();
            vector<double> X_opt;
            int count_successful;
            int countIters, countEvals;

            double start_time = omp_get_wtime();
            for (int k = 0; k < count_func; k++) {
                if (problems[i].type == type_constraned::CONSTR) {
                    func_constr.current_func = k;
                    X_opt = (*func_constr.constr_opt_problem_family)[k]->GetOptimumPoint();
                    mggsa.setF(func_constr);
                } else {
                    func.current_func = k;
                    X_opt = (*func.opt_problem_family)[k]->GetOptimumPoint();
                    mggsa.setF(func);
                }
                if (mggsa.solve_test(X_opt, countIters, countEvals)) {
                    count_iters[k] = countIters;
                } else {
                    count_iters[k] = countIters + 1;
                }
            }
            for (int k = K[i][0]; k <= K[i][1]; k += K[i][2]) {
                count_successful = (int)count_if(count_iters.begin(), count_iters.end(), [k](double elem){ return elem <= k; });
                success_rate[i][j + 2][k / K[i][2]] = (double)count_successful / count_func;
            }
            double end_time = clock();
            double work_time = ((double)end_time - start_time) / CLOCKS_PER_SEC;

            string str_input = "MGGSA: " + problems[i].name + " r = " + to_string(r_array[i][j]) + " key = " + to_string(key_array[j]) + 
                               " time: " + to_string(work_time) + " t_num: " + to_string(omp_get_thread_num()) + "\n";
            cout << str_input;
        }
    }

    for (int i = 0; i < number_family; i++) {
        for (int j = 2; j < r_size + 2; j++) {
            for (int k = K[i][0]; k <= K[i][1]; k += K[i][2]) {
                ofstr_mggsa << k << " " << success_rate[i][j][k / K[i][2]] << endl;
            }
            ofstr_mggsa << endl << endl;
        }
    }
    ofstr_mggsa.close();
#endif
    int total_end_time = clock();
    double total_work_time = ((double)total_end_time - total_start_time) / CLOCKS_PER_SEC;
    cout << "Total time: " << total_work_time << endl;

    int size = (int)key_array.size();
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

    // Plotting operational characteristics(works with gnuplot)
    int error;
#if defined(__linux__)
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/direct_mggsa_operational_characteristics.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
#endif

    char str[100];
    sprintf(str, "gnuplot -c scripts/direct_mggsa_operational_characteristics.gp %d", family_number);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

#if defined(_MSC_VER)
    cin.get();
#endif

	return 0;
}
