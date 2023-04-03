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
#include <mggsa.h>
#include <task.h>

using namespace std;

#define CALC

const int family_number = 2; // 0 - Grishagin, 1 - GKLS,
                             // 2 - Grishagin(constrained), 3 - GKLS(constrained)

const int number_family = 4;

const int chunk = 2;

int main() {
#if defined(CALC)
    ofstream ofstr("output_data/mggsa_operational_characteristics_key_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("output_data/mggsa_operational_characteristics_key_test_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    // vector<vector<int>> K{ {0, 700, 25},
    //                        {0, 1200, 25},
    //                        {0, 2200, 25},
    //                        {0, 4500, 25} };

    vector<vector<int>> K{ {0, 700, 25},
                           {0, 1500, 25},
                           {0, 5000, 25},
                           {0, 7000, 25} };

    int den = 10, incr = 30;
    double eps = 0.01;
    vector<vector<double>> r_array{ {3.0, 2.4, 1.6, 1.0},
                                    {4.2, 3.8, 2.0, 1.0},
                                    {3.5, 2.6, 1.8, 1.0},
                                    {4.7, 3.0, 2.0, 1.0} };
    vector<int> key_array{1, 3, 3, 3};
    vector<double> d_array{0.0, 0.0, 0.01, 0.01};

    int r_size = r_array[0].size();
    vector<vector<vector<double>>> success_rate(number_family);
    for (int i = 0; i < number_family; i++) {
        success_rate[i].resize(r_size);
        for (int j = 0; j < r_size; j++) {
            success_rate[i][j].resize((K[i][1] - K[i][0]) / K[i][2] + 1);
        }
    }

    TGrishaginProblemFamily grishaginProblems;
    TGrishaginConstrainedProblemFamily grishaginConstrainedProblems;
    TGKLSProblemFamily gklsProblems;
    TGKLSConstrainedProblemFamily gklsConstrainedProblems;

    vector<problem_family> problems{ problem_family("GrishaginProblemFamily", &grishaginProblems, type_constraned::NONCONSTR,
                                                    "Grishagin"),
                                     problem_family("GKLSProblemFamily", &gklsProblems, type_constraned::NONCONSTR, "GKLS"),
                                     problem_family("GrishaginProblemConstrainedFamily", &grishaginConstrainedProblems, 
                                                    type_constraned::CONSTR, "GrishaginConstrained"),
                                     problem_family("GKLSProblemConstrainedFamily", &gklsConstrainedProblems, 
                                                    type_constraned::CONSTR, "GKLSConstrained") };

    mggsa_method mggsa(nullptr, -1, -1, vector<double>{}, vector<double>{}, -1.0, -1.0, den, -1, eps, -1, -1, incr);

    int total_start_time = clock();
#pragma omp parallel for schedule(static, chunk) PROC_BIND num_threads(omp_get_num_procs()) collapse(2) \
        shared(number_family, problems, r_array, r_size, success_rate) \
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
            mggsa.setMaxEvals(K[i][1]);
            mggsa.setAB(A, B);
            mggsa.setD(d_array[i]);
            mggsa.setKey(key_array[j]);
            mggsa.setR(r_array[i][j]);

            vector<int> count_evals(problems[i].optProblemFamily->GetFamilySize());
            int count_func = problems[i].optProblemFamily->GetFamilySize();
            vector<double> X_opt;
            int count_successful;
            int countIters, countTrials, countEvals;

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
                if (mggsa.solve_test(X_opt, countIters, countTrials, countEvals)) {
                    count_evals[k] = countEvals;
                } else {
                    count_evals[k] = countEvals + 1;
                }
            }
            for (int k = K[i][0]; k <= K[i][1]; k += K[i][2]) {
                count_successful = (int)count_if(count_evals.begin(), count_evals.end(), [k](double elem){ return elem <= k; });
                success_rate[i][j][k / K[i][2]] = (double)count_successful / count_func;
            }
            double end_time = clock();
            double work_time = ((double)end_time - start_time) / CLOCKS_PER_SEC;

            string str_input = problems[i].name + " r = " + to_string(r_array[i][j]) + " key = " + to_string(key_array[j]) + 
                               " time: " + to_string(work_time) + " t_num: " + to_string(omp_get_thread_num()) + "\n";
            cout << str_input;
        }
    }
    for (int i = 0; i < number_family; i++) {
        for (int j = 0; j < r_size; j++) {
            for (int k = K[i][0]; k <= K[i][1]; k += K[i][2]) {
                ofstr << k << " " << success_rate[i][j][k / K[i][2]] << endl;
            }
            ofstr << endl << endl;
        }
    }
    ofstr.close();
    int total_end_time = clock();
    double total_work_time = ((double)total_end_time - total_start_time) / CLOCKS_PER_SEC;
    cout << "Total time: " << total_work_time << endl;

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
