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
#include <direct_method.h>
#include <task.h>

using namespace std;

// #define CALC

const int family_number = 3; // 0 - Grishagin, 1 - GKLS,
                             // 2 - Grishagin(constrained), 3 - GKLS(constrained),
                             // 4 - comparison Grishagin and GKLS, 5 - comparison Grishagin and GKLS (constrained)

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
#if defined(CALC)
    ofstream ofstr("output_data/direct_operational_characteristics.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("output_data/direct_operational_characteristics_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    vector<vector<int>> K{ {0, 700, 25},
                           {0, 1600, 25},
                           {0, 500, 25},
                           {0, 4000, 25} };

    data_direct_oper_character f_data;
    vector<double> A, B;
    int max_iter = 10000;
    double magic_eps = 1.0e-4, magic_eps_abs = 1.0e-4, eps = 0.01;
    double volume_reltol = 0.0;
    double sigma_reltol  = 0.0;
    vector<direct_algorithm> algorithms{ DIRECT_ORIGINAL, DIRECT_GABLONSKY };

    const int number_family = 4;

    TGrishaginProblemFamily grishaginProblems;
    TGrishaginConstrainedProblemFamily grishaginConstrainedProblems;
    TGKLSProblemFamily gklsProblems;
    TGKLSConstrainedProblemFamily gklsConstrainedProblems;

    vector<vector<int>> count_evals_vec(number_family);
    count_evals_vec[0].resize(grishaginProblems.GetFamilySize(), 0);
    count_evals_vec[1].resize(gklsProblems.GetFamilySize(), 0);
    count_evals_vec[2].resize(grishaginConstrainedProblems.GetFamilySize(), 0);
    count_evals_vec[3].resize(gklsConstrainedProblems.GetFamilySize(), 0);

    vector<problem_family> problems{ problem_family("GrishaginProblemFamily", &grishaginProblems, type_constraned::NONCONSTR, "Grishagin"),
                                     problem_family("GKLSProblemFamily", &gklsProblems, type_constraned::NONCONSTR, "GKLS"),
                                     problem_family("GrishaginProblemConstrainedFamily", &grishaginConstrainedProblems, 
                                                    type_constraned::CONSTR, "GrishaginConstrained"),
                                     problem_family("GKLSProblemConstrainedFamily", &gklsConstrainedProblems, 
                                                    type_constraned::CONSTR, "GKLSConstrained") };

    direct_method direct(f, &f_data, -1, A, B, -1, max_iter, magic_eps, magic_eps_abs, volume_reltol, sigma_reltol, nullptr, DIRECT_ORIGINAL);

    functor_family func;
    functor_family_constr func_constr;
    int count_func, count_successful;
    int start_time, end_time;
    double work_time;

    f_data.eps = eps;

    int total_start_time = clock();
    for (int i = 0; i < number_family; i++) {
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

        count_func = problems[i].optProblemFamily->GetFamilySize();

        cout << problems[i].name << endl;
        for (int j = 0; j < algorithms.size(); j++) {
            cout << "Type of algorithm: " << ((algorithms[j] == DIRECT_ORIGINAL) ? "DIRECT_ORIGINAL" : "DIRECT_GABLONSKY") << endl;
            direct.setAlghorithm(algorithms[j]);

            start_time = clock();
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
                    count_evals_vec[i][k] = f_data.min_count_evals;
                } else {
                    count_evals_vec[i][k] = K[i][1] + 1;
                }
            }
            for (int k = K[i][0]; k <= K[i][1]; k += K[i][2]) {
                count_successful = (int)count_if(count_evals_vec[i].begin(), count_evals_vec[i].end(), [k] (double elem) { return elem <= k; });
                cout << "K = " << k << " success rate = " << (double)count_successful / count_func << endl;
                ofstr << k << " " << (double)count_successful / count_func << endl;
            }
            ofstr << endl << endl;
            end_time = clock();
            work_time = ((double)end_time - start_time) / CLOCKS_PER_SEC;
            cout << "Time: " << work_time << endl;
        }
    }
    ofstr.close();
    int total_end_time = clock();
    int total_work_time = ((double)total_end_time - total_start_time) / CLOCKS_PER_SEC;
    cout << "Total time: " << total_work_time << endl;

    ofstr_opt << "array Name[" << number_family << "]" << endl;
    for (int i = 0; i < number_family; i++) {
        ofstr_opt << "Name[" << i + 1 << "] = \"" << problems[i].short_name << "\"" << endl;
    }
    ofstr_opt.close();
#endif

    // Plotting operational characteristics(works with gnuplot)
    int error;
#if defined(__linux__)
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/direct_operational_characteristics.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
#endif

    char str[100];
    sprintf(str, "gnuplot -c scripts/direct_operational_characteristics.gp %d", family_number);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

#if defined(_MSC_VER)
    cin.get();
#endif

	return 0;
}
