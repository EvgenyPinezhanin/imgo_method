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

#define CALC

const int family_number = 5; // 0 - Grishagin, 1 - GKLS,
                             // 2 - Grishagin(constrained), 3 - GKLS(constrained),
                             // 4 - comparison Grishagin and GKLS, 5 - comparison Grishagin and GKLS (constrained)

double f(int n, const double *X, int *undefined_flag, void *data) {
    data_direct_family *f_data = static_cast<data_direct_family*>(data);
    f_data->count_trials++;
    vector<double> point(n);
    for (int i = 0; i < n; i++) {
        point[i] = X[i];
    }
    f_data->points.push_back(point);

    if (f_data->type == type_constraned::CONSTR) {
        functor_family_constr *problem = static_cast<functor_family_constr*>(f_data->functor);

        bool is_constr = true;
        int constr = (*problem->constr_opt_problem_family)[problem->current_func]->GetConstraintsNumber();
        for (int i = 0; i < constr; i++) {
            if ((*problem)(point, i) > 0.0) is_constr = false;
        }
        if (is_constr) *undefined_flag = 1;

        return (*problem)(point, constr + 1);
    } else {
        functor_family *problem = static_cast<functor_family*>(f_data->functor);

        return (*problem)(point, 1);
    }
}

int main() {
#if defined(CALC)
    ofstream ofstr("output_data/direct_operational_characteristics.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("output_data/direct_operational_characteristics_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    int count_func, count_successful, count_trials;

    int start_time, end_time, total_start_time, total_end_time;
    double work_time, total_work_time;

    vector<vector<int>> K{ {0, 700, 25},
                           {0, 1500, 25},
                           {0, 3000, 25},
                           {0, 4500, 25} };

    data_direct_family f_data;
    vector<double> A, B;
    int max_feval = 10000, max_iter = 10000;
    double magic_eps = 1.0e-4, eps = 0.01;
    double sigma_reltol  = 0.0;
    vector<direct_algorithm> algorithms{ DIRECT_ORIGINAL, DIRECT_GABLONSKY };

    double minf;
    vector<double> X_opt, X;

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

    vector<problem_family> problems{ problem_family("GrishaginProblemFamily", &grishaginProblems, type_constraned::NONCONSTR, "Grishagin"),
                                     problem_family("GKLSProblemFamily", &GKLSProblems, type_constraned::NONCONSTR, "GKLS"),
                                     problem_family("GrishaginProblemConstrainedFamily", &grishaginConstrainedProblems, 
                                                    type_constraned::CONSTR, "GrishaginConstrained"),
                                     problem_family("GKLSProblemConstrainedFamily", &GKLSConstrainedProblems, 
                                                    type_constraned::CONSTR, "GKLSConstrained") };

    direct_method direct(f, &f_data, -1, A, B, max_feval, max_iter, magic_eps, -1.0, sigma_reltol, nullptr, DIRECT_ORIGINAL);

    double volume_reltol;
    functor_family func;
    functor_family_constr func_constr;
    for (int i = 0; i < number_family; i++) {
        if (problems[i].type == type_constraned::CONSTR) {
            func_constr.constr_opt_problem_family = static_cast<IConstrainedOptProblemFamily*>(problems[i].optProblemFamily);
            (*func_constr.constr_opt_problem_family)[0]->GetBounds(A, B);
            direct.setN((*func_constr.constr_opt_problem_family)[0]->GetDimension());
            direct.setVolumeReltol(eps / sqrt((*func_constr.constr_opt_problem_family)[0]->GetDimension()));
            f_data.functor = &func_constr;
        } else {
            func.opt_problem_family = static_cast<IOptProblemFamily*>(problems[i].optProblemFamily);
            (*func.opt_problem_family)[0]->GetBounds(A, B);
            direct.setN((*func.opt_problem_family)[0]->GetDimension());
            direct.setVolumeReltol(eps / sqrt((*func.opt_problem_family)[0]->GetDimension()));
            f_data.functor = &func;
        }

        count_trials_vec.resize(problems[i].optProblemFamily->GetFamilySize());
        count_func = problems[i].optProblemFamily->GetFamilySize();
        direct.setAB(A, B);

        cout << problems[i].name << endl;
        for (int j = 0; j < algorithms.size(); j++) {
            cout << "Type of algorithm: " << ((algorithms[j] == DIRECT_ORIGINAL) ? "DIRECT_ORIGINAL" : "DIRECT_GABLONSKY") << endl;
            direct.setAlghorithm(algorithms[j]);
            start_time = clock();
            for (int k = 0; k < count_func; k++) {
                count_trials = K[i][1];
                if (problems[i].type == type_constraned::CONSTR) {
                    func_constr.current_func = k;
                    X_opt = (*func_constr.constr_opt_problem_family)[k]->GetOptimumPoint();
                } else {
                    func.current_func = k;
                    X_opt = (*func.opt_problem_family)[k]->GetOptimumPoint();
                }
                // if (direct.solve_test(X_opt, count_trials, Stop::ACCURNUMBER)) {
                //     count_trials_vec[i][k] = count_trials;
                // } else {
                //     count_trials_vec[i][k] = count_trials + 1;
                // }
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

    ofstr_opt << "array Name[" << number_family << "]" << endl;
    for (int i = 0; i < number_family; i++) {
        ofstr_opt << "Name[" << i + 1 << "]=\"" << problems[i].short_name << "\"" << endl;
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
