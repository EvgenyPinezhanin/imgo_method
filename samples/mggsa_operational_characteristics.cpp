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

const int family_number = 5; // 0 - Grishagin, 1 - GKLS,
                             // 2 - Grishagin(constrained), 3 - GKLS(constrained),
                             // 4 - comparison Grishagin and GKLS, 5 - comparison Grishagin and GKLS (constrained)

int main() {
#if defined(CALC)
    ofstream ofstr("output_data/mggsa_operational_characteristics.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("output_data/mggsa_operational_characteristics_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    int count_func, count_successful, count_trials;

    int start_time, end_time;
    double work_time;

    vector<vector<int>> K{ {0, 700, 25},
                           {0, 1500, 25},
                           {0, 3000, 25},
                           {0, 4500, 25} };

    int den = 10, key = 1, m = 0;
    double eps = 0.01, r = 0.0, d = 0.0;
    int Nmax = 5000;

    vector<double> A, B, X_opt;
    vector<vector<double>> r_array{ {2.7, 3.0, 3.3},
                                    {4.0, 4.3, 4.6},
                                    {2.7, 3.0, 3.3},
                                    {4.2, 4.5, 4.8} };

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

    mggsa_method mggsa(nullptr, -1, -1, A, B, -1.0, d, den, key, eps, Nmax);

    functor_family func;
    functor_family_constr func_constr;
    for (int i = 0; i < number_family; i++) {
        if (problems[i].type == type_constraned::CONSTR) {
            func_constr.constr_opt_problem_family = static_cast<IConstrainedOptProblemFamily*>(problems[i].optProblemFamily);
            (*func_constr.constr_opt_problem_family)[0]->GetBounds(A, B);
            mggsa.setN((*func_constr.constr_opt_problem_family)[0]->GetDimension());
            mggsa.setM((*func_constr.constr_opt_problem_family)[0]->GetConstraintsNumber());
        } else {
            func.opt_problem_family = static_cast<IOptProblemFamily*>(problems[i].optProblemFamily);
            (*func.opt_problem_family)[0]->GetBounds(A, B);
            mggsa.setN((*func.opt_problem_family)[0]->GetDimension());
            mggsa.setM(0);
        }

        count_trials_vec.resize(problems[i].optProblemFamily->GetFamilySize());
        count_func = problems[i].optProblemFamily->GetFamilySize();
        mggsa.setAB(A, B);

        cout << problems[i].name << endl;
        ofstr << "# " << problems[i].name << endl;
        for (int j = 0; j < r_array[i].size(); j++) {
            cout << "r = " << r_array[i][j] << endl;
            mggsa.setR(r_array[i][j]);
            start_time = clock();
            for (int k = 0; k < count_func; k++) {
                count_trials = K[i][1];
                if (problems[i].type == type_constraned::CONSTR) {
                    func_constr.current_func = k;
                    X_opt = (*func_constr.constr_opt_problem_family)[k]->GetOptimumPoint();
                    mggsa.setF(func_constr);
                } else {
                    func.current_func = k;
                    X_opt = (*func.opt_problem_family)[k]->GetOptimumPoint();
                    mggsa.setF(func);
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

    ofstr_opt << "array Name[" << number_family << "]" << endl;
    ofstr_opt << "array R[" << r_array.size() * 3 << "]" << endl;
    for (int i = 0; i < number_family; i++) {
        ofstr_opt << "Name[" << i + 1 << "]=\"" << problems[i].short_name << "\"" << endl;
        for (int j = 0; j < r_array[i].size(); j++) {
            ofstr_opt << "R[" << (i * 3) + j + 1 << "]=\"" << r_array[i][j] << "\"" << endl; 
        }
    }
    ofstr_opt.close();
#endif

    // Plotting operational characteristics(works with gnuplot)
    int error;
#if defined(__linux__)
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/mggsa_operational_characteristics.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
#endif

    char str[100];
    sprintf(str, "gnuplot -c scripts/mggsa_operational_characteristics.gp %d", family_number);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

#if defined(_MSC_VER)
    cin.get();
#endif

	return 0;
}
