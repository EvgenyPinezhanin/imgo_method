#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS
#endif

#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>

#include <Hill/HillProblemFamily.hpp>
#include <Shekel/ShekelProblemFamily.hpp>
#include <imgo.h>
#include <task.h>

using namespace std;

#define CALC

const int chart_number = 0; // 0 - Hill, 1 - Shekel, 2 - comparsion Hill and Shekel
const int count_family = 2;
int current_func;

THillProblemFamily hillProblems;
double f_hill(double x, int j) {
    switch (j) {
        case 1: return hillProblems[current_func]->ComputeFunction({x});
        default: return std::numeric_limits<double>::quiet_NaN();
    }
}

TShekelProblemFamily shekelProblems;
double f_shekel(double x, int j) {
    switch (j) {
        case 1: return shekelProblems[current_func]->ComputeFunction({x});
        default: return std::numeric_limits<double>::quiet_NaN();
    }
}

int main() {
#if defined(CALC)
    ofstream ofstr("output_data/imgo_operational_characteristics.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("output_data/imgo_operational_characteristics_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    int count_func;
    int count_successful;
    int count_trials;

    int K0 = 0, Kmax = 500, Kstep = 10;

    vector<double> A, B;
    vector<vector<double>> r_array{ {2.4, 3.2, 3.7},
                                    {2.5, 3.4, 4.0} };
    double eps = 0.0001, d = 0.0;
    int m = 0;

    vector<vector<int>> count_trials_vec(count_family);
    count_trials_vec[0].resize(hillProblems.GetFamilySize(), 0);
    count_trials_vec[1].resize(shekelProblems.GetFamilySize(), 0);

    vector<class_problems_fs> problems{ class_problems_fs("HillProblemFamily", &hillProblems, type_constraned::NONCONSTR, f_hill, "hill"),
                                        class_problems_fs("ShekelProblemFamily", &shekelProblems, 
                                                          type_constraned::NONCONSTR, f_shekel, "shekel") };

    imgo_method imgo(nullptr, m, 0.0, 0.0, -1.0, d, eps);

    IOptProblemFamily *opt_problem_family;
    for (int i = 0; i < count_family; i++) {
        opt_problem_family = static_cast<IOptProblemFamily*>(problems[i].problem);

        for (int j = 0; j < r_array[0].size(); j++) {
            ofstr_opt << "r" << j + 1 << "_" << problems[i].short_name << " = \"" << r_array[i][j] << "\"" << endl; 
        }

        imgo.setF(problems[i].f);
        count_func = problems[i].problem->GetFamilySize();
        cout << problems[i].name << endl;
        ofstr << "# " << problems[i].name << endl;
        for (int j = 0; j < r_array[i].size(); j++) {
            cout << "r = " << r_array[i][j] << endl;
            imgo.setR(r_array[i][j]);
            for (int k = 0; k < count_func; k++) {
                current_func = k;
                count_trials = Kmax;
                (*opt_problem_family)[k]->GetBounds(A, B);
                imgo.setAB(A[0], B[0]);
                if (imgo.solve_test((*opt_problem_family)[k]->GetOptimumPoint()[0], count_trials, Stop::ACCURNUMBER)) {
                    count_trials_vec[i][k] = count_trials;
                } else {
                    count_trials_vec[i][k] = count_trials + 1;
                }
            }
            for (int k = K0; k <= Kmax; k += Kstep) {
                count_successful = (int)count_if(count_trials_vec[i].begin(), count_trials_vec[i].end(), [k](double elem){ return elem <= k; });
                cout << "K = " << k << " success rate = " << (double)count_successful / count_func << endl;
                ofstr << k << " " << (double)count_successful / count_func << endl;
            }
            ofstr << endl << endl;
        }
    }

    ofstr.close();
    ofstr_opt.close();
#endif

    // Plotting operational characteristics(works with gnuplot)
    int error;
#if defined(__linux__)
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/imgo_operational_characteristics.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
#endif

    char str[100];
    sprintf(str, "gnuplot -c scripts/imgo_operational_characteristics.gp %d", chart_number);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

#if defined(_MSC_VER)
    cin.get();
#endif
	return 0;
}
