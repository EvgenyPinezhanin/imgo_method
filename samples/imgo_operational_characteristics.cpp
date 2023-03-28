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

// #define CALC

const int family_number = 2; // 0 - Hill, 1 - Shekel, 2 - comparsion Hill and Shekel

int main() {
#if defined(CALC)
    ofstream ofstr("output_data/imgo_operational_characteristics.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("output_data/imgo_operational_characteristics_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    int count_func, count_successful;
    int countIters, countTrials, countEvals;

    int K0 = 0, Kmax = 500, Kstep = 10;

    vector<double> A, B;
    vector<vector<double>> r_array{ {3.0, 3.2, 3.4},
                                    {3.1, 3.4, 3.7} };
    double eps = 0.0001, d = 0.0;
    int m = 0;

    const int number_family = 2;
    THillProblemFamily hillProblems;
    TShekelProblemFamily shekelProblems;

    vector<vector<int>> count_evals(number_family);
    count_evals[0].resize(hillProblems.GetFamilySize(), 0);
    count_evals[1].resize(shekelProblems.GetFamilySize(), 0);

    vector<problem_family> problems{ problem_family("HillProblemFamily", &hillProblems, type_constraned::NONCONSTR, "Hill"),
                                     problem_family("ShekelProblemFamily", &shekelProblems, type_constraned::NONCONSTR, "Shekel") };

    imgo_method imgo(nullptr, m, 0.0, 0.0, -1.0, d, eps);

    functor_family functor;
    for (int i = 0; i < number_family; i++) {
        functor.opt_problem_family = static_cast<IOptProblemFamily*>(problems[i].optProblemFamily);

        count_func = problems[i].optProblemFamily->GetFamilySize();
        cout << problems[i].name << endl;
        ofstr << "# " << problems[i].name << endl;
        for (int j = 0; j < r_array[i].size(); j++) {
            cout << "r = " << r_array[i][j] << endl;
            imgo.setR(r_array[i][j]);
            for (int k = 0; k < count_func; k++) {
                functor.current_func = k;
                imgo.setF(function<double(double, int)>(functor));
                imgo.setMaxIters(Kmax);
                imgo.setMaxEvals(Kmax);
                (*functor.opt_problem_family)[k]->GetBounds(A, B);
                imgo.setAB(A[0], B[0]);
                if (imgo.solve_test((*functor.opt_problem_family)[k]->GetOptimumPoint()[0], countIters, countTrials, countEvals)) {
                    count_evals[i][k] = countEvals;
                } else {
                    count_evals[i][k] = countEvals + 1;
                }
            }
            for (int k = K0; k <= Kmax; k += Kstep) {
                count_successful = (int)count_if(count_evals[i].begin(), count_evals[i].end(), [k](double elem){ return elem <= k; });
                cout << "K = " << k << " success rate = " << (double)count_successful / count_func << endl;
                ofstr << k << " " << (double)count_successful / count_func << endl;
            }
            ofstr << endl << endl;
        }
    }
    ofstr.close();

    ofstr_opt << "array Name[" << number_family << "]" << endl;
    ofstr_opt << "array R[" << r_array.size() * 3 << "]" << endl;
    for (int i = 0; i < number_family; i++) {
        ofstr_opt << "Name[" << i + 1 << "] = \"" << problems[i].short_name << "\"" << endl;
        for (int j = 0; j < r_array[i].size(); j++) {
            ofstr_opt << "R[" << (i * 3 + 1) + j << "] = \"" << r_array[i][j] << "\"" << endl; 
        }
    }
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
    sprintf(str, "gnuplot -c scripts/imgo_operational_characteristics.gp %d", family_number);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

#if defined(_MSC_VER)
    cin.get();
#endif

	return 0;
}
