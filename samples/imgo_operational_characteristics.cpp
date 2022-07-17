#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>

#include <Hill/HillProblemFamily.hpp>
#include <Shekel/ShekelProblemFamily.hpp>
#include <imgo.h>

using namespace std;

#define CALC
const int chart_number = 0; // 0 - Hill, 1 - Shekel, 2 - comparsion Hill and Shekel
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
    vector<double> hill_r_array{2.4, 3.2, 3.7};
    vector<double> shekel_r_array{2.5, 3.4, 4.0};
    double eps = 0.0001, d = 0.0;
    int m = 0;

    vector<int> count_Hill(hillProblems.GetFamilySize(), 0);
    vector<int> count_Shekel(shekelProblems.GetFamilySize(), 0);

    imgo_method imgo(nullptr);
    imgo.setEps(eps);
    imgo.setD(d);

    for (int i = 0; i < hill_r_array.size(); i++) {
        ofstr_opt << "r" << i + 1 << "_hill = \"" << hill_r_array[i] << "\"" << endl; 
    }

    imgo.setF(f_hill);
    count_func = hillProblems.GetFamilySize();
    cout << "HillProblemFamily" << endl;
    ofstr << "# HillProblemFamily" << endl;
    for (int i = 0; i < hill_r_array.size(); i++) {
        cout << "r = " << hill_r_array[i] << endl;
        imgo.setR(hill_r_array[i]);
        for (int i = 0; i < count_func; i++) {
            current_func = i;
            count_trials = Kmax;
            hillProblems[current_func]->GetBounds(A, B);
            imgo.setAB(A[0], B[0]);
            if (imgo.solve_test(hillProblems[i]->GetOptimumPoint()[0], count_trials, ACCURNUMBER)) {
                count_Hill[i] = count_trials;
            } else {
                count_Hill[i] = count_trials + 1;
            }
        }
        for (int i = K0; i <= Kmax; i += Kstep) {
            count_successful = count_if(count_Hill.begin(), count_Hill.end(), [i](double elem){ return elem <= i; });
            cout << "K = " << i << " success rate = " << (double)count_successful / count_func << endl;
            ofstr << i << " " << (double)count_successful / count_func << endl;
        }
        ofstr << endl << endl;
    }

    for (int i = 0; i < shekel_r_array.size(); i++) {
        ofstr_opt << "r" << i + 1 << "_shekel = \"" << shekel_r_array[i] << "\"" << endl; 
    }

    imgo.setF(f_shekel);
    count_func = shekelProblems.GetFamilySize();
    cout << "ShekelProblemFamily" << endl;
    ofstr << "# ShekelProblemFamily" << endl;
    for (int i = 0; i < shekel_r_array.size(); i++) {
        cout << "r = " << shekel_r_array[i] << endl;
        imgo.setR(shekel_r_array[i]);
        for (int i = 0; i < count_func; i++) {
            current_func = i;
            count_trials = Kmax;
            shekelProblems[current_func]->GetBounds(A, B);
            imgo.setAB(A[0], B[0]);
            if (imgo.solve_test(shekelProblems[i]->GetOptimumPoint()[0], count_trials, ACCURNUMBER)) {
                count_Shekel[i] = count_trials;
            } else {
                count_Shekel[i] = count_trials + 1;
            }
        }
        for (int i = K0; i <= Kmax; i += Kstep) {
            count_successful = count_if(count_Shekel.begin(), count_Shekel.end(), [i](double elem){ return elem <= i; });
            cout << "K = " << i << " success rate = " << (double)count_successful / count_func << endl;
            ofstr << i << " " << (double)count_successful / count_func << endl;
        }
        ofstr << endl << endl;
    }

    ofstr.close();
    ofstr_opt.close();
#endif

    // Drawing graphs of operational characteristics
    int error;
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/imgo_operational_characteristics.gp");
    if (error != 0) {
        cerr << "Error chmod" << std::endl;
    }

    char str[100];
    sprintf(str, "gnuplot -p -c scripts/imgo_operational_characteristics.gp %d", chart_number);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << std::endl;
    }

#if defined(_MSC_VER)
    cin.get();
#endif
	return 0;
}
