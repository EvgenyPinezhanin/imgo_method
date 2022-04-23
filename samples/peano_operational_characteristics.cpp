#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <Grishagin/GrishaginProblemFamily.hpp>
#include <GKLS/GKLSProblemFamily.hpp>
#include <imgo.h>
#include <omp.h>

using namespace std;

int current_func;

TGrishaginProblemFamily grishaginProblems;
double f_grishagin(vector<double> x, int j) {
    switch (j) {
        case 1: return grishaginProblems[current_func]->ComputeFunction({ x });
        default: return numeric_limits<double>::quiet_NaN();
    }
};

TGKLSProblemFamily GKLSProblems;
double f_gkls(vector<double> x, int j) {
    switch (j) {
        case 1: return GKLSProblems[current_func]->ComputeFunction({ x });
        default: return numeric_limits<double>::quiet_NaN();
    }
};

int main() {
    ofstream ofstr("peano_operational_characteristics.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("peano_operational_characteristics_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    int count_func;

    int start_time, end_time;
    double work_time;

    int K0 = 0, Kmax = 2200, Kstep = 25;
    int count_successful = 0;
    int count_trials;

    vector<double> A, B;
    vector<double> grish_r_array{2.8, 2.9, 3.0};
    vector<double> gkls_r_array{4.2, 4.3, 4.7};
    int n = 2, den = 10, key = 1, m = 0;
    double eps = 0.001, r = 2.9, d = 0.0;

    vector<int> count_Grishagin(grishaginProblems.GetFamilySize(), 0);
    vector<int> count_GKLS(GKLSProblems.GetFamilySize(), 0);

    imgo_method imgo(nullptr, n, m, A, B, r, d, eps, 0, den, key);

    for (int i = 0; i < grish_r_array.size(); i++) {
        ofstr_opt << "r" << i + 1 << "_grish = \"" << grish_r_array[i] << "\"" << endl; 
    }

    imgo.setFunc(f_grishagin);
    count_func = grishaginProblems.GetFamilySize();
    imgo.setN(grishaginProblems[0]->GetDimension());
    cout << "GrishaginProblemFamily" << endl;
    ofstr << "# GrishaginProblemFamily" << endl;
    for (int i = 0; i < grish_r_array.size(); i++) {
        cout << "r = " << grish_r_array[i] << endl;
        imgo.setR(grish_r_array[i]);
        start_time = clock();
        for (int i = 0; i < count_func; i++) {
            current_func = i;
            count_trials = Kmax;
            grishaginProblems[current_func]->GetBounds(A, B);
            imgo.setAB(A, B);
            if (imgo.solve_test(grishaginProblems[i]->GetOptimumPoint(), count_trials, ACCURNUMBER)) {
                count_Grishagin[i] = count_trials;
            } else {
                count_Grishagin[i] = count_trials + 1;
            }
            // cout << i << endl;
        }
        for (int i = K0; i <= Kmax; i += Kstep) {
            count_successful = count_if(count_Grishagin.begin(), count_Grishagin.end(), [i](double elem){ return elem <= i; });
            cout << "K = " << i << " success rate = " << (double)count_successful / count_func << endl;
            ofstr << i << " " << (double)count_successful / count_func << endl;
        }
        ofstr << endl << endl;
        end_time = clock();
        work_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
        cout << "time: " << work_time << endl;
    }

    n = 2; den = 20; key = 1; m = 0;
    eps = 0.001; r = 4.3; d = 0.0;
    Kmax = 4200;

    for (int i = 0; i < gkls_r_array.size(); i++) {
        ofstr_opt << "r" << i + 1 << "_gkls = \"" << gkls_r_array[i] << "\"" << endl; 
    }

    imgo.setFunc(f_gkls);
    imgo.setDen(den);
    count_func = GKLSProblems.GetFamilySize();
    imgo.setN(GKLSProblems[0]->GetDimension());
    cout << "GKLSProblemFamily" << endl;
    ofstr << "# GKLSProblemFamily" << endl;

    start_time = clock();
    for (int i = 0; i < gkls_r_array.size(); i++) {
        cout << "r = " << gkls_r_array[i] << endl;
        imgo.setR(gkls_r_array[i]);
        start_time = clock();
        for (int i = 0; i < count_func; i++) {
            current_func = i;
            count_trials = Kmax;
            GKLSProblems[current_func]->GetBounds(A, B);
            imgo.setAB(A, B);
            if (imgo.solve_test(GKLSProblems[i]->GetOptimumPoint(), count_trials, ACCURNUMBER)) {
                count_GKLS[i] = count_trials;
            } else {
                count_GKLS[i] = count_trials + 1;
            }
            // cout << i << endl;
        }
        for (int i = K0; i <= Kmax; i += Kstep) {
            count_successful = count_if(count_GKLS.begin(), count_GKLS.end(), [i](double elem){ return elem <= i; });
            cout << "K = " << i << " success rate = " << (double)count_successful / count_func << endl;
            ofstr << i << " " << (double)count_successful / count_func << endl;
        }
        ofstr << endl << endl;
        end_time = clock();
        work_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
        cout << "time: " << work_time << endl;
    }

    ofstr.close();
    ofstr_opt.close();

    // Рисование графиков операционной характеристики
#if defined(__linux__)
    int error;
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x peano_oper_characteristics.gp");
    if (error != 0) {
        cerr << "Error chmod" << std::endl;
    }
    error = system("gnuplot -p -c peano_oper_characteristics.gp");
    if (error != 0) {
        cerr << "Error gnuplot" << std::endl;
    }
#endif

#if defined( _MSC_VER )
    cin.get();
#endif
	return 0;
}

// #pragma omp parallel for schedule(dynamic) num_threads(NUM_THREADS) shared(count_func) private(current_func)
// cout << omp_get_thread_num() << endl;
// const int NUM_THREADS = 8;