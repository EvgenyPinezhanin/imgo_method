#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <Grishagin/GrishaginProblemFamily.hpp>
#include <GKLS/GKLSProblemFamily.hpp>
#include <imgo.h>

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

    int count_func;

    int start_time, end_time;
    double work_time;

    int K0 = 0, Kmax = 1500, Kstep = 50;
    int count_successful = 0;
    int count_trials;

    vector<double> A, B;
    int n = 2, den = 10, key = 1, m = 0;
    double eps = 0.001, r = 2.9, d = 0.0;

    vector<int> count_Grishagin(grishaginProblems.GetFamilySize(), 0);
    vector<int> count_GKLS(GKLSProblems.GetFamilySize(), 0);

    imgo_method imgo(nullptr, n, m, A, B, r, d, eps, 0, den, key);

    imgo.setFunc(f_grishagin);
    count_func = grishaginProblems.GetFamilySize();
    ofstr << "# " << count_func << " " << eps << endl;
    cout << "GrishaginProblemFamily" << endl;

    start_time = clock();
    for (int i = 0; i < count_func; i++) {
        current_func = i;
        grishaginProblems[current_func]->GetBounds(A, B);
        imgo.setA(A);
        imgo.setB(B);
        imgo.setN(grishaginProblems[i]->GetDimension());
        imgo.solve_test(grishaginProblems[i]->GetOptimumPoint(), count_trials);
        count_Grishagin[i] = count_trials;
    }
    for (int i = K0; i <= Kmax; i += Kstep) {
        count_successful = count_if(count_Grishagin.begin(), count_Grishagin.end(), [i](double elem){ return elem <= i; });
        cout << "K = " << i << " success rate = " << (double)count_successful / count_func << endl;
        ofstr << i << " " << (double)count_successful / count_func << endl;
    }
    end_time = clock();
    work_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    cout << "time: " << work_time << endl;

    imgo.setFunc(f_gkls);
    r = 5.0;
    imgo.setR(r);
    count_func = GKLSProblems.GetFamilySize();
    ofstr << "# " << count_func << " " << eps << endl;
    cout << "GKLSProblemFamily" << endl;

    start_time = clock();
    for (int i = 0; i < count_func; i++) {
        current_func = i;
        GKLSProblems[current_func]->GetBounds(A, B);
        imgo.setA(A);
        imgo.setB(B);
        imgo.setN(GKLSProblems[i]->GetDimension());
        imgo.solve_test(GKLSProblems[i]->GetOptimumPoint(), count_trials);
        count_GKLS[i] = count_trials;
    }
    for (int i = K0; i <= Kmax; i += Kstep) {
        count_successful = count_if(count_GKLS.begin(), count_GKLS.end(), [i](double elem){ return elem <= i; });
        cout << "K = " << i << " success rate = " << (double)count_successful / count_func << endl;
        ofstr << i << " " << (double)count_successful / count_func << endl;
    }
    end_time = clock();
    work_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    cout << "time: " << work_time << endl;

    ofstr.close();

    // Рисование графиков операционной характеристики
#if defined(__linux__)
    int error;
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x chart.gp");
    if (error != 0) {
        cerr << "Error chmod" << std::endl;
    }
    error = system("gnuplot -p -c peano_oper_char.gp");
    if (error != 0) {
        cerr << "Error gnuplot" << std::endl;
    }
#endif

#if defined( _MSC_VER )
    cin.get();
#endif
	return 0;
}
