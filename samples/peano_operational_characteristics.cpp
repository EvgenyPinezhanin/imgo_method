#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <Grishagin/GrishaginProblemFamily.hpp>
#include <imgo.h>

using namespace std;

TGrishaginProblemFamily grishaginProblems;
int current_func;
const int count_func = 100;

double f(vector<double> x, int j) {
    switch (j) {
        case 1: return grishaginProblems[current_func]->ComputeFunction({ x });
        default: return numeric_limits<double>::quiet_NaN();
    }
};

int main() {
    ofstream ofstr("peano_operational_characteristics.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";

    int start_time, end_time;
    double work_time;

    vector<double> A, B;
    int K0 = 0, Kmax = 1500, Kstep = 50;
    int count_successful = 0;
    int n = 2, den = 10, key = 1, m = 0;
    double eps = 0.001, r = 2.9, d = 0.0;
    int count;

    imgo_method imgo(nullptr, n, m, A, B, r, d, eps, 0, den, key);
    vector<int> count_Grishagin(grishaginProblems.GetFamilySize(), 0);

    imgo.setFunc(f);
    ofstr << count_func << " " << eps << endl;

    start_time = clock();
    cout << "GrishaginProblemFamily" << endl;
    for (int i = 0; i < count_func; i++) {
        current_func = i;
        grishaginProblems[current_func]->GetBounds(A, B);
        imgo.setA(A);
        imgo.setB(B);
        imgo.setN(grishaginProblems[i]->GetDimension());
        imgo.solve_test(grishaginProblems[i]->GetOptimumPoint(), count);
        count_Grishagin[i] = count;
    }
    for (int i = K0; i <= Kmax; i += Kstep) {
        count_successful = count_if(count_Grishagin.begin(), count_Grishagin.end(), [i](double elem){ return elem <= i; });
        cout << "K = " << i << " success rate = " << (double)count_successful / count_func << endl;
        ofstr << i << " " << (double)count_successful / count_func << endl;
    }
    end_time = clock();

    work_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    cout << "time: " << work_time << endl;

    ofstr.close();
    #if defined( _MSC_VER )
        cin.get();
    #endif
	return 0;
}
