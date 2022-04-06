#include<iostream>
#include<fstream>
#include<Grishagin/GrishaginProblemFamily.hpp>
#include<Shekel/ShekelProblemFamily.hpp>
#include<imgo.h>

using namespace std;

TGrishaginProblemFamily grishaginFam;
int current_func;
const int count_func = 100;

double f(vector<double> x, int j) {
    return grishaginFam[current_func]->ComputeFunction({ x });
};

int main() {
    ofstream ofstr("peano_operational_characteristics.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";

    vector<double> A, B;
    int K0 = 0, Kmax = 500, Kstep = 10;
    int count_successful = 0;
    int n = 2, den = 10, key = 1, m = 0;
    double eps = 0.0001, r = 2.0, d = 0.0;

    auto f_null = [](vector<double> x, int j) -> double { return 0.0; };
    imgo_method imgo(f_null, n, m, A, B, eps, r, d, den, key);
    vector<bool> bool_Grishagin(grishaginFam.GetFamilySize(), false);

    imgo.setFunc(f);
    ofstr << count_func << " " << eps << endl;
    for (int i = K0; i <= Kmax; i += Kstep) {
        current_func = 0;
        for (int j = 0; j < grishaginFam.GetFamilySize(); j++) {
            if (!bool_Grishagin[j]) {
                grishaginFam[current_func]->GetBounds(A, B);
                imgo.setA(A);
                imgo.setB(B);
                if (imgo.solve_test(grishaginFam[j]->GetOptimumPoint(), i)) {
                    count_successful++;
                    bool_Grishagin[j] = true;
                }
            }
            current_func++;
        }
        cout << "K = " << i << " success rate = " << (double)count_successful / count_func << endl;
        ofstr << i << " " << (double)count_successful / count_func << endl;
    }

    ofstr.close();
    #if defined( _MSC_VER )
        cin.get();
    #endif
	return 0;
}