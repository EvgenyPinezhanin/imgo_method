#include <iostream>
#include <fstream>
#include <Grishagin/GrishaginProblemFamily.hpp>
#include <imgo.h>

using namespace std;

TGrishaginProblemFamily grishaginProblems;
int current_func;
const int count_func = 100;

double f(vector<double> x, int j) {
    return grishaginProblems[current_func]->ComputeFunction({ x });
};

int main() {
    ofstream ofstr("peano_operational_characteristics.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";

    vector<double> A, B;
    int K0 = 0, Kmax = 1500, Kstep = 100;
    int count_successful = 0;
    int n = 2, den = 10, key = 1, m = 0;
    double eps = 0.001, r = 2.9, d = 0.0;

    auto f_null = [](vector<double> x, int j) -> double { return 0.0; };
    imgo_method imgo(f_null, n, m, A, B, eps, r, d, den, key);
    vector<bool> bool_Grishagin(grishaginProblems.GetFamilySize(), false);

    imgo.setFunc(f);
    ofstr << count_func << " " << eps << endl;
    cout << "GrishaginProblemFamily" << endl;
    for (int i = K0; i <= Kmax; i += Kstep) {
        current_func = 0;
        for (int j = 0; j < grishaginProblems.GetFamilySize(); j++) {
            if (!bool_Grishagin[j]) {
                grishaginProblems[current_func]->GetBounds(A, B);
                imgo.setA(A);
                imgo.setB(B);
                imgo.setN(grishaginProblems[j]->GetDimension());
                if (imgo.solve_test(grishaginProblems[j]->GetOptimumPoint(), i)) {
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

//std::cout << "Parameters for constructing the Peano curve:" << std::endl;
//    std::cout << "n = " << n << " m = " << den << " key = " << key << std::endl;
//    std::cout << "Trials result:" << std::endl;
//    std::cout << "Number of trials = " << number_trials << std::endl;
//    std::cout << "x*_min = " << x_opt << " y*_min = " << y_opt << std::endl;
//    std::cout << "x_min = " << X[0] << " y_min = " << X[1] << std::endl;
//    std::cout << "||X* - X|| = " << sqrt((x_opt - X[0]) * (x_opt - X[0]) + (y_opt - X[1]) * (y_opt - X[1])) << std::endl;