#include<iostream>
#include<fstream>
#include<Hill/HillProblem.hpp>
#include<Hill/HillProblemFamily.hpp>
#include<Shekel/ShekelProblem.hpp>
#include<Shekel/ShekelProblemFamily.hpp>
#include<imgo.h>

using namespace std;

THillProblemFamily hillFam;
TShekelProblemFamily shekelFam;
int current_func;
const int count_func = 2000;

double f_Hill(double x, int j) {
    double res = hillFam[current_func]->ComputeFunction({ x });
    switch (j) {
        case 1: return res;
    }
    return -1.0;
}

double f_Shekel(double x, int j) {
    double res = shekelFam[current_func]->ComputeFunction({ x });
    switch (j) {
        case 1: return res;
    }
    return -1.0;
}

int main() {
    ofstream ofstr("operational_characteristics.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    vector<double> a, b;

    auto f_null = [](double x, int j)->double {return 0.0; };
    imgo_method imgo(f_null, 0, 0, 0, 0, 3.0);

    int K0 = 0;
    int Kmax = 500;
    int Kstep = 10;
    int count_successful = 0;
    double eps = 0.0001;

    vector<bool> bool_Hill(hillFam.GetFamilySize());
    for (int i = 0; i < hillFam.GetFamilySize(); i++) {
        bool_Hill[i] = false;
    }
    vector<bool> bool_Shekel(shekelFam.GetFamilySize());
    for (int i = 0; i < shekelFam.GetFamilySize(); i++) {
        bool_Shekel[i] = false;
    }

    imgo.setEps(eps);
    imgo.setM(0);
    ofstr << count_func << " " << eps << endl;
    for (int i = K0; i <= Kmax; i += Kstep) {
        current_func = 0;
        imgo.setFunc(f_Hill);
        hillFam[0]->GetBounds(a, b);
        imgo.setA(a[0]);
        imgo.setB(b[0]);
        for (int j = 0; j < hillFam.GetFamilySize(); j++) {
            if (!bool_Hill[j]) {
                if (imgo.solve_test(hillFam[j]->GetOptimumPoint()[0], i)) {
                    count_successful++;
                    bool_Hill[j] = true;
                }
            }
            current_func++;
        }

        current_func = 0;
        imgo.setFunc(f_Shekel);
        shekelFam[0]->GetBounds(a, b);
        imgo.setA(a[0]);
        imgo.setB(b[0]);
        for (int j = 0; j < shekelFam.GetFamilySize(); j++) {
            if (!bool_Shekel[j]) {
                if (imgo.solve_test(shekelFam[j]->GetOptimumPoint()[0], i)) {
                    count_successful++;
                    bool_Shekel[j] = true;
                }
            }
            current_func++;
        }
        cout << "K = " << i << " success rate = " << (double)count_successful / count_func << endl;
        ofstr << i << " " << (double)count_successful / count_func << endl;
    }

    ofstr.close();
    cin.get();

	return 0;
}