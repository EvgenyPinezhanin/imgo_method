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
#include <output_results.h>

using namespace std;

#define CALC

const int familyNumber = 2; // 0 - Hill, 1 - Shekel, 2 - comparsion Hill and Shekel

int main() {

#if defined(CALC)
    ofstream ofstr("output_data/imgo_operational_characteristics.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstrOpt("output_data/imgo_operational_characteristics_opt.txt");
    if (!ofstrOpt.is_open()) cerr << "File opening error\n";

    const int numberFamily = 2;

    THillProblemFamily hillProblems;
    TShekelProblemFamily shekelProblems;

    vector<ProblemFamily> problems{ ProblemFamily("HillProblemFamily", &hillProblems, TypeConstraints::NoConstraints, "Hill"),
                                    ProblemFamily("ShekelProblemFamily", &shekelProblems, TypeConstraints::NoConstraints, "Shekel") };

    double eps = 0.0001, d = 0.0;
    int numberConstraints = 0;
    int maxTrials = 100000, maxFevals = 100000;

    ImgoMethod imgo(nullptr, numberConstraints, 0.0, 0.0, -1.0, d, eps, maxTrials, maxFevals);

    vector<vector<double>> r{ {3.0, 3.2, 3.4},
                              {3.1, 3.4, 3.7} };

    int K0 = 0, Kmax = 500, Kstep = 10;

    vector<vector<int>> numberTrialsArray(numberFamily);
    numberTrialsArray[0].resize(hillProblems.GetFamilySize(), 0);
    numberTrialsArray[1].resize(shekelProblems.GetFamilySize(), 0);

    vector<double> A, B;
    int numberFunctions, numberSuccessful;
    int numberTrials, numberFevals;
    FunctorFamily functor;

    for (int i = 0; i < numberFamily; i++) {
        functor.optProblemFamily = static_cast<IOptProblemFamily*>(problems[i].optProblemFamily);
        numberFunctions = problems[i].optProblemFamily->GetFamilySize();

        cout << problems[i].name << endl;
        for (int j = 0; j < r[i].size(); j++) {
            cout << "r = " << r[i][j] << endl;
            imgo.setR(r[i][j]);

            for (int k = 0; k < numberFunctions; k++) {
                functor.currentFunction = k;
                imgo.setF(function<double(double, int)>(functor));
                imgo.setMaxTrials(Kmax);
                (*functor.optProblemFamily)[k]->GetBounds(A, B);
                imgo.setAB(A[0], B[0]);

                if (imgo.solveTest((*functor.optProblemFamily)[k]->GetOptimumPoint()[0], numberTrials, numberFevals)) {
                    numberTrialsArray[i][k] = numberTrials;
                } else {
                    numberTrialsArray[i][k] = Kmax + 1;
                }
            }
            for (int k = K0; k <= Kmax; k += Kstep) {
                numberSuccessful = (int)count_if(numberTrialsArray[i].begin(), numberTrialsArray[i].end(),
                                                 [k] (double elem) { return elem <= k; });
                cout << "K = " << k << " success rate = " << (double)numberSuccessful / numberFunctions << endl;
                ofstr << k << " " << (double)numberSuccessful / numberFunctions << endl;
            }
            ofstr << endl << endl;
        }
    }
    ofstr.close();

    initArrayGnuplot(ofstrOpt, "familyNames", numberFamily);
    initArrayGnuplot(ofstrOpt, "r", r.size() * 3);
    for (int i = 0; i < numberFamily; i++) {
        setValueInArrayGnuplot(ofstrOpt, "familyNames", i + 1, problems[i].shortName);
        for (int j = 0; j < r[i].size(); j++) {
            setValueInArrayGnuplot(ofstrOpt, "r", i * 3 + 1 + j, r[i][j]); 
        }
    }
    ofstrOpt.close();
#endif

    drawGraphGnuplot("scripts/imgo_operational_characteristics.gp", familyNumber);

#if defined(_MSC_VER)
    cin.get();
#endif

	return 0;
}
