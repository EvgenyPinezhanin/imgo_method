#include <iostream>
#include <string>
#include <vector>

#include <Hill/HillProblemFamily.hpp>
#include <Shekel/ShekelProblemFamily.hpp>
// #include <tasks/FamilyOptProblemsTask.h>
#include <Solver.h>
#include <opt_methods/GsaMethod.h>
#include <opt_methods/PiyavskyMethod.h>
#include <opt_methods/ScanningMethod.h>
#include <gnuplot/OutputFile.h>
#include <gnuplot/Script.h>
#include <omp.h>

#define CALC
#define DRAW

const int displayType = 1; // 0 - application, 1 - png, 2 - png(notitle)
const int familyNumber = 0; // 0 - Hill, 1 - Shekel

int main() {
/*     const int numberFamily = 2;

    THillProblemFamily hillProblems;
    TShekelProblemFamily shekelProblems;

    vector<> problems{ ProblemFamily("HillProblemFamily", &hillProblems, TypeConstraints::NoConstraints, "Hill"),
                                    ProblemFamily("ShekelProblemFamily", &shekelProblems, TypeConstraints::NoConstraints, "Shekel") };

const std::vector<std::string> methodNames{ "scanning", "piyavsky", "gsa" };

#if defined( CALC )
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

    double totalStartTime = omp_get_wtime();
    for (int i = 0; i < numberFamily; i++) {
        functor.optProblemFamily = static_cast<IOptProblemFamily*>(problems[i].optProblemFamily);
        numberFunctions = problems[i].optProblemFamily->GetFamilySize();

        cout << problems[i].name << endl;
        for (int j = 0; j < r[i].size(); j++) {
            cout << "r = " << r[i][j] << endl;
            imgo.setR(r[i][j]);

            double startTime = omp_get_wtime();
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
            double endTime = omp_get_wtime();
            double workTime = endTime - startTime;
            cout << "Time: " << workTime << endl;
        }
    }
    ofstr.close();
    double totalEndTime = omp_get_wtime();
    double totalWorkTime = totalEndTime - totalStartTime;
    cout << "Total time: " << totalWorkTime << endl;
#endif

#if defined( DRAW )
    VariablesFile variablesFile("output_data/onedimensional_operational_characteristics/vars.txt");
    if (!variablesFile.isOpen()) std::cerr << "Variables file opening error\n";

    
    initArrayGnuplot(ofstrOpt, "familyName", numberFamily);
    initArrayGnuplot(ofstrOpt, "r", r.size() * 3);
    for (int i = 0; i < numberFamily; i++) {
        setValueInArrayGnuplot(ofstrOpt, "familyName", i + 1, problems[i].shortName);
        for (int j = 0; j < r[i].size(); j++) {
            setValueInArrayGnuplot(ofstrOpt, "r", i * 3 + 1 + j, r[i][j]);
        }
    }

    variablesFile.closeFile();

    Script script("scripts/onedimensional_operational_characteristics.gp");
    script.addArgs(std::vector<int>{ displayType, familyNumber });
    script.start();
    if (script.isError() == 2) std::cerr << "Error gnuplot\n";
    if (script.isError() == 1) std::cerr << "Error chmod\n";
#endif

#if defined(_MSC_VER)
    cin.get();
#endif */

	return 0;
}
