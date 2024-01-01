#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS
    #define PROC_BIND
#else
    #define PROC_BIND proc_bind(master)
#endif

#include <iostream>
#include <fstream>
#include <iostream>
#include <algorithm>

#include <Grishagin/GrishaginProblemFamily.hpp>
#include <Grishagin/GrishaginConstrainedProblemFamily.hpp>
#include <GKLS/GKLSProblemFamily.hpp>
#include <GKLS/GKLSConstrainedProblemFamily.hpp>
#include <opt_methods/MggsaMethod.h>
#include <opt_problems/FamilyProblem.h>
#include <gnuplot/OutputFile.h>
#include <gnuplot/Script.h>
#include <omp.h>

using namespace std;

#define CALC
// #define DRAW

const int familyNumber = 0; // 0 - Grishagin(with constraints), 2 - GKLS(with constraints)
const int displayType = 0; // 0 - application, 1 - png, 2 - png(notitle)

int main() {
    OutputFile vars_file("output_data/mggsa_operational_characteristics/vars.txt");
    if (!vars_file.isOpen()) std::cerr << "vars_file opening error\n";

    const int chunk = 1;

    const int numberFamily = 2;
    TGrishaginConstrainedProblemFamily grishaginConstrainedProblems;
    TGKLSConstrainedProblemFamily gklsConstrainedProblems;

    std::vector<FamilyProblem> problems{ FamilyProblem("GrishaginProblemConstrainedFamily", &grishaginConstrainedProblems, 
                                                       "GrishaginConstrained"),
                                         FamilyProblem("GKLSProblemConstrainedFamily", &gklsConstrainedProblems, 
                                                       "GKLSConstrained") };

    vector<vector<int>> K{ { 0, 2500, 25 },
                           { 0, 4500, 25 } };

    vector<vector<double>> r{ { 3.0, 2.6, 2.2, 1.8 },
                              { 4.5, 4.1, 3.7, 3.3 } };
    vector<int> key{ 1, 3, 3, 3 };
    vector<double> d{ 0.0, 0.0, 0.01, 0.01 };

#if defined(CALC)
    ofstream ofstr("output_data/mggsa_operational_characteristics/operational_characteristics.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";

    int den = 10, incr = 30;
    double eps = 0.01;
    int maxFevals = 100000;

    MggsaMethod mggsa(nullptr, -1, -1, vector<double>{}, vector<double>{}, -1.0, -1.0, den, -1, eps, -1, maxFevals, incr);

    int sizeR = r[0].size();
    vector<vector<vector<double>>> successRate(numberFamily);
    for (int i = 0; i < numberFamily; i++) {
        successRate[i].resize(sizeR);
        for (int j = 0; j < sizeR; j++) {
            successRate[i][j].resize((K[i][1] - K[i][0]) / K[i][2] + 1);
        }
    }

    double totalStartTime = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, chunk) PROC_BIND num_threads(omp_get_num_procs()) collapse(2) \
        shared(numberFamily, problems, r, key, d, sizeR, successRate) \
        firstprivate(mggsa)
    for (int i = 0; i < numberFamily; i++) {
        for (int j = 0; j < sizeR; j++) {
            vector<double> A, B;
            FunctorFamilyConstrained functor;

            functor.constrainedOptProblemFamily = static_cast<IConstrainedOptProblemFamily*>(problems[i].optProblemFamily);
            (*functor.constrainedOptProblemFamily)[0]->GetBounds(A, B);
            mggsa.setN((*functor.constrainedOptProblemFamily)[0]->GetDimension());
            mggsa.setNumberConstraints((*functor.constrainedOptProblemFamily)[0]->GetConstraintsNumber());
            
            mggsa.setMaxTrials(K[i][1]);
            mggsa.setAB(A, B);
            mggsa.setD(d[i]);
            mggsa.setKey(key[j]);
            mggsa.setR(r[i][j]);

            vector<double> XOpt;
            int numberSuccessful, numberTrials, numberFevals;
            int numberFunctions = problems[i].optProblemFamily->GetFamilySize();
            vector<int> numberTrialsArray(numberFunctions);

            double startTime = omp_get_wtime();
            for (int k = 0; k < numberFunctions; k++) {
                functor.currentFunction = k;
                XOpt = (*functor.constrainedOptProblemFamily)[k]->GetOptimumPoint();
                mggsa.setF(functor);

                if (mggsa.solveTest(XOpt, numberTrials, numberFevals)) {
                    numberTrialsArray[k] = numberTrials;
                } else {
                    numberTrialsArray[k] = K[i][1] + 1;
                }
            }

            for (int k = K[i][0]; k <= K[i][1]; k += K[i][2]) {
                numberSuccessful = (int)count_if(numberTrialsArray.begin(), numberTrialsArray.end(), [k](double elem){ return elem <= k; });
                successRate[i][j][k / K[i][2]] = (double)numberSuccessful / numberFunctions;
            }
            double endTime = omp_get_wtime();
            double workTime = endTime - startTime;

            string strOutput = problems[i].name + " r = " + to_string(r[i][j]) + " key = " + to_string(key[j]) + 
                               " time: " + to_string(workTime) + " thread number: " + to_string(omp_get_thread_num()) + "\n";
            cout << strOutput;
        }
    }

    for (int i = 0; i < numberFamily; i++) {
        for (int j = 0; j < sizeR; j++) {
            for (int k = K[i][0]; k <= K[i][1]; k += K[i][2]) {
                ofstr << k << " " << successRate[i][j][k / K[i][2]] << endl;
            }
            ofstr << endl << endl;
        }
    }
    ofstr.close();
    double totalEndTime = omp_get_wtime();
    double totalWorkTime = totalEndTime - totalStartTime;
    cout << "Total time: " << totalWorkTime << endl;
#endif

    int sizeKey = key.size();
    vars_file.setVariable("numberKey", sizeKey, false);
    vars_file.initArray("familyName", numberFamily);
    vars_file.initArray("r", sizeKey * numberFamily);
    vars_file.initArray("key", sizeKey);
    for (int i = 0; i < numberFamily; i++) {
        vars_file.setValueInArray("familyName", i + 1, problems[i].shortName);
        for (int j = 0; j < r[i].size(); j++) {
            vars_file.setValueInArray("r", (i * sizeKey) + j + 1, r[i][j]);
        }
    }
    for (int i = 0; i < sizeKey; i++) {
        vars_file.setValueInArray("key", i + 1, key[i]);
    }
    vars_file.close();

#if defined( DRAW )
    Script script("scripts/mggsa_operational_characteristics.gp");
    script.addArgs(std::vector<int>{ displayType, familyNumber });
    script.start();
    if (script.isError() == 2) std::cerr << "Error gnuplot\n";
    if (script.isError() == 1) std::cerr << "Error chmod\n";
#endif

#if defined( _MSC_VER )
    cin.get();
#endif

	return 0;
}
