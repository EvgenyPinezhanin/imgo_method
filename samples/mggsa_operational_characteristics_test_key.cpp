#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS
    #define PROC_BIND
#else
    #define PROC_BIND proc_bind(master)
#endif

#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>

#include <Grishagin/GrishaginProblemFamily.hpp>
#include <Grishagin/GrishaginConstrainedProblemFamily.hpp>
#include <GKLS/GKLSProblemFamily.hpp>
#include <GKLS/GKLSConstrainedProblemFamily.hpp>
#include <mggsa.h>
#include <task.h>
#include <output_results.h>
#include <omp.h>

using namespace std;

// #define CALC

const int familyNumber = 3; // 0 - Grishagin, 1 - GKLS,
                            // 2 - Grishagin(constrained), 3 - GKLS(constrained)

int main() {
    ofstream ofstrOpt("output_data/mggsa_operational_characteristics_test_key_opt.txt");
    if (!ofstrOpt.is_open()) cerr << "File opening error\n";

    const int chunk = 2;

    const int numberFamily = 4;

    TGrishaginProblemFamily grishaginProblems;
    TGrishaginConstrainedProblemFamily grishaginConstrainedProblems;
    TGKLSProblemFamily gklsProblems;
    TGKLSConstrainedProblemFamily gklsConstrainedProblems;

    vector<ProblemFamily> problems{ ProblemFamily("GrishaginProblemFamily", &grishaginProblems, TypeConstraints::NoConstraints, "Grishagin"),
                                    ProblemFamily("GKLSProblemFamily", &gklsProblems, TypeConstraints::NoConstraints, "GKLS"),
                                    ProblemFamily("GrishaginProblemConstrainedFamily", &grishaginConstrainedProblems, 
                                                  TypeConstraints::Constraints, "GrishaginConstrained"),
                                    ProblemFamily("GKLSProblemConstrainedFamily", &gklsConstrainedProblems, 
                                                  TypeConstraints::Constraints, "GKLSConstrained") };

    vector<vector<int>> K{ { 0,  700, 25 },
                           { 0, 1200, 25 },
                           { 0, 2500, 25 },
                           { 0, 4500, 25 } };

    vector<vector<double>> r{ { 3.0, 2.8, 2.6, 2.4 },
                              { 4.3, 4.1, 3.9, 3.7 },
                              { 3.0, 2.6, 2.2, 1.8 },
                              { 4.5, 4.1, 3.7, 3.3 } };
    vector<int> key{ 1, 3, 3, 3 };
    vector<double> d{ 0.0, 0.0, 0.01, 0.01 };

#if defined(CALC)
    ofstream ofstr("output_data/mggsa_operational_characteristics_test_key.txt");
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
            FunctorFamily functor;
            FunctorFamilyConstrained functorConstrained;

            if (problems[i].type == TypeConstraints::Constraints) {
                functorConstrained.constrainedOptProblemFamily = static_cast<IConstrainedOptProblemFamily*>(problems[i].optProblemFamily);
                (*functorConstrained.constrainedOptProblemFamily)[0]->GetBounds(A, B);
                mggsa.setN((*functorConstrained.constrainedOptProblemFamily)[0]->GetDimension());
                mggsa.setNumberConstraints((*functorConstrained.constrainedOptProblemFamily)[0]->GetConstraintsNumber());
            } else {
                functor.optProblemFamily = static_cast<IOptProblemFamily*>(problems[i].optProblemFamily);
                (*functor.optProblemFamily)[0]->GetBounds(A, B);
                mggsa.setN((*functor.optProblemFamily)[0]->GetDimension());
                mggsa.setNumberConstraints(0);
            }
            
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
                if (problems[i].type == TypeConstraints::Constraints) {
                    functorConstrained.currentFunction = k;
                    XOpt = (*functorConstrained.constrainedOptProblemFamily)[k]->GetOptimumPoint();
                    mggsa.setF(functorConstrained);
                } else {
                    functor.currentFunction = k;
                    XOpt = (*functor.optProblemFamily)[k]->GetOptimumPoint();
                    mggsa.setF(functor);
                }
                if (mggsa.solveTest(XOpt, numberTrials, numberFevals)) {
                    numberTrialsArray[k] = numberFevals;
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
    setVariableGnuplot(ofstrOpt, "numberKey", to_string(sizeKey));
    initArrayGnuplot(ofstrOpt, "familyNames", numberFamily);
    initArrayGnuplot(ofstrOpt, "r", sizeKey * numberFamily);
    initArrayGnuplot(ofstrOpt, "key", sizeKey);
    for (int i = 0; i < numberFamily; i++) {
        setValueInArrayGnuplot(ofstrOpt, "familyNames", i + 1, problems[i].shortName);
        setValueInArrayGnuplot(ofstrOpt, "key", i + 1, key[i], false);
        for (int j = 0; j < r[i].size(); j++) {
            setValueInArrayGnuplot(ofstrOpt, "r", (i * sizeKey) + j + 1, r[i][j]);
        }
    }
    ofstrOpt.close();

    drawGraphGnuplot("scripts/mggsa_operational_characteristics_test_key.gp", familyNumber);

#if defined(_MSC_VER)
    cin.get();
#endif

	return 0;
}
