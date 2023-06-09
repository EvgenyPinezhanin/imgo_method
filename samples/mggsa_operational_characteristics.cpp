#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS    
#endif

#include <iostream>
#include <fstream>
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

const int familyNumber = 1; // 1 - Grishagin, 2 - GKLS,
                            // 3 - Grishagin(constrained), 4 - GKLS(constrained),
                            // 5 - comparison Grishagin and GKLS, 6 - comparison Grishagin and GKLS(constrained)
const int displayType = 0; // 0 - application, 1 - png

int main() {
    ofstream ofstrOpt("output_data/mggsa_operational_characteristics_opt.txt");
    if (!ofstrOpt.is_open()) cerr << "File opening error\n";

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

    vector<vector<double>> r{ { 2.7, 3.0, 3.3 },
                              { 4.0, 4.3, 4.6 },
                              { 2.7, 3.0, 3.3 },
                              { 4.2, 4.5, 4.8 } };

    vector<vector<int>> K{ { 0, 700, 25 },
                           { 0, 1500, 25 },
                           { 0, 3000, 25 },
                           { 0, 4500, 25 } };

#if defined(CALC)
    ofstream ofstr("output_data/mggsa_operational_characteristics.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";

    int den = 10, key = 1;
    double eps = 0.01, d = 0.0;
    int maxFevals = 100000;

    MggsaMethod mggsa(nullptr, -1, -1, vector<double>{}, vector<double>{}, -1.0, d, den, key, eps, -1, maxFevals);

    vector<double> A, B, XOpt;
    vector<int> numberTrialsArray;
    int numberFunctions, numberSuccessful;
    int numberTrials, numberFevals;
    double startTime, endTime, workTime;
    FunctorFamily functor;
    FunctorFamilyConstrained functorConstrained;

    double totalStartTime = omp_get_wtime();
    for (int i = 0; i < numberFamily; i++) {
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
        numberTrialsArray.resize(problems[i].optProblemFamily->GetFamilySize());
        numberFunctions = problems[i].optProblemFamily->GetFamilySize();
        mggsa.setAB(A, B);

        cout << problems[i].name << endl;
        for (int j = 0; j < r[i].size(); j++) {
            cout << "r = " << r[i][j] << endl;
            mggsa.setR(r[i][j]);

            startTime = omp_get_wtime();
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
                    numberTrialsArray[k] = numberTrials;
                } else {
                    numberTrialsArray[k] = K[i][1] + 1;
                }
            }
            for (int k = K[i][0]; k <= K[i][1]; k += K[i][2]) {
                numberSuccessful = (int)count_if(numberTrialsArray.begin(), numberTrialsArray.end(), [k] (double elem) { return elem <= k; });
                cout << "K = " << k << " success rate = " << (double)numberSuccessful / numberFunctions << endl;
                ofstr << k << " " << (double)numberSuccessful / numberFunctions << endl;
            }
            ofstr << endl << endl;
            endTime = omp_get_wtime();
            workTime = endTime - startTime;
            cout << "Time: " << workTime << endl;
        }
    }
    ofstr.close();
    double totalEndTime = omp_get_wtime();
    double totalWorkTime = totalEndTime - totalStartTime;
    cout << "Total time: " << totalWorkTime << endl;
#endif

    initArrayGnuplot(ofstrOpt, "familyNames", numberFamily);
    initArrayGnuplot(ofstrOpt, "r", r.size() * 3);
    for (int i = 0; i < numberFamily; i++) {
        setValueInArrayGnuplot(ofstrOpt, "familyNames", i + 1, problems[i].shortName);
        for (int j = 0; j < r[i].size(); j++) {
            setValueInArrayGnuplot(ofstrOpt, "r", (i * 3) + j + 1, r[i][j]);
        }
    }
    ofstrOpt.close();

    vector<int> args{ displayType, familyNumber };
    drawGraphGnuplot("scripts/mggsa_operational_characteristics.gp", args);

#if defined(_MSC_VER)
    cin.get();
#endif

	return 0;
}
