#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS
    #define PROC_BIND
#else
    #define PROC_BIND proc_bind(master)
#endif

#include <iostream>
#include <fstream>
#include <string>
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

const int type = 2; // 0 - number of trials, 1 - number of trial points, 
                    // 2 - number of trial points / number of trials
const int familyNumber = 3; // 1 - Grishagin, 2 - GKLS,
                            // 3 - constrained Grishagin, 4 - constrained GKLS
const int displayType = 2; // 0 - application, 1 - png, 2 - png(notitle)

int main() {
    ofstream ofstrOpt("output_data/mggsa_operational_characteristics_test_incr_opt.txt");
    if (!ofstrOpt.is_open()) cerr << "File opening error\n";

    const int chunk = 4;

    const int numberFamily = 4;

    TGrishaginProblemFamily grishaginProblems;
    TGKLSProblemFamily gklsProblems;
    TGrishaginConstrainedProblemFamily grishaginConstrainedProblems;
    TGKLSConstrainedProblemFamily gklsConstrainedProblems;

    vector<ProblemFamily> problems{ ProblemFamily("GrishaginProblemFamily", &grishaginProblems, TypeConstraints::NoConstraints, "Grishagin"),
                                    ProblemFamily("GKLSProblemFamily", &gklsProblems, TypeConstraints::NoConstraints, "GKLS"),
                                    ProblemFamily("GrishaginProblemConstrainedFamily", &grishaginConstrainedProblems, 
                                                  TypeConstraints::Constraints, "GrishaginConstrained"),
                                    ProblemFamily("GKLSProblemConstrainedFamily", &gklsConstrainedProblems, 
                                                  TypeConstraints::Constraints, "GKLSConstrained") };

    vector<vector<double>> r{ { 3.0, 2.8, 2.6, 2.4 },
                              { 4.3, 4.1, 3.9, 3.7 },
                              { 3.0, 2.6, 2.2, 1.8 },
                              { 4.5, 4.1, 3.7, 3.3 } };
    vector<int> incr{ 0, 20, 40, 60 };
    vector<double> d{ 0.0, 0.0, 0.01, 0.01 };

#if defined(CALC)
    ofstream ofstr("output_data/mggsa_operational_characteristics_test_incr.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";

    int den = 10, key = 3;
    double eps = 0.01;
    int maxTrials = 30000, maxFevals = 100000;

    MggsaMethod mggsa(nullptr, -1, -1, vector<double>{}, vector<double>{}, -1.0, -1.0, den, key, eps, maxTrials, maxFevals, -1);

    vector<vector<vector<int>>> maxNumberTrials(numberFamily), maxNumberTrialPoints(numberFamily);
    for (int i = 0; i < numberFamily; i++) {
        maxNumberTrials[i].resize(r[i].size());
        maxNumberTrialPoints[i].resize(r[i].size());
        for (int j = 0; j < r[i].size(); j++) {
            maxNumberTrials[i][j].resize(incr.size());
            maxNumberTrialPoints[i][j].resize(incr.size());
        }
    }

    int sizeR = r[0].size();

    double totalStartTime = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, chunk) PROC_BIND num_threads(omp_get_num_procs()) collapse(2) \
        shared(numberFamily, problems, d, r, sizeR, incr, maxNumberTrials, maxNumberTrialPoints) \
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

            mggsa.setAB(A, B);
            mggsa.setD(d[i]);
            mggsa.setR(r[i][j]);

            int numberFunctions = problems[i].optProblemFamily->GetFamilySize();
            vector<bool> converge(numberFunctions);
            vector<int> numberTrialsArray(numberFunctions), numberTrialPointsArray(numberFunctions);
            vector<double> XOpt;
            int numberTrials, numberFevals;
            double startTime, endTime, workTime;

            for (int k = 0; k < incr.size(); k++) {
                mggsa.setIncr(incr[k]);

                startTime = omp_get_wtime();
                for (int l = 0; l < numberFunctions; l++) {
                    if (problems[i].type == TypeConstraints::Constraints) {
                        functorConstrained.currentFunction = l;
                        XOpt = (*functorConstrained.constrainedOptProblemFamily)[l]->GetOptimumPoint();
                        mggsa.setF(functorConstrained);
                    } else {
                        functor.currentFunction = l;
                        XOpt = (*functor.optProblemFamily)[l]->GetOptimumPoint();
                        mggsa.setF(functor);
                    }
                    converge[l] = mggsa.solveTest(XOpt, numberTrials, numberFevals);
                    numberTrialsArray[l] = numberTrials;
                    numberTrialPointsArray[l] = mggsa.getNumberTrialPoints();
                }
                endTime = omp_get_wtime();
                workTime = endTime - startTime;

                int maxTmp = -1, index = -1;
                for (int i = 0; i < numberFunctions; i++) {
                    if (numberTrialsArray[i] > maxTmp && converge[i]) {
                        maxTmp = numberTrialsArray[i];
                        index = i;
                    }
                }
                if (index != -1) {
                    maxNumberTrials[i][j][k] = maxTmp;
                    maxNumberTrialPoints[i][j][k] = numberTrialPointsArray[index];
                } else {
                    maxNumberTrials[i][j][k] = 1;
                    maxNumberTrialPoints[i][j][k] = 0;
                }

                string strOutput = problems[i].name + " r = " + to_string(r[i][j]) + " incr = " + to_string(incr[k]) + 
                                   " number of trials = " + to_string(maxNumberTrials[i][j][k]) + " number of trial points = " + 
                                   to_string(maxNumberTrialPoints[i][j][k]) + " time: " + to_string(workTime) + " thread number: " +
                                   to_string(omp_get_thread_num()) + "\n";
                cout << strOutput;
            }
        }
    }
    double totalEndTime = omp_get_wtime();
    double totalWorkTime = totalEndTime - totalStartTime;
    cout << "Total time: " << totalWorkTime << endl;

    for (int i = 0; i < numberFamily; i++) {
        for (int j = 0; j < r[i].size(); j++) {
            for (int k = 0; k < incr.size(); k++) {
                ofstr << incr[k] << " " << maxNumberTrials[i][j][k] << " " << maxNumberTrialPoints[i][j][k] << endl;
            }
            ofstr << endl << endl;
        }
        ofstr << endl << endl;
    }
    ofstr.close();
#endif

    size_t sizeIncr = incr.size();
    setVariableGnuplot(ofstrOpt, "numberKey", sizeIncr, false);
    initArrayGnuplot(ofstrOpt, "familyNames", numberFamily);
    initArrayGnuplot(ofstrOpt, "r", sizeIncr * numberFamily);
    for (int i = 0; i < numberFamily; i++) {
        setValueInArrayGnuplot(ofstrOpt, "familyNames", i + 1, problems[i].shortName);
        for (int j = 0; j < r[i].size(); j++) {
            setValueInArrayGnuplot(ofstrOpt, "r", (i * sizeIncr) + j + 1, r[i][j]);
        }
    }
    ofstrOpt.close();

    vector<int> args{ displayType, type, familyNumber };
    drawGraphGnuplot("scripts/mggsa_operational_characteristics_test_incr.gp", args);

#if defined(_MSC_VER)
    cin.get();
#endif

    return 0;
}
