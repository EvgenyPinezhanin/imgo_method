#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS
    #define PROC_BIND
#else
    #define PROC_BIND proc_bind(spread)
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <algorithm>
#include <functional>

#include <omp.h>
#include <Grishagin/GrishaginProblemFamily.hpp>
#include <Grishagin/GrishaginConstrainedProblemFamily.hpp>
#include <GKLS/GKLSProblemFamily.hpp>
#include <GKLS/GKLSConstrainedProblemFamily.hpp>
#include <mggsa.h>
#include <task.h>
#include <output_results.h>

using namespace std;

// #define CALC

const int type = 1; // 0 - P_max, 1 - count_trials, 2 - count points, 3 - c_points / c_trials
const int familyNumber = 3; // 0 - Grishagin, 1 - GKLS,
                            // 2 - Grishagin(constrained), 3 - GKLS(constrained)

int main() {
    ofstream ofstr("output_data/mggsa_operational_characteristics_r_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstrOpt("output_data/mggsa_operational_characteristics_r_test_opt.txt");
    if (!ofstrOpt.is_open()) cerr << "File opening error\n";

    const int chunk = 10;

    const int numberFamily = 4;

    TGrishaginProblemFamily grishaginProblems;
    TGrishaginConstrainedProblemFamily grishaginConstrainedProblems;
    TGKLSProblemFamily gklsProblems;
    TGKLSConstrainedProblemFamily gklsConstrainedProblems;

    vector<ProblemFamily> problems{ ProblemFamily("GrishaginProblemFamily", &grishaginProblems, TypeConstrants::NoConstraints,
                                                    "Grishagin"),
                                    ProblemFamily("GKLSProblemFamily", &gklsProblems, TypeConstrants::NoConstraints, "GKLS"),
                                    ProblemFamily("GrishaginProblemConstrainedFamily", &grishaginConstrainedProblems,
                                                    TypeConstrants::Constraints, "GrishaginConstrained"),
                                    ProblemFamily("GKLSProblemConstrainedFamily", &gklsConstrainedProblems,
                                                    TypeConstrants::Constraints, "GKLSConstrained") };

    vector<int> K{ 700, 1200, 2500, 4500 };

    vector<double> rMin{ 1.0, 1.0, 1.0, 1.0 };
    vector<double> rMax{ 5.0, 5.0, 5.0, 5.0 };
    vector<double> step{ 0.05, 0.05, 0.01, 0.01 };
    vector<vector<double>> r(numberFamily);
    vector<int> sizeR(numberFamily);
    for (int i = 0; i < numberFamily; i++) {
        for (double j = rMin[i]; j <= rMax[i]; j += step[i]) {
            r[i].push_back(j);
        }
        sizeR[i] = (int)r[i].size();
    }
    vector<int> key{ 1, 3 };
    size_t sizeKey = key.size();
    vector<double> d{ 0.0, 0.0, 0.01, 0.01 };
    vector<int> maxTrials{ 10000, 15000, 25000, 30000 };

    int numberFunctions;
    vector<vector<vector<vector<int>>>> numberTrialsData(numberFamily), numberTrialPointsData(numberFamily);
    for (int i = 0; i < numberFamily; i++) {
        numberFunctions = problems[i].optProblemFamily->GetFamilySize();
        numberTrialsData[i].resize(sizeKey);
        numberTrialPointsData[i].resize(sizeKey);
        for (int j = 0; j < sizeKey; j++) {
            numberTrialsData[i][j].resize(sizeR[i]);
            numberTrialPointsData[i][j].resize(sizeR[i]);
            for (int k = 0; k < sizeR[i]; k++) {
                numberTrialsData[i][j][k].resize(numberFunctions);
                numberTrialPointsData[i][j][k].resize(numberFunctions);
            }
        }
    }

    vector<vector<vector<double>>> P(numberFamily);
    vector<vector<vector<int>>> maxNumberTrials(numberFamily), maxNumberTrialPoints(numberFamily);
    for (int i = 0; i < numberFamily; i++) {
        P[i].resize(sizeKey);
        maxNumberTrials[i].resize(sizeKey);
        maxNumberTrialPoints[i].resize(sizeKey);
        for (int j = 0; j < sizeKey; j++) {
            P[i][j].resize(r[i].size());
            maxNumberTrials[i][j].resize(r[i].size());
            maxNumberTrialPoints[i][j].resize(r[i].size());
        }
    }

    double totalStartTime = omp_get_wtime();
#if defined(CALC)
    ofstream ofstrData("output_data/mggsa_operational_characteristics_test_r_data.txt");
    if (!ofstrData.is_open()) cerr << "File opening error\n";

    int den = 10, incr = 20;
    double eps = 0.01;
    int maxFevals = 1000000;

    MggsaMethod mggsa(nullptr, -1, -1, vector<double>{}, vector<double>{}, -1.0, -1.0, den, -1, eps, -1, maxFevals, incr);

    int numberR = 0;
    for (int i = 0; i < sizeR.size(); i++) {
        numberR += sizeR[i];
    }

#pragma omp parallel for schedule(static, chunk) PROC_BIND num_threads(omp_get_num_procs()) collapse(2) \
        shared(numberR, d, r, P, numberTrialsData, numberTrialPointsData) \
        firstprivate(K, sizeR, key, problems, mggsa)
    for (int t = 0; t < numberR; t++) {
        for (int j = 0; j < sizeKey; j++) {
            int i, k, tTmp = t;
            for (int l = 0; l < sizeR.size(); l++) {
                if (tTmp < sizeR[l]) {
                    i = l; k = tTmp;
                    break;
                } else {
                    tTmp -= sizeR[l];
                }
            }

            vector<double> A, B;
            FunctorFamily functor;
            FunctorFamilyConstrained functorConstrained;

            if (problems[i].type == TypeConstrants::Constraints) {
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
            mggsa.setR(r[i][k]);
            mggsa.setD(d[i]);
            mggsa.setKey(key[j]);
            mggsa.setMaxTrials(maxTrials[i]);

            vector<double> XOpt;
            int numberTrials, numberFevals;
            int numberFunctions = problems[i].optProblemFamily->GetFamilySize();
            vector<int> numberTrialsArray(numberFunctions), numberTrialPointsArray(numberFunctions);

            double startTime = omp_get_wtime();
            for (int l = 0; l < numberFunctions; l++) {
                if (problems[i].type == TypeConstrants::Constraints) {
                    functorConstrained.currentFunction = l;
                    XOpt = (*functorConstrained.constrainedOptProblemFamily)[l]->GetOptimumPoint();
                    mggsa.setF(functorConstrained);
                } else {
                    functor.currentFunction = l;
                    XOpt = (*functor.optProblemFamily)[l]->GetOptimumPoint();
                    mggsa.setF(functor);
                }
                if (mggsa.solveTest(XOpt, numberTrials, numberFevals)) {
                    numberTrialsArray[l] = numberTrials;
                } else {
                    numberTrialsArray[l] = maxTrials[i] + 1;
                }
                numberTrialPointsArray[l] = mggsa.getNumberTrialPoints();
            }

            numberTrialsData[i][j][k] = numberTrialsArray;
            numberTrialPointsData[i][j][k] = numberTrialPointsArray;
            double endTime = omp_get_wtime();
            double workTime = endTime - startTime;

            string strOutput = problems[i].name + " key = " + to_string(key[j]) + " r = " + to_string(r[i][k]) + 
                               " time: " + to_string(workTime) + " thread number: " + to_string(omp_get_thread_num()) + "\n";
            cout << strOutput;
        }
    }

    for (int i = 0; i < numberFamily; i++) {
        for (int j = 0; j < sizeKey; j++) {
            for (int k = 0; k < r[i].size(); k++) {
                for (int l = 0; l < numberTrialsData[i][j][k].size(); l++) {
                    ofstrData << numberTrialsData[i][j][k][l] << " ";
                }
                ofstrData << endl;
                for (int l = 0; l < numberTrialPointsData[i][j][k].size(); l++) {
                    ofstrData << numberTrialPointsData[i][j][k][l] << " ";
                }
                ofstrData << endl;
            }
        }
    }
    ofstrData.close();
#else
    ifstream ifstrData("output_data/mggsa_operational_characteristics_test_r_data.txt");
    if (!ifstrData.is_open()) cerr << "File opening error\n";

    size_t sizeNumberTrialsData;
    for (int i = 0; i < numberFamily; i++) {
        for (int j = 0; j < sizeKey; j++) {
            for (int k = 0; k < r[i].size(); k++) {
                sizeNumberTrialsData = numberTrialsData[i][j][k].size();
                for (int l = 0; l < sizeNumberTrialsData; l++) {
                    ifstrData >> numberTrialsData[i][j][k][l];
                }
                for (int l = 0; l < sizeNumberTrialsData; l++) {
                    ifstrData >> numberTrialPointsData[i][j][k][l];
                }
            }
        }
    }
    ifstrData.close();
#endif
    double totalEndTime = omp_get_wtime();
    double totalWorkTime = totalEndTime - totalStartTime;
    cout << "Total time: " << totalWorkTime << endl;

    int numberSuccessful, maxStep, index;

    for (int i = 0; i < numberFamily; i++) {
        numberFunctions = problems[i].optProblemFamily->GetFamilySize();
        maxStep = K[i];

        for (int j = 0; j < sizeKey; j++) {
            for (int k = 0; k < r[i].size(); k++) {
                numberSuccessful = (int)count_if(numberTrialsData[i][j][k].begin(), numberTrialsData[i][j][k].end(),
                                                 [maxStep] (double elem) { return elem <= maxStep; });
                P[i][j][k] = (double)numberSuccessful / numberFunctions;

                auto iter = max_element(numberTrialsData[i][j][k].begin(), numberTrialsData[i][j][k].end());
                maxNumberTrialPoints[i][j][k] = *iter;
                for (int l = 0; l < numberTrialsData[i][j][k].size(); l++) {
                    if (*iter == numberTrialsData[i][j][k][l]) index = l;
                }
                maxNumberTrialPoints[i][j][k] = numberTrialPointsData[i][j][k][index];
            }
        }
    }

    for (int i = 0; i < numberFamily; i++) {
        for (int j = 0; j < key.size(); j++) {
            for (int k = 0; k < r[i].size(); k++) {
                ofstr << r[i][k] << " " << P[i][j][k] << " " <<
                         maxNumberTrials[i][j][k] << " " << maxNumberTrialPoints[i][j][k] << endl;
            }
            ofstr << endl << endl;
        }
        ofstr << endl << endl;
    }
    ofstr.close();

    initArrayGnuplot(ofstrOpt, "familyNames", numberFamily);
    initArrayGnuplot(ofstrOpt, "K", numberFamily);
    for (int i = 0; i < numberFamily; i++) {
        setValueInArrayGnuplot(ofstrOpt, "familyNames", i + 1, "\"" + problems[i].shortName + "\"");
        setValueInArrayGnuplot(ofstrOpt, "K", i + 1, to_string(K[i]));
    }
    ofstrOpt.close();

    vector<int> args{ type, familyNumber };
    drawGraphGnuplot("scripts/mggsa_operational_characteristics_test_r.gp", args);

#if defined(_MSC_VER)
    cin.get();
#endif

	return 0;
}
