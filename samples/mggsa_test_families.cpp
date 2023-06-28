#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS
    #define PROC_BIND
#else
    #define PROC_BIND proc_bind(spread)
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

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

const int familyNumber = 3; // 1 - Grishagin, 2 - GKLS
                            // 3 - constrained Grishagin, 4 - constrained GKLS
const int displayType = 2; // 0 - application, 1 - png, 2 - png(notitle)

void functionGrid(vector<ofstream> &ofstr, Functor *functorFamily, const vector<double> &A, const vector<double> &B,
                  double gridStep, TypeConstraints type, double minValue, int numberConstraints = 0);

int main() {
#if defined( CALC )
    ofstream ofstr("output_data/mggsa_test_families.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstrOpt("output_data/mggsa_test_families_opt.txt");
    if (!ofstrOpt.is_open()) cerr << "File opening error\n";

    const int chunk = 1;

    const int numberFamily = 4;

    int numThreads = min(numberFamily, omp_get_num_threads());

    TGrishaginProblemFamily grishaginProblems;
    TGrishaginConstrainedProblemFamily grishaginConstrainedProblems;
    TGKLSProblemFamily gklsProblems;
    TGKLSConstrainedProblemFamily gklsConstrainedProblems;

    vector<ProblemFamily> problems{ ProblemFamily("GrishaginProblemFamily", &grishaginProblems, TypeConstraints::NoConstraints,
                                                  "Grishagin"),
                                    ProblemFamily("GKLSProblemFamily", &gklsProblems, TypeConstraints::NoConstraints, "GKLS"),
                                    ProblemFamily("GrishaginProblemConstrainedFamily", &grishaginConstrainedProblems,
                                                  TypeConstraints::Constraints, "GrishaginConstrained"),
                                    ProblemFamily("GKLSProblemConstrainedFamily", &gklsConstrainedProblems,
                                                  TypeConstraints::Constraints, "GKLSConstrained") };

    vector<vector<ofstream>> ofstrFamily(numberFamily);
    for (int i = 0; i < numberFamily; i++) {
        for (int j = 0; j < 2; j++) {
            ofstrFamily[i].push_back(ofstream("output_data/mggsa_test_families_" + problems[i].shortName +
                                              "_" + to_string(j + 1) + ".txt"));
        }
    }

    vector<int> functionNumber{ 0, 0, 0, 0 };

    double gridStep = 0.01;

    vector<double> eps{ 0.01, 0.01, 0.01, 0.01 };
    vector<double> r{ 2.2, 2.2, 2.2, 3.0 };
    vector<double> d{ 0.0, 0.0, 0.02, 0.001 };

    vector<int> den{ 10, 10, 10, 10 };
    vector<int> key{ 3, 3, 3, 3 };
    vector<int> incr{ 1, 1, 1, 1 };

    vector<int> maxTrials{ 100000, 100000, 100000, 100000 };
    vector<int> maxFevals{ 100000, 100000, 100000, 100000 };

    MggsaMethod mggsa;

    vector<vector<vector<double>>> trialPoints(numberFamily);
    vector<vector<double>> pointOptArray(numberFamily), pointArray(numberFamily);
    vector<vector<double>> AArray(numberFamily), BArray(numberFamily);
    vector<double> minValue(numberFamily);

#pragma omp parallel for schedule(dynamic, chunk) PROC_BIND num_threads(numThreads) \
        shared(problems, functionNumber, ofstrFamily, eps, r, d, den, key, incr, maxTrials, \
               maxFevals, trialPoints, pointOptArray, pointArray, AArray, BArray, minValue) \
        firstprivate(mggsa)
    for (int i = 0; i < numberFamily; i++) {
        vector<double> A, B, X, pointOpt;
        FunctorFamily functor;
        FunctorFamilyConstrained functorConstrained;
        int numberTrials, numberFevals, numberConstraints;

        if (problems[i].type == TypeConstraints::Constraints) {
            functorConstrained.constrainedOptProblemFamily = static_cast<IConstrainedOptProblemFamily*>(problems[i].optProblemFamily);
            (*functorConstrained.constrainedOptProblemFamily)[0]->GetBounds(A, B);
            pointOpt = (*functorConstrained.constrainedOptProblemFamily)[functionNumber[i]]->GetOptimumPoint();
            minValue[i] = (*functorConstrained.constrainedOptProblemFamily)[functionNumber[i]]->GetOptimumValue();
            pointOpt.push_back(minValue[i]);
            pointOptArray[i] = pointOpt;
            mggsa.setF(functorConstrained);
            mggsa.setN((*functorConstrained.constrainedOptProblemFamily)[0]->GetDimension());
            numberConstraints = (*functorConstrained.constrainedOptProblemFamily)[0]->GetConstraintsNumber();
            mggsa.setNumberConstraints(numberConstraints);
        } else {
            functor.optProblemFamily = static_cast<IOptProblemFamily*>(problems[i].optProblemFamily);
            (*functor.optProblemFamily)[0]->GetBounds(A, B);
            pointOpt = (*functor.optProblemFamily)[functionNumber[i]]->GetOptimumPoint();
            minValue[i] = (*functor.optProblemFamily)[functionNumber[i]]->GetOptimumValue();
            pointOpt.push_back(minValue[i]);
            pointOptArray[i] = pointOpt;
            mggsa.setF(functor);
            mggsa.setN((*functor.optProblemFamily)[0]->GetDimension());
            numberConstraints = 0;
            mggsa.setNumberConstraints(numberConstraints);
        }

        mggsa.setAB(A, B);
        mggsa.setEps(eps[i]);
        mggsa.setR(r[i]);
        mggsa.setD(d[i]);
        mggsa.setDen(den[i]);
        mggsa.setKey(key[i]);
        mggsa.setMaxTrials(maxTrials[i]);
        mggsa.setMaxFevals(maxFevals[i]);

        mggsa.solve(numberTrials, numberFevals, X);

        mggsa.getPoints(trialPoints[i]);

        if (problems[i].type == TypeConstraints::Constraints) {
            X.push_back(functorConstrained(X, numberConstraints + 1));
            functionGrid(ofstrFamily[i], &functorConstrained, A, B, gridStep, problems[i].type, minValue[i], numberConstraints);
        } else {
            X.push_back(functor(X, numberConstraints + 1));
            functionGrid(ofstrFamily[i], &functor, A, B, gridStep, problems[i].type, minValue[i]);
        }
        pointArray[i] = X;
        AArray[i] = A;
        BArray[i] = B;
    }

    for (int i = 0; i < numberFamily; i++) {
        for (int j = 0; j < ofstrFamily[i].size(); j++) {
            ofstrFamily[i][j].close();
        }
    }

    for (int i = 0; i < numberFamily; i++) {
        addPointGnuplot(ofstr, pointOptArray[i]);
        addPointGnuplot(ofstr, pointArray[i]);
        addPointsGnuplot(ofstr, trialPoints[i]);
    }
    ofstr.close();

    initArrayGnuplot(ofstrOpt, "familyNames", numberFamily);
    initArrayGnuplot(ofstrOpt, "functionNumber", numberFamily);
    initArrayGnuplot(ofstrOpt, "constrained", numberFamily);
    initArrayGnuplot(ofstrOpt, "minValue", numberFamily);
    initArrayGnuplot(ofstrOpt, "A", numberFamily * 2);
    initArrayGnuplot(ofstrOpt, "B", numberFamily * 2);
    for (int i = 0; i < numberFamily; i++) {
        setValueInArrayGnuplot(ofstrOpt, "familyNames", i + 1, problems[i].shortName);
        setValueInArrayGnuplot(ofstrOpt, "functionNumber", i + 1, functionNumber[i] + 1, false);
        setValueInArrayGnuplot(ofstrOpt, "constrained", i + 1, problems[i].type == TypeConstraints::Constraints ? 1 : 0, false);
        setValueInArrayGnuplot(ofstrOpt, "minValue", i + 1, minValue[i], false);
        for (int j = 0; j < 2; j++) {
            setValueInArrayGnuplot(ofstrOpt, "A", 2 * i + j + 1, AArray[i][j], false);
            setValueInArrayGnuplot(ofstrOpt, "B", 2 * i + j + 1, BArray[i][j], false);
        }
    }
    ofstrOpt.close();
#endif

    vector<int> args{ displayType, familyNumber };
    drawGraphGnuplot("scripts/mggsa_test_families.gp", args);

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}

void functionGrid(vector<ofstream> &ofstr, Functor *functorFamily, const vector<double> &A, const vector<double> &B,
                  double gridStep, TypeConstraints type, double minValue, int numberConstraints) {
    int N = floor((B[0] - A[0]) / gridStep);
    double eps = gridStep / 2.0;

    ofstr[0] << N + 1 << " ";
    for (double i = A[0]; i <= B[0] + eps; i += gridStep) {
        ofstr[0] << i << " ";
    }
    ofstr[0] << "\n";
    for (double i = A[1]; i <= B[1] + eps; i += gridStep) {
	    ofstr[0] << i << " ";
        for (double j = A[0]; j <= B[0] + eps; j += gridStep) {
	        ofstr[0] << (*functorFamily)(vector<double>{ j, i }, numberConstraints + 1) << " ";
        }
        ofstr[0] << "\n";
    }
    ofstr[0] << endl;

    bool constrained;

    ofstr[1] << N + 1 << " ";
    for (double i = A[0]; i <= B[0] + eps; i += gridStep) {
        ofstr[1] << i << " ";
    }
    ofstr[1] << "\n";
    for (double i = A[1]; i <= B[1] + eps; i += gridStep) {
	    ofstr[1] << i << " ";
        for (double j = A[0]; j <= B[0] + eps; j += gridStep) {
            constrained = type == TypeConstraints::Constraints ? true : false;
            for (int k = 0; k < numberConstraints; k++) {
                if ((*functorFamily)(vector<double>{ j, i }, k + 1) > 0.0) constrained = false;
            }
            if (constrained) {
	            ofstr[1] << (*functorFamily)(vector<double>{ j, i }, numberConstraints + 1) << " ";
            } else {
                ofstr[1] << minValue - 1.0 << " ";
            }
        }
        ofstr[1] << "\n";
    }
    ofstr[1] << endl;
}
