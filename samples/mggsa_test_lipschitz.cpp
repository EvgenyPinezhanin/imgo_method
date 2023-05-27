#if defined( _MSC_VER )
    #define PROC_BIND
#else
    #define PROC_BIND proc_bind(master)
#endif

#include <iostream>
#include <fstream>
#include <vector>

#include <Grishagin/grishagin_function.hpp>
#include <Grishagin/GrishaginConstrainedProblem.hpp>
#include <GKLS/GKLSProblem.hpp>
#include <GKLS/GKLSConstrainedProblem.hpp>
#include <mggsa.h>
#include <map.h>
#include <task.h>
#include <output_results.h>
#include <omp.h>

using namespace std;

using vector4d = vector<vector<vector<vector<double>>>>;

// #define CALC
// #define OUTPUT_INFO

void calculation(MggsaMethod &mggsa, vector4d &lipschitzConst, ProblemSingle problem, int numberFunctions, int maxTrials,
                 int maxFevals, double eps, double r, double d, int key, int den, int incr, const int incrMin, const int mMin,
                 const int keyMin);

const int familyNumber = 3; // 0 - Grishagin, 1 - GKLS
                            // 2 - Grishagin(constrained), 3 - GKLS(constrained)

int main() {
    ofstream ofstrOpt("output_data/mggsa_test_lipschitz_opt.txt");
    if (!ofstrOpt.is_open()) cerr << "File opening error\n";

    const int chunk = 10;

    const int numberFunctions = 4;

    TGrishaginProblem grishaginProblem;
    GrishaginConstrainedProblem grishaginConstrainedProblem;
    TGKLSProblem gklsProblem;
    TGKLSConstrainedProblem gklsConstrainedProblem;

    vector<ProblemSingle> problems{ ProblemSingle("Grishagin", &grishaginProblem, TypeConstraints::NoConstraints),
                                    ProblemSingle("GKLS", &gklsProblem, TypeConstraints::NoConstraints),
                                    ProblemSingle("GrishaginConstrained", &grishaginConstrainedProblem, TypeConstraints::Constraints),
                                    ProblemSingle("GKLSConstrained", &gklsConstrainedProblem, TypeConstraints::Constraints) };

    const int incrMin = 1, incrMax = 40;
    const int mMin = 8, mMax = 12;
    const int keyMin = 1, keyMax = 3;

    vector<double> r{ 3.0, 4.3, 3.0, 4.5 };

#if defined(CALC)
    ofstream ofstr("output_data/mggsa_test_lipschitz.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";

    double eps = 0.01, d = 0.0;
    int maxTrials = 100000, maxFevals = 100000;

    MggsaMethod mggsa(nullptr, -1, -1, vector<double>(), vector<double>(), -1.0, d, -1, -1, eps, maxTrials, maxFevals, -1);

    vector4d lipschitzConst(numberFunctions);
    for (int i = 0; i < numberFunctions; i++) {
        lipschitzConst[i].resize(keyMax - keyMin + 1);
        for (int j = 0; j < keyMax - keyMin + 1; j++) {
            lipschitzConst[i][j].resize(mMax - mMin + 1);
            for (int k = 0; k < mMax - mMin + 1; k++) {
                lipschitzConst[i][j][k].resize(incrMax - incrMin + 1);
            }
        }
    }

    double totalStartTime = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, chunk) PROC_BIND num_threads(omp_get_num_procs()) collapse(4) \
        shared(numberFunctions, keyMin, keyMax, mMin, mMax, incrMin, incrMax, problems, r, lipschitzConst) \
        firstprivate(mggsa, eps, d, maxTrials, maxFevals)
    for (int i = 0; i < numberFunctions; i++) {
        for (int j = keyMin; j <= keyMax; j++) {
            for (int k = mMin; k <= mMax; k++) {
                for (int l = incrMin; l <= incrMax; l++) {
                    if (j != keyMax) {
                        if (l == incrMin) {
                            calculation(mggsa, lipschitzConst, problems[i], i, maxTrials, maxFevals, eps, r[i],
                                        d, j, k, l, incrMin, mMin, keyMin);
                            for (int m = incrMin + 1; m <= incrMax; m++) {
                                lipschitzConst[i][j - keyMin][k - mMin][m - incrMin] = lipschitzConst[i][j - keyMin][k - mMin][0];
                            }
                        }
                    } else {
                        calculation(mggsa, lipschitzConst, problems[i], i, maxTrials, maxFevals, eps, r[i],
                                    d, j, k, l, incrMin, mMin, keyMin);
                    }
                }
            }
        }
    }
    double totalEndTime = omp_get_wtime();
    double totalWorkTime = totalEndTime - totalStartTime;
    cout << "Total time: " << totalWorkTime << endl;

    for (int i = 0; i < numberFunctions; i++) {
        for (int j = keyMin; j <= keyMax; j++) {
            for (int k = mMin; k <= mMax; k++) {
                for (int l = incrMin; l <= incrMax; l++) {
                    ofstr << l << " " << lipschitzConst[i][j - keyMin][k - mMin][l - incrMin] << endl;
                }
                ofstr << endl << endl;
            }
        }
    }
    ofstr.close();
#endif

    initArrayGnuplot(ofstrOpt, "familyNames", numberFunctions);
    for (int i = 0; i < numberFunctions; i++) {
        setValueInArrayGnuplot(ofstrOpt, "familyNames", i + 1, problems[i].name);
    }
    ofstrOpt.close();

    vector<int> args{ familyNumber, keyMin, keyMax, mMin, mMax, incrMin, incrMax };
    drawGraphGnuplot("scripts/mggsa_test_lipschitz.gp", args);

    return 0;
}

void calculation(MggsaMethod &mggsa, vector4d &lipschitzConst, ProblemSingle problem, int numberFunctions, int maxTrials,
                 int maxFevals, double eps, double r, double d, int key, int den, int incr, const int incrMin, const int mMin,
                 const int keyMin) {
    double accuracy, fXOpt, fX;
    int numberConstraints, n;
    int numberTrials, numberFevals;
    vector<double> A, B, XOpt, X, L;
    FunctorSingle functor;
    FunctorSingleConstrained functorConstrained;

    double startTime = omp_get_wtime();

    if (problem.type == TypeConstraints::Constraints) {
        functorConstrained.constrainedOptProblem = static_cast<IConstrainedOptProblem*>(problem.optProblem);
        functorConstrained.constrainedOptProblem->GetBounds(A, B);
        n = functorConstrained.constrainedOptProblem->GetDimension();
        numberConstraints = functorConstrained.constrainedOptProblem->GetConstraintsNumber();
        XOpt = functorConstrained.constrainedOptProblem->GetOptimumPoint();
        mggsa.setF(functorConstrained);
        fXOpt = functorConstrained(XOpt, numberConstraints + 1);
    } else {
        functor.optProblem = static_cast<IOptProblem*>(problem.optProblem);
        functor.optProblem->GetBounds(A, B);
        n = functor.optProblem->GetDimension();
        numberConstraints = 0;
        XOpt = functor.optProblem->GetOptimumPoint();
        mggsa.setF(functor);
        fXOpt = functor(XOpt, 1);
    }

    mggsa.setNumberConstraints(numberConstraints);
    mggsa.setN(n);
    mggsa.setAB(A, B);
    mggsa.setKey(key);
    mggsa.setDen(den);
    mggsa.setIncr(incr);
    mggsa.setR(r);

    mggsa.solve(numberTrials, numberFevals, X);
    mggsa.getL(L);

    if (problem.type == TypeConstraints::Constraints) {
        fX = functorConstrained(X, numberConstraints + 1);
    } else {
        fX = functor(X, 1);
    }
    lipschitzConst[numberFunctions][(size_t)key - keyMin][(size_t)den - mMin][(size_t)incr - incrMin] = L[0];
    accuracy = sqrt((XOpt[0] - X[0]) * (XOpt[0] - X[0]) + (XOpt[1] - X[1]) * (XOpt[1] - X[1]));

    double endTime = omp_get_wtime();
    double workTime = endTime - startTime;

#if defined(OUTPUT_INFO)
    printResultMggsa("Rastrigin function", n, numberConstraints, A, B, vector<double>(), XOpt, fXOpt, maxTrials,
                     maxFevals, eps, r, d, den, key, incr, numberTrials, numberFevals, L, X, fX);
#else
    string strOutput = problem.name + " key = " + to_string(key) + " m = " + to_string(den) + " incr = " + to_string(incr) +
                        " time: " + to_string(workTime) + " thread number = " + to_string(omp_get_thread_num()) + "\n";
    cout << strOutput;
#endif
}
