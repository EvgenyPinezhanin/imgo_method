#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS    
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <Grishagin/grishagin_function.hpp>
#include <Grishagin/GrishaginConstrainedProblem.hpp>
#include <GKLS/GKLSProblem.hpp>
#include <GKLS/GKLSConstrainedProblem.hpp>
#include <mggsa.h>
#include <task.h>
#include <output_results.h>
#include <omp.h>

using namespace std;

const int familyNumber = 0; // 0 - grishagin, 1 - gkls
                            // 2 - constrained grishagin, 3 - constrained GKLS

int main() {
    ofstream ofstr("output_data/mggsa_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstrOpt("output_data/mggsa_test_opt.txt");
    if (!ofstrOpt.is_open()) cerr << "File opening error\n";

    const int chunk = 4;

    const int numberFunctions = 4;

    TGrishaginProblem grishaginProblem;
    GrishaginConstrainedProblem grishaginConstrainedProblem;
    TGKLSProblem gklsProblem;
    TGKLSConstrainedProblem gklsConstrainedProblem;

    vector<ProblemSingle> problems{ ProblemSingle("Grishagin", &grishaginProblem, TypeConstrants::NoConstraints),
                                    ProblemSingle("GKLS", &gklsProblem, TypeConstrants::NoConstraints),
                                    ProblemSingle("GrishaginConstrained", &grishaginConstrainedProblem, TypeConstrants::Constraints),
                                    ProblemSingle("GKLSConstrained", &gklsConstrainedProblem, TypeConstrants::Constraints) };

    const int incrMin = 1, incrMax = 40;
    const int mMin = 8, mMax = 12;
    const int keyMin = 1, keyMax = 3;

    double eps = 0.001, r = 2.2, d = 0.05;
    int n = 2, den = 10, key = 1;
    int maxTrials = 100000, maxFevals = 100000;

    MggsaMethod mggsa;

    vector<double> X, lambdas;
    vector<vector<double>> points;
    vector<TrialConstrained> trials;
    int numberTrials, numberFevals;

/*     for (int i = 0; i < taskArray.size(); i++) {
        if (taskArray[i].used) {
            mggsa.setF(taskArray[i].f);
            mggsa.setN(taskArray[i].n);
            mggsa.setNumberConstraints(taskArray[i].numberConstraints);
            mggsa.setAB(taskArray[i].A, taskArray[i].B);
            mggsa.setEps(taskArray[i].eps);
            mggsa.setMaxTrials(taskArray[i].maxTrials);
            mggsa.setMaxFevals(taskArray[i].maxFevals);
            mggsa.setR(taskArray[i].r);
            mggsa.setD(taskArray[i].d);
            mggsa.setDen(taskArray[i].den);
            mggsa.setKey(taskArray[i].key);

            mggsa.solve(numberTrials, numberFevals, X);
            mggsa.getLambda(lambdas);

            printResultMggsa(taskArray[i].name, taskArray[i].n, taskArray[i].numberConstraints, taskArray[i].A, taskArray[i].B,
                             taskArray[i].L, taskArray[i].XOpt, taskArray[i].f(taskArray[i].XOpt, taskArray[i].numberConstraints + 1),
                             taskArray[i].maxTrials, taskArray[i].maxFevals, taskArray[i].eps, taskArray[i].r, taskArray[i].d,
                             taskArray[i].den, taskArray[i].key, taskArray[i].incr, numberTrials, numberFevals, lambdas, X,
                             taskArray[i].f(X, taskArray[i].numberConstraints + 1));

            addPointGnuplot(ofstr, taskArray[i].XOpt, taskArray[i].f(taskArray[i].XOpt, taskArray[i].numberConstraints + 1));
            addPointGnuplot(ofstr, X, taskArray[i].f(X, taskArray[i].numberConstraints + 1));

            mggsa.getPoints(points);
            mggsa.getTrialPoints(trials);
            addPointsGnuplot(ofstr, points, trials);
        }
    }
    ofstr.close();

    size_t sizeTaskArray = taskArray.size();
    string arraysName[] = { "AX", "BX", "AY", "BY" };
    for (int i = 0; i < 4; i++) {
        initArrayGnuplot(ofstrOpt, arraysName[i], sizeTaskArray);
    }
    for (int i = 0; i < sizeTaskArray; i++) {
        for (int j = 0; j < 2; j++) {
            setValueInArrayGnuplot(ofstrOpt, arraysName[2 * j], i + 1, taskArray[i].A[j], false);
            setValueInArrayGnuplot(ofstrOpt, arraysName[2 * j + 1], i + 1, taskArray[i].B[j], false);
        }
    }
    ofstrOpt.close();

    drawGraphGnuplot("scripts/mggsa_test.gp", functionNumber); */

#if defined(_MSC_VER)
    cin.get();
#endif

    return 0;
}
