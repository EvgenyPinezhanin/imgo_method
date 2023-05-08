#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS
#endif

#include <fstream>
#include <iostream>
#include <vector>
#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include <math.h>
#else
    #include <cmath>
#endif

#include <mggsa.h>
#include <map.h>
#include <output_results.h>
#include <omp.h>

using namespace std;

double euclideanDistance(const vector<double> &X, const vector<double> &Y) {
    double res = 0.0;
    size_t dimension = X.size();
    for (int i = 0; i < dimension; i++) {
        res += (X[i] - Y[i]) * (X[i] - Y[i]);
    }
    return sqrt(res);
}

#define CALC
// #define OUTPUT_INFO

const int type = 2; // 0 - number trials, 1 - number trial points,
                    // 2 - accuracy, 3 - number trial points / number trials
const int nType = 3; // nMin ... nMax

int main() {
    ofstream ofstrOpt("output_data/mggsa_test_incr_opt.txt");
    if (!ofstrOpt.is_open()) cerr << "File opening error\n";

    const int nMin = 2, nMax = 3;
    const int numberN = nMax - nMin + 1;
    const int incrArray[numberN][2] = { {1, 60},
                                        {1, 100} };
    const int mMin = 8, mMax = 10;

#if defined(CALC)
    ofstream ofstr("output_data/mggsa_test_incr.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";

    const int chunk = 10;

    double eps = 0.01, r = 2.1, d = 0.0;
    int numberConstraints = 0, key = 3;
    int maxTrials = 100000, maxFevals = 100000;

    vector<double> XOpt, A, B;
    for (int i = 0; i < nMin - 1; i++) {
        XOpt.push_back(0.0);
        // A.push_back(-5.12);
        // B.push_back(5.12);
        A.push_back(-0.5);
        B.push_back(1.0);
    }

    MggsaMethod mggsa(nullptr, -1, numberConstraints, A, B, r, d, -1, key, eps, maxTrials, maxFevals, -1);

    vector<vector<vector<double>>> accuracyArray(numberN);
    vector<vector<vector<int>>> numberTrialsArray(numberN), numberTrialPointsArray(numberN);
    for (int i = nMin; i <= nMax; i++) {
        accuracyArray[i - nMin].resize(mMax - mMin + 1);
        numberTrialsArray[i - nMin].resize(mMax - mMin + 1);
        numberTrialPointsArray[i - nMin].resize(mMax - mMin + 1);
        for (int j = mMin; j <= mMax; j++) {
            accuracyArray[i - nMin][j - mMin].resize(incrArray[i - nMin][1] - incrArray[i - nMin][0] + 1);
            numberTrialsArray[i - nMin][j - mMin].resize(incrArray[i - nMin][1] - incrArray[i - nMin][0] + 1);
            numberTrialPointsArray[i - nMin][j - mMin].resize(incrArray[i - nMin][1] - incrArray[i - nMin][0] + 1);
        }
    }

    vector<double> lambdas, X;
    int numberTrials, numberFevals;
    double accuracy;

    double totalStartTime = omp_get_wtime();
    for (int i = nMin; i <= nMax; i++) {
        mggsa.setN(i);
        XOpt.push_back(0.0);
        A.push_back(-0.5);
        B.push_back(1.0);
        mggsa.setAB(A, B);

    #pragma omp parallel for schedule(dynamic, chunk) proc_bind(spread) num_threads(omp_get_num_procs()) collapse(2) \
            shared(accuracyArray, numberTrialsArray, numberTrialPointsArray) \
            firstprivate(mggsa, nMin, mMin, incrArray) \
            private(lambdas, X, XOpt, accuracy, numberTrials, numberFevals)
        for (int j = mMin; j <= mMax; j++) {
            for (int k = incrArray[i - nMin][0]; k <= incrArray[i - nMin][1]; k++) {
                double startTime = omp_get_wtime();

                mggsa.setDen(j);
                mggsa.setIncr(k);
                mggsa.setF([n = i] (vector<double> x, int j) -> double {
                                double sum = 0.0;
                                for (int i = 0; i < n; i++) {
                                    sum += x[i] * x[i] - 10.0 * cos(2.0 * M_PI * x[i]);
                                }
                                switch (j) {
                                    case 1: return 10.0 * n + sum;
                                    default: return numeric_limits<double>::quiet_NaN();
                                }
                            });
 
                mggsa.solve(numberTrials, numberFevals, X);
                mggsa.getLambda(lambdas);
 
                accuracy = euclideanDistance(XOpt, X);
                accuracyArray[i - nMin][j - mMin][k - incrArray[i - nMin][0]] = accuracy;
                numberTrialsArray[i - nMin][j - mMin][k - incrArray[i - nMin][0]] = numberTrials;
                numberTrialPointsArray[i - nMin][j - mMin][k - incrArray[i - nMin][0]] = mggsa.getNumberTrialPoints();

                double endTime = omp_get_wtime();
                double workTime = endTime - startTime;

            #if defined(OUTPUT_INFO)
                printResultMggsa("Rastrigin function", N, numberConstraints, A, B, vector<double>(), XOpt,
                                 fRastrigin(XOpt, numberConstraints + 1), maxTrials, maxFevals, eps, r, d, j, key, k,
                                 numberTrials, numberFevals, lambdas, X, fRastrigin(X, numberConstraints + 1));
            #else
                string strOutput = "Rastrigin: n = " + to_string(i) + " m = " + to_string(j) + " incr = " + to_string(k) +
                                   " number of trials = " + to_string(numberTrials) + " number of trial points = " +
                                   to_string(mggsa.getNumberTrialPoints()) + " time: " + to_string(workTime) + " thread number = " +
                                   to_string(omp_get_thread_num()) + "\n";
                cout << strOutput;
            #endif
            }
        }
    }
    double totalEndTime = omp_get_wtime();
    double totalWorkTime = totalEndTime - totalStartTime;
    cout << "Total time: " << totalWorkTime << endl;

    for (int i = nMin; i <= nMax; i++) {
        for (int j = mMin; j <= mMax; j++) {
            for (int k = incrArray[i - nMin][0]; k <= incrArray[i - nMin][1]; k++) {
                ofstr << k << " " << numberTrialsArray[i - nMin][j - mMin][k - incrArray[i - nMin][0]]
                           << " " << numberTrialPointsArray[i - nMin][j - mMin][k - incrArray[i - nMin][0]] 
                           << " " << accuracyArray[i - nMin][j - mMin][k - incrArray[i - nMin][0]] << endl;
            }
            ofstr << endl << endl;
        }
        ofstr << endl << endl;
    }
    ofstr.close();
#endif

    string arraysName[] = { "incrMin", "incrMax" };
    for (int i = 0; i < numberN; i++) {
        initArrayGnuplot(ofstrOpt, arraysName[i], numberN);
    }
    for (int i = 0; i < numberN; i++) {
        for (int j = 0; j < 2; j++) {
            setValueInArrayGnuplot(ofstrOpt, arraysName[j], i + 1, to_string(incrArray[i][j]));
        }
    }
    ofstrOpt.close();

    vector<int> args{ type, nType, nMin, nMax, mMin, mMax };
    drawGraphGnuplot("scripts/mggsa_test_incr.gp", args);

    return 0;
}
