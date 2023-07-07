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

/* #include <Grishagin/GrishaginProblemFamily.hpp>
#include <Grishagin/GrishaginConstrainedProblemFamily.hpp>
#include <GKLS/GKLSProblemFamily.hpp>
#include <GKLS/GKLSConstrainedProblemFamily.hpp>
#include <direct_method.h>
#include <StronginMethods/mggsa.h>
#include <task.h>
#include <output_results.h>
#include <omp.h>

using namespace std; */

// #define CALC_DIRECT
// #define CALC_MGGSA

/* const int familyNumber = 0; // 0 - Grishagin, 1 - GKLS,
                            // 2 - Grishagin(with constraints), 3 - GKLS(with constraints)
const int displayType = 0; // 0 - application, 1 - png, 2 - png(notitle)

double f(int n, const double *X, int *undefinedFlag, void *data) {
    DataDirectOperationalCharacteristics *fData = static_cast<DataDirectOperationalCharacteristics*>(data);
    fData->numberFevals++;
    vector<double> point(n), optPoint(n);
    for (int i = 0; i < n; i++) {
        point[i] = X[i];
    }

    double f;
    if (fData->type == TypeConstraints::Constraints) {
        FunctorFamilyConstrained *problem = static_cast<FunctorFamilyConstrained*>(fData->functor);

        optPoint = (*problem->constrainedOptProblemFamily)[problem->currentFunction]->GetOptimumPoint();

        bool constrained = true;
        int constraints = (*problem->constrainedOptProblemFamily)[problem->currentFunction]->GetConstraintsNumber();
        for (int i = 0; i < constraints; i++) {
            if ((*problem)(point, i) > 0.0) constrained = false;
        }
        if (constrained) *undefinedFlag = 1;

        f = (*problem)(point, constraints + 1);
    } else {
        FunctorFamily *problem = static_cast<FunctorFamily*>(fData->functor);

        optPoint = (*problem->optProblemFamily)[problem->currentFunction]->GetOptimumPoint();

        f = (*problem)(point, 1);
    }
    point.push_back(f);
    fData->points.push_back(point);

    if (!fData->converge) {
        double distance = 0.0;
        for (int i = 0; i < n; i++) {
            distance += (point[i] - optPoint[i]) * (point[i] - optPoint[i]);
        }
        if (distance <= fData->eps) {
            fData->converge = true;
            fData->minNumberFevals = fData->numberFevals;
        }
    }

    return f;
} */

int main() {
/*     ofstream ofstrOpt("output_data/direct_mggsa_operational_characteristics_opt.txt");
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

    vector<vector<int>> K{ { 0, 700,  25 },
                           { 0, 1500, 25 },
                           { 0, 3000, 25 },
                           { 0, 4500, 25 } };

    double eps = 0.01;
    DataDirectOperationalCharacteristics fData(nullptr, TypeConstraints::NoConstraints, eps);

    // Parameters of DIRECT
    int maxIters = 10000;
    double magicEps = 1.0e-4;
    double volumeReltol = 0.0;
    double sigmaReltol  = 0.0;
    vector<direct_algorithm> algorithms{ DIRECT_ORIGINAL, DIRECT_GABLONSKY };

    DirectMethod direct(f, &fData, -1, vector<double>{}, vector<double>{}, -1, maxIters, magicEps,
                        volumeReltol, sigmaReltol, nullptr, DIRECT_ORIGINAL);

    // Parameters of mggsa
    int den = 10, incr = 30;
    int maxFevals = 100000;
    vector<vector<double>> r{ { 3.0, 2.6 },
                              { 4.3, 3.9 },
                              { 3.0, 2.2 },
                              { 4.5, 4.1 } };
    vector<int> key{ 1, 3 };
    vector<double> d{ 0.0, 0.0, 0.01, 0.01 };

    MggsaMethod mggsa(nullptr, -1, -1, vector<double>{}, vector<double>{}, -1.0, -1.0, den, -1, eps, -1, maxFevals, incr);

    int sizeR = r[0].size();
    vector<vector<vector<double>>> successRate(numberFamily);
    for (int i = 0; i < numberFamily; i++) {
        successRate[i].resize(sizeR + 2);
        for (int j = 0; j < sizeR + 2; j++) {
            successRate[i][j].resize((K[i][1] - K[i][0]) / K[i][2] + 1);
        }
    }

    double totalStartTime = omp_get_wtime();
#if defined(CALC_DIRECT)
    ofstream ofstrDirect("output_data/direct_mggsa_operational_characteristics_direct.txt");
    if (!ofstrDirect.is_open()) cerr << "File opening error\n";

    for (int i = 0; i < numberFamily; i++) {
        for (int j = 0; j < algorithms.size(); j++) {
            vector<double> A, B;
            FunctorFamily functor;
            FunctorFamilyConstrained functorConstrained;

            if (problems[i].type == TypeConstraints::Constraints) {
                functorConstrained.constrainedOptProblemFamily = static_cast<IConstrainedOptProblemFamily*>(problems[i].optProblemFamily);
                (*functorConstrained.constrainedOptProblemFamily)[0]->GetBounds(A, B);
                direct.setN((*functorConstrained.constrainedOptProblemFamily)[0]->GetDimension());
                fData.functor = &functorConstrained;
                fData.type = TypeConstraints::Constraints;
            } else {
                functor.optProblemFamily = static_cast<IOptProblemFamily*>(problems[i].optProblemFamily);
                (*functor.optProblemFamily)[0]->GetBounds(A, B);
                direct.setN((*functor.optProblemFamily)[0]->GetDimension());
                fData.functor = &functor;
                fData.type = TypeConstraints::NoConstraints;
            }

            direct.setAB(A, B);
            direct.setMaxFevals(K[i][1]);
            direct.setAlghorithm(algorithms[j]);

            int numberFunctions = problems[i].optProblemFamily->GetFamilySize();
            vector<int> numberFevals(numberFunctions);
            int numberSuccessful;

            double startTime = omp_get_wtime();
            for (int k = 0; k < numberFunctions; k++) {
                fData.numberFevals = 0;
                fData.converge = false;

                if (problems[i].type == TypeConstraints::Constraints) {
                    functorConstrained.currentFunction = k;
                } else {
                    functor.currentFunction = k;
                }
                direct.solveTest();

                if (fData.converge) {
                    numberFevals[k] = fData.minNumberFevals;
                } else {
                    numberFevals[k] = K[i][1] + 1;
                }
            }
            for (int k = K[i][0]; k <= K[i][1]; k += K[i][2]) {
                numberSuccessful = (int)count_if(numberFevals.begin(), numberFevals.end(), [k] (double elem) { return elem <= k; });
                successRate[i][j][k / K[i][2]] = (double)numberSuccessful / numberFunctions;
            }
            double endTime = omp_get_wtime();
            double workTime = endTime - startTime;

            string typeDirect = (algorithms[j] == DIRECT_ORIGINAL) ? "ORIGINAL" : "GABLONSKY";
            string strOutput = "DIRECT " + typeDirect + " " + problems[i].name + " time: " + to_string(workTime) + "\n";
            cout << strOutput;
        }
    }
    for (int i = 0; i < numberFamily; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = K[i][0]; k <= K[i][1]; k += K[i][2]) {
                ofstrDirect << k << " " << successRate[i][j][k / K[i][2]] << endl;
            }
            ofstrDirect << endl << endl;
        }
    }
    ofstrDirect.close();
#endif

#if defined(CALC_MGGSA)
    ofstream ofstrMggsa("output_data/direct_mggsa_operational_characteristics_mggsa.txt");
    if (!ofstrMggsa.is_open()) cerr << "File opening error\n";

    const int chunkMggsa = 2;

#pragma omp parallel for schedule(static, chunkMggsa) PROC_BIND num_threads(omp_get_num_procs()) collapse(2) \
        shared(numberFamily, problems, r, sizeR, successRate, K, key, d) \
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
                    numberTrialsArray[k] = numberTrials;
                } else {
                    numberTrialsArray[k] = K[i][1] + 1;
                }
            }
            
            for (int k = K[i][0]; k <= K[i][1]; k += K[i][2]) {
                numberSuccessful = (int)count_if(numberTrialsArray.begin(), numberTrialsArray.end(), [k](double elem){ return elem <= k; });
                successRate[i][j + 2][k / K[i][2]] = (double)numberSuccessful / numberFunctions;
            }
            double endTime = omp_get_wtime();
            double workTime = endTime - startTime;

            string strOutput = "MGGSA: " + problems[i].name + " r = " + to_string(r[i][j]) + " key = " + to_string(key[j]) + 
                               " time: " + to_string(workTime) + " thread num: " + to_string(omp_get_thread_num()) + "\n";
            cout << strOutput;
        }
    }

    for (int i = 0; i < numberFamily; i++) {
        for (int j = 2; j < sizeR + 2; j++) {
            for (int k = K[i][0]; k <= K[i][1]; k += K[i][2]) {
                ofstrMggsa << k << " " << successRate[i][j][k / K[i][2]] << endl;
            }
            ofstrMggsa << endl << endl;
        }
    }
    ofstrMggsa.close();
#endif
    double totalEndTime = omp_get_wtime();
    double totalWorkTime = totalEndTime - totalStartTime;
    cout << "Total time: " << totalWorkTime << endl;

    int sizeKey = key.size();
    setVariableGnuplot(ofstrOpt, "numberKey", to_string(sizeKey), false);
    initArrayGnuplot(ofstrOpt, "familyName", numberFamily);
    initArrayGnuplot(ofstrOpt, "r", sizeKey * numberFamily);
    initArrayGnuplot(ofstrOpt, "key", sizeKey);
    for (int i = 0; i < sizeKey; i++) {
        setValueInArrayGnuplot(ofstrOpt, "key", i + 1, key[i]);
    }
    for (int i = 0; i < numberFamily; i++) {
        setValueInArrayGnuplot(ofstrOpt, "familyName", i + 1, problems[i].shortName);

        for (int j = 0; j < r[i].size(); j++) {
            setValueInArrayGnuplot(ofstrOpt, "r", (i * sizeKey) + j + 1, r[i][j]);
        }
    }
    ofstrOpt.close();

    vector<int> args{ displayType, familyNumber };
    drawGraphGnuplot("scripts/direct_mggsa_operational_characteristics.gp", args);

#if defined( _MSC_VER )
    cin.get();
#endif
 */
	return 0;
}
