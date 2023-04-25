#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS    
#endif

#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>

#include <Grishagin/GrishaginProblemFamily.hpp>
#include <Grishagin/GrishaginConstrainedProblemFamily.hpp>
#include <GKLS/GKLSProblemFamily.hpp>
#include <GKLS/GKLSConstrainedProblemFamily.hpp>
#include <direct_method.h>
#include <task.h>
#include <output_results.h>
#include <omp.h>

using namespace std;

// #define CALC

const int familyNumber = 3; // 0 - Grishagin, 1 - GKLS,
                            // 2 - Grishagin(constrained), 3 - GKLS(constrained),
                            // 4 - comparison Grishagin and GKLS, 5 - comparison Grishagin and GKLS (constrained)

double f(int n, const double *X, int *undefinedFlag, void *data) {
    DataDirectOperationalCharacteristics *fData = static_cast<DataDirectOperationalCharacteristics*>(data);
    fData->numberFevals++;
    vector<double> point(n), optPoint(n);
    for (int i = 0; i < n; i++) {
        point[i] = X[i];
    }
    fData->points.push_back(point);

    double f;
    if (fData->type == TypeConstrants::Constraints) {
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
}

int main() {
    ofstream ofstrOpt("output_data/direct_operational_characteristics_opt.txt");
    if (!ofstrOpt.is_open()) cerr << "File opening error\n";

    const int numberFamily = 4;

    TGrishaginProblemFamily grishaginProblems;
    TGrishaginConstrainedProblemFamily grishaginConstrainedProblems;
    TGKLSProblemFamily gklsProblems;
    TGKLSConstrainedProblemFamily gklsConstrainedProblems;

    vector<ProblemFamily> problems{ ProblemFamily("GrishaginProblemFamily", &grishaginProblems, TypeConstrants::NoConstraints, "Grishagin"),
                                    ProblemFamily("GKLSProblemFamily", &gklsProblems, TypeConstrants::NoConstraints, "GKLS"),
                                    ProblemFamily("GrishaginProblemConstrainedFamily", &grishaginConstrainedProblems, 
                                                    TypeConstrants::Constraints, "GrishaginConstrained"),
                                    ProblemFamily("GKLSProblemConstrainedFamily", &gklsConstrainedProblems, 
                                                    TypeConstrants::Constraints, "GKLSConstrained") };

#if defined(CALC)
    ofstream ofstr("output_data/direct_operational_characteristics.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";

    vector<vector<int>> K{ { 0, 700, 25 },
                           { 0, 1600, 25 },
                           { 0, 500, 25 },
                           { 0, 4000, 25 } };

    DataDirectOperationalCharacteristics fData;

    vector<double> A, B;
    int maxIters = 10000;
    double magicEps = 1.0e-4, eps = 0.01;
    double volumeReltol = 0.0;
    double sigmaReltol  = 0.0;
    vector<direct_algorithm> algorithms{ DIRECT_ORIGINAL, DIRECT_GABLONSKY };

    DirectMethod direct(f, &fData, -1, A, B, -1, maxIters, magicEps, volumeReltol, sigmaReltol, nullptr, DIRECT_ORIGINAL);

    vector<vector<int>> numberFevalsArray(numberFamily);
    numberFevalsArray[0].resize(grishaginProblems.GetFamilySize(), 0);
    numberFevalsArray[1].resize(gklsProblems.GetFamilySize(), 0);
    numberFevalsArray[2].resize(grishaginConstrainedProblems.GetFamilySize(), 0);
    numberFevalsArray[3].resize(gklsConstrainedProblems.GetFamilySize(), 0);

    FunctorFamily functor;
    FunctorFamilyConstrained functorConstrained;
    int numberFunctions, numberSuccessful;
    double startTime, endTime, workTime;

    fData.eps = eps;

    double totalStartTime = omp_get_wtime();
    for (int i = 0; i < numberFamily; i++) {
        if (problems[i].type == TypeConstrants::Constraints) {
            functorConstrained.constrainedOptProblemFamily = static_cast<IConstrainedOptProblemFamily*>(problems[i].optProblemFamily);
            (*functorConstrained.constrainedOptProblemFamily)[0]->GetBounds(A, B);
            direct.setN((*functorConstrained.constrainedOptProblemFamily)[0]->GetDimension());
            fData.functor = &functorConstrained;
            fData.type = TypeConstrants::Constraints;
        } else {
            functor.optProblemFamily = static_cast<IOptProblemFamily*>(problems[i].optProblemFamily);
            (*functor.optProblemFamily)[0]->GetBounds(A, B);
            direct.setN((*functor.optProblemFamily)[0]->GetDimension());
            fData.functor = &functor;
            fData.type = TypeConstrants::NoConstraints;
        }
        direct.setAB(A, B);
        direct.setMaxFevals(K[i][1]);

        numberFunctions = problems[i].optProblemFamily->GetFamilySize();

        cout << problems[i].name << endl;
        for (int j = 0; j < algorithms.size(); j++) {
            cout << "Type of algorithm: " << ((algorithms[j] == DIRECT_ORIGINAL) ? "DIRECT_ORIGINAL" : "DIRECT_GABLONSKY") << endl;
            direct.setAlghorithm(algorithms[j]);

            startTime = omp_get_wtime();
            for (int k = 0; k < numberFunctions; k++) {
                fData.numberFevals = 0;
                fData.converge = false;

                if (problems[i].type == TypeConstrants::Constraints) {
                    functorConstrained.currentFunction = k;
                } else {
                    functor.currentFunction = k;
                }
                direct.solveTest();

                if (fData.converge) {
                    numberFevalsArray[i][k] = fData.minNumberFevals;
                } else {
                    numberFevalsArray[i][k] = K[i][1] + 1;
                }
            }
            for (int k = K[i][0]; k <= K[i][1]; k += K[i][2]) {
                numberSuccessful = (int)count_if(numberFevalsArray[i].begin(), numberFevalsArray[i].end(), [k] (double elem) { return elem <= k; });
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
    for (int i = 0; i < numberFamily; i++) {
        setValueInArrayGnuplot(ofstrOpt, "familyNames", i + 1, "\"" + problems[i].shortName + "\"");
    }
    ofstrOpt.close();

    drawGraphGnuplot("scripts/direct_operational_characteristics.gp", familyNumber);

#if defined(_MSC_VER)
    cin.get();
#endif

	return 0;
}
