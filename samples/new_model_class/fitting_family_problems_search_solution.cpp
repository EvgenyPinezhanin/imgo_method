#if defined( _MSC_VER )
    #define PROC_BIND
#else
    #define PROC_BIND proc_bind(master)
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <memory>

#include <opt_methods/MggsaMethod.h>
#include <opt_methods/GsaMethod.h>
#include <opt_methods/ScanningMethod.h>
#include <test_opt_problems/FittingFamilyOptProblems.h>
#include <gnuplot/OutputFile.h>
#include <gnuplot/Script.h>
#include <MyMath.h>
#include <omp.h>

// #define MGGSA_WITH_GSA
#define MGGSA_WITH_SCANNING

#define SEARCH_OPTIMUM_MGGSA
// #define SEARCH_OPTIMUM_SCAN
#define DRAW

#if defined( MGGSA_WITH_GSA )
    using OptMethod = GsaMethod<OneDimensionalSupportiveOptProblem>;
#elif defined( MGGSA_WITH_SCANNING )
    using OptMethod = ScanningMethod<OneDimensionalSupportiveOptProblem>;
#endif

using MggsaParameters = MggsaMethod<FittingFamilyOptProblems<OptMethod>>::Parameters;
using Result = MggsaMethod<FittingFamilyOptProblems<OptMethod>>::GeneralNumericalMethod::Result;
using Report = MggsaMethod<FittingFamilyOptProblems<OptMethod>>::Report;
using OptMethodParameters = OptMethod::Parameters;
using TypeSolve = MggsaMethod<FittingFamilyOptProblems<OptMethod>>::TypeSolve;

const std::string methodName = "mggsa";

const size_t displayType = 2; // 0 - application(only u(t)), 1 - png, 2 - png(notitle)
const size_t graphType = 0; // 0 - u(t), 1 - slices
const size_t problemNumber = 0; // 0, 1, ..., familySize - 1

int main() {
    // Default parameters //
    double accuracyOptMethod = 0.000005, reliabilityOptMethod = 3.0;
    size_t maxTrialsOptMethod = 10000, maxFevalsOptMethod = 10000;

#if defined( MGGSA_WITH_GSA )
    OptMethodParameters optMethodParameters(accuracyOptMethod, 0.0, maxTrialsOptMethod,
                                            maxFevalsOptMethod, reliabilityOptMethod);
#elif defined( MGGSA_WITH_SCANNING )
    OptMethodParameters optMethodParameters(accuracyOptMethod, 0.0, maxTrialsOptMethod, maxFevalsOptMethod);
#endif

    OptMethod optMethod;
    optMethod.setParameters(optMethodParameters);

    FittingFamilyOptProblems fittingFamilyOptProblems(optMethod);

    size_t dimension = fittingFamilyOptProblems.getDimension();
    size_t familySize = fittingFamilyOptProblems.getFamilySize();

    double accuracy = 0.005, error = 0.0, d = 0.01;
    std::vector<double> reliability(dimension + 4, 3.0);
    size_t maxTrials = 60000, maxFevals = 1000000000;
    size_t density = 12, key = 1, increment = 0;
    TypeSolve typeSolve = TypeSolve::SOLVE;
    MggsaParameters mggsaParameters(accuracy, error, maxTrials, maxFevals, reliability,
                                    d, density, key, increment, typeSolve);

    MggsaMethod<FittingFamilyOptProblems<OptMethod>> mggsa;
    mggsa.setParameters(mggsaParameters);
    mggsa.setProblem(fittingFamilyOptProblems);

    int chunk = 2;
    double step = 0.01;

    opt::MultiDimensionalSearchArea area;
    // Default parameters //

    double startTime, endTime, workTime;

    double totalStartTime = omp_get_wtime();
#if defined( SEARCH_OPTIMUM_MGGSA )
    std::vector<std::vector<double>> optimalPointsMGGSA(familySize);
    std::vector<double> optimalValuesMGGSA(familySize);

#pragma omp parallel for schedule(dynamic, chunk) PROC_BIND num_threads(omp_get_num_procs()) \
        shared(optimalPointsMGGSA, optimalValuesMGGSA) firstprivate(fittingFamilyOptProblems, mggsa)
    for (size_t i = 0; i < familySize; ++i) {
        fittingFamilyOptProblems.setProblemNumber(i);
        mggsa.setProblem(fittingFamilyOptProblems);

        std::unique_ptr<Result> result(static_cast<Result*>(mggsa.createResult()));
        std::unique_ptr<Report> report(static_cast<Report*>(mggsa.createReport()));

        double startTime = omp_get_wtime();
        mggsa.solve(*result);
        double endTime = omp_get_wtime();
        double workTime = endTime - startTime;

        optimalPointsMGGSA[i] = result->point;
        optimalValuesMGGSA[i] = result->value;

        std::vector<std::vector<double>> optimalPoints;
        fittingFamilyOptProblems.getOptimalPoints(optimalPoints);
        double optimalValue = fittingFamilyOptProblems.getOptimalValue();

        std::ostringstream output;
        output << "MGGSAMethod, problem number = " << i << ", number trials = " << result->numberTrials << ", X = ";
        report->printPoint(output, result->point);
        output << ", f(X) = " << result->value << ", X_opt = ";
        report->printPoint(output, optimalPoints[0]);
        output << ", f(X_opt) = " << optimalValue << ", f(X_opt) - f(X) = " << optimalValue - result->value
               << ", work time = " << workTime << ", thread number = " << omp_get_thread_num() << "\n";
        std::cout << output.str();
    }

    std::ofstream optimalPointsFileMGGSA("output_data/new_model_class/fitting_family_problems_search_solution/opt_points_mggsa.txt");
    if (!optimalPointsFileMGGSA.is_open()) std::cerr << "opt_points_mggsa.txt opening error\n";
    for (size_t i = 0; i < familySize; ++i) {
        if (i % 10 == 0) {
            optimalPointsFileMGGSA << "\n";
        }
        optimalPointsFileMGGSA << "std::vector<std::vector<double>>{ std::vector<double>{ ";
        optimalPointsFileMGGSA << optimalPointsMGGSA[i][0];
        for (size_t j = 1; j < dimension; ++j) {
            optimalPointsFileMGGSA << ", " << optimalPointsMGGSA[i][j];
        }
        optimalPointsFileMGGSA << " } },\n";
    }
    optimalPointsFileMGGSA.close();

    std::ofstream optimalValuesFileMGGSA("output_data/new_model_class/fitting_family_problems_search_solution/opt_values_mggsa.txt");
    if (!optimalValuesFileMGGSA.is_open()) std::cerr << "opt_values_mggsa.txt opening error\n";
    for (size_t i = 0; i < familySize; ++i) {
        if (i % 5 == 0) {
            optimalValuesFileMGGSA << "\n";
        }
        optimalValuesFileMGGSA << optimalValuesMGGSA[i] << ", ";
    }
    optimalValuesFileMGGSA.close();
#endif

#if defined( SEARCH_OPTIMUM_SCAN )
    fittingFamilyOptProblems.setProblemNumber(0);
    area = fittingFamilyOptProblems.getSearchArea();
    step = (area.upBound[0] - area.lowerBound[0]) / 100;
    size_t startProblem = 0, finishProblem = 99;
    fittingFamilyOptProblems.setIsSortX(false);
    size_t numberOptimalPoints;

    std::vector<std::vector<std::vector<double>>> optimalPointsScan(familySize,
        std::vector<std::vector<double>>(numberOptimalPoints, std::vector<double>(dimension)));
    std::vector<std::vector<double>> optimalValuesScan(familySize, std::vector<double>(numberOptimalPoints));

#pragma omp parallel for schedule(dynamic, chunk) PROC_BIND num_threads(omp_get_num_procs()) \
        shared(optimalPointsScan, optimalValuesScan, startProblem, finishProblem) \
        firstprivate(fittingFamilyOptProblems, step, familySize, dimension)
    for (size_t number = startProblem; number < finishProblem + 1; ++number) { 
        double minValue = std::numeric_limits<double>::infinity(), tmpValue;
        std::vector<double> minPoint(dimension), currPoint(dimension);
        bool constrained;
        size_t count = 0;

        fittingFamilyOptProblems.setProblemNumber(number);
        opt::MultiDimensionalSearchArea area = fittingFamilyOptProblems.getSearchArea();
        size_t numberConstraints;

        double startTime = omp_get_wtime();
        for (double i = area.lowerBound[0]; i <= area.upBound[0]; i+= step) {
            for (double j = area.lowerBound[1]; j <= area.upBound[1]; j+= step) {
                for (double k = area.lowerBound[2]; k <= area.upBound[2]; k+= step) {
                    for (double l = area.lowerBound[3]; l <= area.upBound[3]; l+= step) {
                        currPoint = std::vector<double>{ i, j, k, l };
                        constrained = true;

                        numberConstraints = fittingFamilyOptProblems.getNumberConstraints();
                        for (size_t constraint = 0; constraint < numberConstraints; ++constraint) {
                            if (fittingFamilyOptProblems.computeConstraintFunction(currPoint, constraint) > 0.0) {
                                constrained = false;
                                break;
                            }
                        }

                        if (constrained) {
                            tmpValue = fittingFamilyOptProblems.computeObjectiveFunction(currPoint);
                            if (tmpValue <= minValue) {
                                minPoint = currPoint;
                                minValue = tmpValue;
                                optimalPointsScan[number][count % 10] = minPoint;
                                optimalValuesScan[number][count % 10] = minValue;
                                ++count;
                            }
                        }
                    }
                }
            }
        }
        double endTime = omp_get_wtime();
        double workTime = endTime - startTime;

        std::ostringstream output;
        output << "MggsaMethod, problem number = " << number << ", work time = " << workTime
               << ", thread number = " << omp_get_thread_num() << "\n";
        std::cout << output.str();
    }

    std::ofstream optimalPointsFileScan("output_data/fitting_family_problems/opt_points_scan.txt");
    if (!optimalPointsFileScan.is_open()) std::cerr << "opt_points_scan.txt opening error\n";
    for (size_t i = startProblem; i < finishProblem + 1; ++i) {
        optimalPointsFileScan << "problem number = " << i << "\n"; 
        for (size_t j = 0; j < 10; ++j) {
            optimalPointsFileScan << "std::vector<std::vector<double>>{ std::vector<double>{ ";
            optimalPointsFileScan << optimalPointsScan[i][j][0];
            for (size_t k = 1; k < dimension; ++k) {
                optimalPointsFileScan << ", " << optimalPointsScan[i][j][k];
            }
            optimalPointsFileScan << " } },\n";
        }
        optimalPointsFileScan << "\n";
    }
    optimalPointsFileScan.close();

    std::ofstream optValuesFileScan("output_data/fitting_family_problems/opt_values_scan.txt");
    if (!optValuesFileScan.is_open()) std::cerr << "opt_values_scan.txt opening error\n";
    for (size_t i = startProblem; i < finishProblem + 1; ++i) {
        optValuesFileScan << "problem number = " << i << "\n"; 
        for (size_t j = 0; j < 10; ++j) {
            optValuesFileScan << optimalValuesScan[i][j] << ", ";
        }
        optValuesFileScan << "\n";
    }
    optValuesFileScan.close();
#endif
    double totalEndTime = omp_get_wtime();
    std::cout << "Total time: " << totalEndTime - totalStartTime << "\n";

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}
