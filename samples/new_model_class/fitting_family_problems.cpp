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
#include <functional>
#include <algorithm>
#include <memory>
#include <limits>

#include <opt_methods/MggsaMethod.h>
#include <opt_methods/GsaMethod.h>
#include <test_opt_problems/FittingFamilyOptProblems.h>
#include <gnuplot/OutputFile.h>
#include <gnuplot/Script.h>
#include <MyMath.h>
#include <omp.h>

// #define SEARCH_OPTIMUM_MGGSA
// #define SEARCH_OPTIMUM_SCAN
// #define SLICES
// #define CALC_EXPERIMENT
#define CALC_COEFFICIENTS
// #define EXAMPLE
#define DRAW

using OptMethod = GsaMethod<OneDimensionalSupportiveOptProblem>;
using MggsaParameters = MggsaMethod<FittingFamilyOptProblems<OptMethod>>::Parameters;
using Result = MggsaMethod<FittingFamilyOptProblems<OptMethod>>::GeneralNumericalMethod::Result;
using Report = MggsaMethod<FittingFamilyOptProblems<OptMethod>>::Report;
using GsaParameters = GsaMethod<OneDimensionalSupportiveOptProblem>::Parameters;
using TypeSolve = MggsaMethod<FittingFamilyOptProblems<OptMethod>>::TypeSolve;

const std::string methodName = "mggsa";

const size_t displayType = 2; // 0 - application(only u(t)), 1 - png, 2 - png(notitle)
const size_t graphType = 0; // 0 - u(t), 1 - slices
const size_t problemNumber = 0; // 0, 1, ..., familySize - 1

int main() {
    // Default parameters //
    double accuracyGsa = 0.001, reliabilityGsa = 3.0;
    size_t maxTrialsGsa = 10000, maxFevalsGsa = 10000;
    GsaParameters gsaParameters(accuracyGsa, 0.0, maxTrialsGsa, maxFevalsGsa, reliabilityGsa);

    GsaMethod<OneDimensionalSupportiveOptProblem> gsa;
    gsa.setParameters(gsaParameters);

    FittingFamilyOptProblems fittingFamilyOptProblems(gsa);

    size_t dimension = fittingFamilyOptProblems.getDimension();
    size_t familySize = fittingFamilyOptProblems.getFamilySize();

    double accuracy = 0.01, error = 0.03, d = 0.01;
    std::vector<double> reliability(dimension + 4, 2.0);
    size_t maxTrials = 10000, maxFevals = 1000000000;
    size_t density = 12, key = 1, increment = 0;
    TypeSolve typeSolve = TypeSolve::SOLVE;
    MggsaParameters mggsaParameters(accuracy, error, maxTrials, maxFevals, reliability,
                                    d, density, key, increment, typeSolve);

    MggsaMethod<FittingFamilyOptProblems<OptMethod>> mggsa;
    mggsa.setParameters(mggsaParameters);
    mggsa.setProblem(fittingFamilyOptProblems);

    int chunk = 1;
    double step = 0.01;

    opt::MultiDimensionalSearchArea area;
    // Default parameters //

    double startTime, endTime, workTime;

    double totalStartTime = omp_get_wtime();
#if defined( CALC_EXPERIMENT )
    std::vector<std::vector<double>> reliabilities
    {
        std::vector<double>(dimension + 4, 10.0),
        std::vector<double>(dimension + 4, 7.0),
        std::vector<double>(dimension + 4, 5.0),
        std::vector<double>(dimension + 4, 3.0),
    };
    size_t reliabilitySize = reliabilities.size();

    std::vector<std::vector<size_t>> numberTrials(reliabilitySize, std::vector<size_t>(familySize));
    std::vector<std::vector<bool>> reached(reliabilitySize, std::vector<bool>(familySize));

    mggsaParameters = MggsaParameters(0.0, error, maxTrials, maxFevals, reliabilities[0],
                                      d, density, key, increment, typeSolve);

    mggsa.setParameters(mggsaParameters);

#pragma omp parallel for schedule(dynamic, chunk) PROC_BIND num_threads(1) collapse(2) \
        shared(reliabilities, reliabilitySize, numberTrials, reached) firstprivate(fittingFamilyOptProblems, mggsa)
    for (size_t i = 0; i < reliabilitySize; ++i) {
        for (size_t j = 0; j < familySize; ++j) {
            mggsa.setReliability(reliabilities[i]);
            fittingFamilyOptProblems.setProblemNumber(j);
            mggsa.setProblem(fittingFamilyOptProblems);

            std::unique_ptr<Result> result(static_cast<Result*>(mggsa.createResult()));
            std::unique_ptr<Report> report(static_cast<Report*>(mggsa.createReport()));

            double startTime = omp_get_wtime();
            reached[i][j] = mggsa.solveTest(*result);
            double endTime = omp_get_wtime();
            double workTime = endTime - startTime;

            numberTrials[i][j] = result->numberTrials;

            std::ostringstream output;
            output << "MggsaMethod, problem number = " << j << ", reliability = " << reliabilities[i][0]
                   << ", number trials = " << result->numberTrials << ", point = ";
            report->printPoint(output, result->point);
            output << ", value = " << result->value << ", work time = "
                   << workTime << ", thread number = " << omp_get_thread_num() << " " << "\n";
            std::cout << output.str();
        }
    }

    std::stringstream fileName;
    for (size_t i = 0; i < reliabilitySize; ++i) {
        fileName << std::setprecision(2) << "output_data/fitting_family_problems/experiment/"
                 << reliabilities[i][0] << ".txt";
        std::ofstream experimentFile(fileName.str());
        if (!experimentFile.is_open()) std::cerr << fileName.str() + " opening error\n";
        fileName.str("");

        for (size_t j = 0; j < familySize; ++j) {
            experimentFile << "problem number = " << j
                           << ", number trials = " << numberTrials[i][j]
                           << ", reached = " << reached[i][j]
                           << "\n";
        }
        experimentFile.close();
    }
#endif
    double totalEndTime = omp_get_wtime();
    std::cout << "Total time: " << totalEndTime - totalStartTime << "\n";

#if defined( DRAW )
    Script script("scripts/fitting_family_problems.gp");
    script.addArgs(std::vector<size_t>{ displayType, graphType, problemNumber });
    script.start();
    if (script.isError() == 2) std::cerr << "Error gnuplot\n";
    if (script.isError() == 1) std::cerr << "Error chmod\n";
#endif

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}
