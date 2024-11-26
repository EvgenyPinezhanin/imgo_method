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
#include <test_opt_problems/FittingFamilyOptProblems.h>
#include <gnuplot/OutputFile.h>
#include <gnuplot/Script.h>
#include <MyMath.h>
#include <omp.h>

// #define SLICES
#define SOLUTION
#define DRAW

using OptMethod = GsaMethod<OneDimensionalSupportiveOptProblem>;
using MggsaParameters = MggsaMethod<FittingFamilyOptProblems<OptMethod>>::Parameters;
using Result = MggsaMethod<FittingFamilyOptProblems<OptMethod>>::GeneralNumericalMethod::Result;
using Report = MggsaMethod<FittingFamilyOptProblems<OptMethod>>::Report;
using GsaParameters = GsaMethod<OneDimensionalSupportiveOptProblem>::Parameters;
using TypeSolve = MggsaMethod<FittingFamilyOptProblems<OptMethod>>::TypeSolve;

const std::string methodName = "mggsa";

const size_t displayType = 1; // 0 - application(only u(t)), 1 - png, 2 - png(notitle)
const size_t graphType = 1; // 0 - slices, 1 - solution, 2 - multi curves
const size_t problemNumber = 2; // 0, 1, ..., familySize - 1 (only for application type)

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
    size_t numberCoefficients = fittingFamilyOptProblems.getNumberCoefficients();

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
#if defined( SLICES )
    step = 0.01;

    area = fittingFamilyOptProblems.getSearchArea();

    int N = floor((area.upBound[0] - area.lowerBound[0]) / step);
    double eps = step / 2.0;
    size_t numberConstraints = fittingFamilyOptProblems.getNumberConstraints();

    familySize = 1;

    for (size_t i = 0; i < familySize; ++i) {
        fittingFamilyOptProblems.setProblemNumber(i);
        for (size_t first = 0; first < 3; ++first) {
            for (size_t second = first + 1; second < 4; ++second) {
                std::ofstream ofstrFunction("output_data/fitting_family_problems/slices/" + std::to_string(i + 1) + "_" +
                                            std::to_string(first + 1) + "_" + std::to_string(second + 1) + "f.txt");
                std::ofstream ofstrConstraints("output_data/fitting_family_problems/slices/" + std::to_string(i + 1) + "_" +
                                               std::to_string(first + 1) + "_" + std::to_string(second + 1) + "g.txt");

                std::vector<std::vector<double>> optimalPoints;
                fittingFamilyOptProblems.getOptimalPoints(optimalPoints);
                std::vector<double> optimalPoint = optimalPoints[0];

                // TODO: fix border
                ofstrFunction << N + 1 << " ";
                for (double j = area.lowerBound[0]; j <= area.upBound[0] + eps; j += step) {
                    ofstrFunction << j << " ";
                }
                ofstrFunction << "\n";
                for (double j = area.lowerBound[0]; j <= area.upBound[0] + eps; j += step) {
	                ofstrFunction << j << " ";
                    optimalPoint[first] = j;
                    for (double k = area.lowerBound[0]; k <= area.upBound[0] + eps; k += step) {
                        optimalPoint[second] = k;
	                    ofstrFunction << fittingFamilyOptProblems.computeObjectiveFunction(optimalPoint) << " ";
                    }
                    ofstrFunction << "\n";
                }
                ofstrFunction << std::endl;
                ofstrFunction.close();

                bool constrained;

                ofstrConstraints << N + 1 << " ";
                for (double j = area.lowerBound[0]; j <= area.upBound[0] + eps; j += step) {
                    ofstrConstraints << j << " ";
                }
                ofstrConstraints << "\n";
                for (double j = area.lowerBound[0]; j <= area.upBound[0] + eps; j += step) {
	                ofstrConstraints << j << " ";
                    optimalPoint[first] = j;
                    for (double k = area.lowerBound[0]; k <= area.upBound[0] + eps; k += step) {
                        optimalPoint[second] = k;
                        constrained = true;
                        for (size_t l = 0; l < numberConstraints; ++l) {
                            if (fittingFamilyOptProblems.computeConstraintFunction(optimalPoint, l) > 0.0) {
                                constrained = false;
                            }
                        }
                        if (constrained) {
	                        ofstrConstraints << fittingFamilyOptProblems.computeObjectiveFunction(optimalPoint) << " ";
                        } else {
                            ofstrConstraints << fittingFamilyOptProblems.getOptimalValue() - 1.0 << " ";
                        }
                    }
                    ofstrConstraints << "\n";
                }
                ofstrConstraints << std::endl;
                ofstrConstraints.close();
            }
        }
    }
#endif

#if defined( SOLUTION )
    OutputFile varsSolutionFile("output_data/new_model_class/fitting_family_problems_plot/vars_solution.txt");
    if (!varsSolutionFile.isOpen()) std::cerr << "vars_answer.txt opening error\n";

    varsSolutionFile.setVariable("A", fittingFamilyOptProblems.getLeftBoundWideWindow(), false);
    varsSolutionFile.setVariable("B", fittingFamilyOptProblems.getRightBoundWideWindow(), false);
    varsSolutionFile.setVariable("delta", fittingFamilyOptProblems.getDelta(), false);
    varsSolutionFile.setVariable("numberCoefficients", fittingFamilyOptProblems.getNumberCoefficients(), false);
    varsSolutionFile.setVariable("familySize", familySize, false);

    varsSolutionFile.initArray("coeffs", numberCoefficients * familySize);
    varsSolutionFile.initArray("x_opt", dimension * familySize);
    varsSolutionFile.initArray("T", 3 * familySize);
    varsSolutionFile.initArray("Q", 3 * familySize);

    std::vector<std::vector<double>> optimalPoints;
    std::vector<double> coefficients;
    std::vector<point> testPoints;

    fittingFamilyOptProblems.getTestPoints(testPoints);
    varsSolutionFile.setVariable("left_bound", fittingFamilyOptProblems.getLeftBoundWideWindow(), false);
    varsSolutionFile.setVariable("right_bound", testPoints[testPoints.size() - 1].x[0], false);

    size_t numberWideWindowPoints = fittingFamilyOptProblems.getNumberWideWindowPoints();

    for (size_t i = 3; i < 4; ++i) {
        fittingFamilyOptProblems.setProblemNumber(i);

        std::cout << i << "\n";
        fittingFamilyOptProblems.getOptimalPoints(optimalPoints);
        fittingFamilyOptProblems.computeObjectiveFunction(optimalPoints[0]);
        fittingFamilyOptProblems.getCoefficients(coefficients);

        for (size_t k = 0; k < numberCoefficients; ++k) {
            std::cout << coefficients[k] << " ";
        }
        std::cout << "\n";

        varsSolutionFile.setValuesInArray("coeffs", i * numberCoefficients + 1, coefficients, false);
        varsSolutionFile.setValuesInArray("x_opt", i * dimension + 1, optimalPoints[0], false);

        fittingFamilyOptProblems.getTestPoints(testPoints);

        for (size_t j = 0; j < 3; ++j) {
            varsSolutionFile.setValueInArray("T", i * 3 + j + 1, testPoints[numberWideWindowPoints + j].x[0], false);
            varsSolutionFile.setValueInArray("Q", i * 3 + j + 1, testPoints[numberWideWindowPoints + j].x[1], false);
        }
    }

    varsSolutionFile.close();
#endif

#if defined( MULTI_CURVES )
    OutputFile varsFile("output_data/fitting_family_problems/vars.txt");
    if (!varsFile.isOpen()) std::cerr << "vars.txt opening error\n";

    size_t numberCurves = 4;
    std::vector<size_t> vectorMaxTrials{ 500, 1000, 5000, 10000 };

    size_t problemNumber = 6;
    // TODO: think about set problem as const ref
    fittingFamilyOptProblems.setProblemNumber(problemNumber);
    mggsa.setProblem(fittingFamilyOptProblems);

    std::vector<double> coefficients;

    size_t numberCoefficients = fittingFamilyOptProblems.getNumberCoefficients();

    varsFile.setVariable("A", fittingFamilyOptProblems.getLeftBoundWideWindow(), false);
    varsFile.setVariable("B", fittingFamilyOptProblems.getRightBoundWideWindow(), false);
    varsFile.setVariable("delta", fittingFamilyOptProblems.getDelta(), false);
    varsFile.setVariable("numberCoefficients", numberCoefficients, false);

    varsFile.initArray("coeffs", numberCoefficients * numberCurves);
    varsFile.initArray("x_opt", dimension * numberCurves);

    std::unique_ptr<Result> result(mggsa.createResult());

    for (size_t i = 0; i < numberCurves; ++i) {
        mggsa.setMaxTrials(vectorMaxTrials[i]);

        mggsa.solve(*result);

        fittingFamilyOptProblems.computeObjectiveFunction(result->point);
        fittingFamilyOptProblems.getCoefficients(coefficients);

        for (size_t i = 0; i < fittingFamilyOptProblems.getNumberConstraints() + 1; ++i) {
            std::cout << i << " " << fittingFamilyOptProblems.computeConstraintFunction(result->point, i) << "\n";
        }
        std::cout << "\n";

        varsFile.setValuesInArray("coeffs", i * numberCoefficients + 1, coefficients, false);
        varsFile.setValuesInArray("x_opt", i * dimension + 1, result->point, false);
    }

    std::vector<point> testPoints;
    size_t numberWideWindowPoints = fittingFamilyOptProblems.getNumberWideWindowPoints();
    fittingFamilyOptProblems.getTestPoints(testPoints);

    varsFile.initArray("T", 3);
    varsFile.initArray("Q", 3);
    for (size_t i = 0; i < 3; ++i) {
        varsFile.setValueInArray("T", i + 1, testPoints[numberWideWindowPoints + i].x[0], false);
        varsFile.setValueInArray("Q", i + 1, testPoints[numberWideWindowPoints + i].x[1], false);
    }

    varsFile.close();
#endif
    double totalEndTime = omp_get_wtime();
    std::cout << "Total time: " << totalEndTime - totalStartTime << "\n";

#if defined( DRAW )
    Script script("scripts/new_model_class/fitting_family_problems_plot.gp");
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
