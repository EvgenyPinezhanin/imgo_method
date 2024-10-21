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

using MggsaParameters = MggsaMethod<FittingFamilyOptProblems>::Parameters;
using Result = MggsaMethod<FittingFamilyOptProblems>::GeneralNumericalMethod::Result;
using Report = MggsaMethod<FittingFamilyOptProblems>::Report;
using GsaParameters = GsaMethod<OneDimensionalSupportiveOptProblem>::Parameters;
using TypeSolve = MggsaMethod<FittingFamilyOptProblems>::TypeSolve;

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

    MggsaMethod<FittingFamilyOptProblems> mggsa;
    mggsa.setParameters(mggsaParameters);
    mggsa.setProblem(fittingFamilyOptProblems);

    int chunk = 1;
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

        std::ostringstream output;
        output << "MGGSAMethod, problem number = " << i << ", number trials = " << result->numberTrials
               << ", work time = " << workTime << ", thread number = " << omp_get_thread_num() << "\n";
        std::cout << output.str();
    }

    std::ofstream optimalPointsFileMGGSA("output_data/fitting_family_problems/opt_points_mggsa.txt");
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

    std::ofstream optimalValuesFileMGGSA("output_data/fitting_family_problems/opt_values_mggsa.txt");
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

#if defined( CALC_COEFFICIENTS )
    OutputFile varsFile("output_data/fitting_family_problems/vars.txt");
    if (!varsFile.isOpen()) std::cerr << "vars.txt opening error\n";

    size_t problemNumber = 0;
    fittingFamilyOptProblems.setProblemNumber(problemNumber);

    varsFile.setVariable("A", fittingFamilyOptProblems.getLeftBoundWideWindow(), false);
    varsFile.setVariable("B", fittingFamilyOptProblems.getRightBoundWideWindow(), false);
    varsFile.setVariable("delta", fittingFamilyOptProblems.getDelta(), false);
    varsFile.setVariable("numberCoefficients", fittingFamilyOptProblems.getNumberCoefficients(), false);

    std::vector<std::vector<double>> points {
        std::vector<double>{ 0.258507, 0.507743, 1.00524, 1.25351 },
        std::vector<double>{ 0.0865198, 0.483937, 0.947428, 1.04022 },
        std::vector<double>{ 0.0233606, 0.63649, 1.00524, 1.1782 },
        std::vector<double>{ 0.0238464, 0.633575, 1.06063, 1.17043 }
    };

    std::vector<std::vector<double>> coefficients_array {
        std::vector<double>{ -0.839262, 3.63611, -2.79527, -1.57913, -0.6276, -1.44098, 0.848332, -0.0322575 },
        std::vector<double>{ -3.97722, 1.12522, -0.706866, -2.31499, 2.1157, -4.21661, 0.163885, 3.73852 },
        std::vector<double>{ -22.6007, 3.11866, -2.45585, 1.6177, -3.35149, -1.44369, 1.79755, -1.21832 },
        std::vector<double>{ -24.0483, 3.53391, -2.30606, 2.01984, -5.65306, 0.916497, 2.92328, -3.13491 }
    };

    varsFile.initArray("coeffs", 8 * 4);
    varsFile.initArray("x_opt", 4 * 4);
    for (size_t i = 0; i < 4; i++) {
        varsFile.setValuesInArray("coeffs", i * 8 + 1, coefficients_array[i], false);
        varsFile.setValuesInArray("x_opt", i * 4 + 1, points[i], false);
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

#if defined( EXAMPLE )
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
