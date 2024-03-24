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

#define SEARCH_OPTIMUM_MGGSA
// #define SLICES
// #define SEARCH_OPTIMUM_SCAN
// #define DRAW

using MggsaParameters = MggsaMethod<FittingFamilyOptProblems>::Parameters;
using Result = MggsaMethod<FittingFamilyOptProblems>::Result;
using Report = MggsaMethod<FittingFamilyOptProblems>::Report;
using GsaParameters = GsaMethod<OneDimensionalSupportiveOptProblem>::Parameters;
using TypeSolve = MggsaMethod<FittingFamilyOptProblems>::TypeSolve;

const std::string methodName = "mggsa";

const size_t displayType = 1; // 0 - application(only u(t)), 1 - png, 2 - png(notitle)
const size_t graphType = 0; // 0 - u(t), 1 - slices
const size_t problemNumber = 0; // 0, 1, ..., familySize - 1

int main() {
    OutputFile varsFile("output_data/fitting_family_problems/vars.txt");
    if (!varsFile.isOpen()) std::cerr << "vars.txt opening error\n";

    double accuracyGsa = 0.01, reliabilityGsa = 2.5;
    size_t maxTrialsGsa = 10000, maxFevalsGsa = 10000;
    GsaParameters gsaParameters(accuracyGsa, 0.0, maxTrialsGsa, maxFevalsGsa, reliabilityGsa);

    GsaMethod<OneDimensionalSupportiveOptProblem> gsa;
    gsa.setParameters(gsaParameters);

    FittingFamilyOptProblems fittingFamilyOptProblems(gsa);

    size_t dimension = fittingFamilyOptProblems.getDimension();
    size_t familySize = fittingFamilyOptProblems.getFamilySize();

    double accuracy = 0.01, d = 0.01;
    std::vector<double> reliability(dimension + 4, 10.0);
    size_t maxTrials = 100000, maxFevals = 100000;
    size_t density = 11, key = 1, increment = 0;
    TypeSolve typeSolve = TypeSolve::SOLVE;
    MggsaParameters mggsaParameters(accuracy, 0.0, maxTrials, maxFevals, reliability,
                                    d, density, key, increment, typeSolve);

    MggsaMethod<FittingFamilyOptProblems> mggsa;
    mggsa.setParameters(mggsaParameters);

    std::vector<std::vector<double>> optPoints(familySize);
    std::vector<double> optValues(familySize);
    double startTime, endTime, workTime;

    varsFile.setVariable("leftBoundWindow", fittingFamilyOptProblems.getLeftBoundWindow(), false);
    varsFile.setVariable("rightBoundWindow", fittingFamilyOptProblems.getRightBoundWindow(), false);
    varsFile.setVariable("delta", fittingFamilyOptProblems.getDelta(), false);
    varsFile.setVariable("numberCoefficients", fittingFamilyOptProblems.getNumberCoefficients(), false);
    varsFile.initArray("optPoints", familySize * dimension);
    varsFile.initArray("optValues", familySize);
    for (size_t i = 0; i < familySize; ++i) {
        fittingFamilyOptProblems.setProblemNumber(i);
        varsFile.setValueInArray("optValues", i + 1, fittingFamilyOptProblems.getOptimalValue(), false);
    }

    double totalStartTime = omp_get_wtime();
#if defined( SEARCH_OPTIMUM_MGGSA )
    int chunk = 1;
#pragma omp parallel for schedule(dynamic, chunk) PROC_BIND num_threads(omp_get_num_procs()) \
        shared(optPoints, optValues) firstprivate(fittingFamilyOptProblems, mggsa)
    for (size_t i = 0; i < familySize; ++i) {
        fittingFamilyOptProblems.setProblemNumber(i);
        mggsa.setProblem(fittingFamilyOptProblems);

        std::unique_ptr<Result> result(static_cast<Result*>(mggsa.createResult()));
        std::unique_ptr<Report> report(static_cast<Report*>(mggsa.createReport()));

        double startTime = omp_get_wtime();
        mggsa.solve(*result);
        double endTime = omp_get_wtime();
        double workTime = endTime - startTime;

        optPoints[i] = result->point;
        optValues[i] = result->value;

        std::ostringstream output;
        output << "MggsaMethod, problem number = " << i << ", number trials = " << result->numberTrials
               << ", work time = " << workTime << ", thread number = " << omp_get_thread_num() << "\n";
        std::cout << output.str();
    }

    std::ofstream optPointsFile("output_data/fitting_family_problems/opt_points_mggsa.txt");
    if (!optPointsFile.is_open()) std::cerr << "opt_points_mggsa.txt opening error\n";
    for (size_t i = 0; i < familySize; ++i) {
        if (i % 10 == 0) {
            optPointsFile << "\n";
        }
        optPointsFile << "std::vector<std::vector<double>>{ std::vector<double>{ ";
        optPointsFile << optPoints[i][0];
        for (size_t j = 1; j < dimension; ++j) {
            optPointsFile << ", " << optPoints[i][j];
        }
        optPointsFile << " } },\n";
    }
    optPointsFile.close();

    std::ofstream optValuesFile("output_data/fitting_family_problems/opt_values_mggsa.txt");
    if (!optValuesFile.is_open()) std::cerr << "opt_values_mggsa.txt opening error\n";
    for (size_t i = 0; i < familySize; ++i) {
        if (i % 5 == 0) {
            optValuesFile << "\n";
        }
        optValuesFile << optValues[i] << ", ";
    }
    optValuesFile.close();
#endif

#if defined( SLICES )
    const double step = 0.01;

    opt::MultiDimensionalSearchArea area = fittingFamilyOptProblems.getSearchArea();

    int N = floor((area.upBound[0] - area.lowerBound[0]) / step);
    double eps = step / 2.0;
    size_t numberConstraints = fittingFamilyOptProblems.getNumberConstraints();

    for (size_t i = 0; i < familySize; ++i) {
        fittingFamilyOptProblems.setProblemNumber(i);
        for (size_t first = 0; first < 3; ++first) {
            for (size_t second = first + 1; second < 4; ++second) {
                std::ofstream ofstrFunction("output_data/fitting_family_problems/" + std::to_string(i + 1) + "_" +
                                            std::to_string(first + 1) + "_" + std::to_string(second + 1) + "f.txt");
                std::ofstream ofstrConstraints("output_data/fitting_family_problems/" + std::to_string(i + 1) + "_" +
                                               std::to_string(first + 1) + "_" + std::to_string(second + 1) + "g.txt");

                std::vector<std::vector<double>> optimalPoints;
                fittingFamilyOptProblems.getOptimalPoints(optimalPoints);
                std::vector<double> optimalPoint = optimalPoints[0];

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

#if defined( SEARCH_OPTIMUM_SCAN )
    double step = 0.2;

    std::vector<std::vector<std::vector<double>>> optimalPoints(familySize, std::vector<std::vector<double>>(10, std::vector<double>(4)));
    std::vector<std::vector<double>> optimalValues(familySize, std::vector<double>(10));

    int chunk = 1;
#pragma omp parallel for schedule(dynamic, chunk) PROC_BIND num_threads(omp_get_num_procs()) \
        shared(optimalPoints, optimalValues) firstprivate(fittingFamilyOptProblems, step, familySize)
    for (size_t number = 0; number < familySize; ++number) { 
        double minValue = std::numeric_limits<double>::infinity(), tmpValue;
        std::vector<double> minPoint(4), point(4);
        bool constrained;
        size_t count = 0;

        fittingFamilyOptProblems.setProblemNumber(number);
        opt::MultiDimensionalSearchArea area = fittingFamilyOptProblems.getSearchArea();

        double startTime = omp_get_wtime();
        for (double i = area.lowerBound[0]; i <= area.upBound[0]; i+= step) {
            for (double j = area.lowerBound[1]; j <= area.upBound[1]; j+= step) {
                for (double k = area.lowerBound[2]; k <= area.upBound[2]; k+= step) {
                    for (double l = area.lowerBound[3]; l <= area.upBound[3]; l+= step) {
                        point = std::vector<double>{ i, j, k, l };
                        constrained = true;

                        for (size_t constraint = 0; constraint < 7; ++constraint) {
                            if (fittingFamilyOptProblems.computeConstraintFunction(point, constraint) > 0.0) {
                                constrained = false;
                                break;
                            }
                        }

                        if (constrained) {
                            tmpValue = fittingFamilyOptProblems.computeObjectiveFunction(point);
                            if (tmpValue <= minValue) {
                                minPoint = point;
                                minValue = tmpValue;
                                optimalPoints[number][count % 10] = minPoint;
                                optimalValues[number][count % 10] = minValue;
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

    std::ofstream optPointsFileScan("output_data/fitting_family_problems/opt_points_scan.txt");
    if (!optPointsFileScan.is_open()) std::cerr << "opt_points_scan.txt opening error\n";
    for (size_t i = 0; i < familySize; ++i) {
        optPointsFileScan << "problem number = " << i << "\n"; 
        for (size_t j = 0; j < 10; ++j) {
            optPointsFileScan << "std::vector<std::vector<double>>{ std::vector<double>{ ";
            optPointsFileScan << optimalPoints[i][j][0];
            for (size_t k = 1; k < dimension; ++k) {
                // std::cout << i << " " << j << " " << k << "\n";
                optPointsFileScan << ", " << optimalPoints[i][j][k];
            }
            optPointsFileScan << " } },\n";
        }
        optPointsFileScan << "\n";
    }
    optPointsFileScan.close();

    std::ofstream optValuesFileScan("output_data/fitting_family_problems/opt_values_scan.txt");
    if (!optValuesFileScan.is_open()) std::cerr << "opt_values_scan.txt opening error\n";
    for (size_t i = 0; i < familySize; ++i) {
        optValuesFileScan << "problem number = " << i << "\n"; 
        for (size_t j = 0; j < 10; ++j) {
            optValuesFileScan << optimalValues[i][j] << ", ";
        }
        optValuesFileScan << "\n";
    }
    optValuesFileScan.close();

#endif
    double totalEndTime = omp_get_wtime();
    std::cout << "Total time: " << totalEndTime - totalStartTime << "\n";

#if defined( DRAW )
    Script script("scripts/fitting_family_problems.gp");
    script.addArgs(std::vector<int>{ displayType, graphType, problemNumber });
    script.start();
    if (script.isError() == 2) std::cerr << "Error gnuplot\n";
    if (script.isError() == 1) std::cerr << "Error chmod\n";
#endif

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}
        // trialsFile.open("output_data/sample_fitting_problem/" + method_names[i] + "_trials.txt");
        // if (!trialsFile.isOpen()) std::cerr << "trialsFile opening error\n";


        // report->print(std::cout, task, *result, workTime);

        // sample_test_problem.computeObjectiveFunction(result->point);
        // std::cout << "c = (" << coefficients[0];
        // for (int i = 1; i < number_coefficients; ++i) {
        //     std::cout << "; " << coefficients[i];
        // }
        // std::cout << ")\n\n";
// 
        // vars_file.setValueInArray("methodNames", i + 1, method_names[i]);
        // vars_file.setValueInArray("xOpt", i + 1, result->point, false);
        // vars_file.setValuesInArray("c", i * number_coefficients + 1, coefficients, false);
// 
        // sample_test_problem.getOptimalPoints(optimal_points);
        // trialsFile.addPoints(optimal_points, sample_test_problem.getOptimalValue());
        // trialsFile.addPoint(result->point, result->value);
// 
        // methods[i]->getTrialPoints(trials);
        // trialsFile.addPoints(trials);
// 
        // trialsFile.close();

/*             int numberTrials, numberFevals;
    std::vector<double> X, L;

    trialsFile.open("output_data/sample_test_problem_family/" + methodName + "_trials.txt");
    if (!trialsFile.is_open()) std::cerr << "trialsFile opening error\n";

    start_time = omp_get_wtime();
    mggsa.solve(numberTrials, numberFevals, X);
    mggsa.getL(L);
    end_time = omp_get_wtime();
    work_time = end_time - start_time;

    std::cout << "Method name: " << methodName << '\n';
    printResultMggsa("sample_test_problem_family", N + 1, 2 * N + 1, A, B, L, std::vector<double>{0, 0, 0, 0}, -1,
                     maxTrials, maxFevals, accuracy, reliability, d, den, key, incr, numberTrials, numberFevals, L, X,
                     problem(X, 2 * N + 1));
    problem(X, 3);
    std::cout << "max |u_der| = " << std::abs(problem(X, 2 * N + 1)) << std::endl; 
    varsFile.setVariable("minValue", problem(X, 2 * N + 1), false);
    functionGrid(problem, A, B, step, problem(X, 2 * N + 1), 2 * N + 1, X);

    problem(X, 3);
    std::cout << "c = (" << coefficients[0];
    for (int i = 1; i < number_coefficients; ++i) {
        std::cout << "; " << coefficients[i];
    }
    std::cout << ")\n\n";
    } */

    /*     varsFile.setValueInArray("method_names", 1, methodName);
    varsFile.initArray("x_opt", X, false);
    varsFile.initArray("c", coefficients, false);
    varsFile.initArray("T", 3);
    varsFile.initArray("Q", 3);
    for (int i = 0; i < 3; ++i) {
        varsFile.setValueInArray("T", i + 1, q[i + 7].x[0], false);
        varsFile.setValueInArray("Q", i + 1, q[i + 7].x[1], false);
    } */

    // sample_test_problem.getOptimalPoints(optimal_points);
    // trialsFile.add_points(optimal_points, sample_test_problem.getOptimalValue());
    // trialsFile.add_point(result.point, result.value);

    // methods[i]->getTrialPoints(trials);
    // trialsFile.add_points(trials);

    // trialsFile.close();

/*     varsFile.close();

    OutputFile test_points_file("output_data/sample_test_problem_family/test_points.txt");
    if (!test_points_file.isOpen()) std::cerr << "test_points_file opening error\n";

    for (int i = 0; i < number_window_points; ++i) {
        test_points_file.addPoint(q[i].x, q[i].x[1], false);
    }

    test_points_file.close(); */

    // output_file function_points_file("output_data/sample_test_problem/function_points.txt");
    // if (!function_points_file.is_open()) std::cerr << "function_points_file opening error\n";
// 
    // double a = sample_test_problem.getSearchArea().lowerBound;
    // double b = sample_test_problem.getSearchArea().upBound;
    // for (double i = a; i < b; i += step) {
    //     function_points_file.add_point(i, sample_test_problem.computeObjFunction(i), false);
    // }

    // function_points_file.close();