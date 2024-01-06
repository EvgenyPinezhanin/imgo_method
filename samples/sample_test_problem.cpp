#include <iostream>
#include <vector>
#include <string>
#include <functional>

#include <opt_methods/ScanningMethod.h>
#include <opt_methods/PiyavskyMethod.h>
#include <opt_methods/GsaMethod.h>
#include <opt_problems/OneDimensionalProblem.h>
#include <gnuplot/OutputFile.h>
#include <gnuplot/Script.h>
#include <MyMath.h>
#include <omp.h>

#define CALC
#define DRAW

const std::vector<std::string> method_names{ "scanning", "piyavsky", "gsa" };

const int display_type = 2; // 0 - application, 1 - png, 2 - png(notitle)
const int graph_type = 1; // 0 - u(t), 1 - trial points
const int method_index = 2; // 0 - scanning, 1 - piyavsky, 2 - gsa

const double step = 0.01;

size_t number_test_points = 17;
const std::vector<point> q{ point( 0.0, 0.74),
                            point( 1.0, 0.67),
                            point( 3.0, 1.1 ),
                            point( 3.9, 1.25),
                            point( 6.0, 1.2 ),
                            point( 8.4, 1.3 ),
                            point(10.0, 1.6 ),
                            point(12.0, 1.4 ),
                            point(15.0, 1.7 ),
                            point(17.0, 1.6 ),
                            point(20.0, 1.5 ),
                            point(22.0, 1.57),
                            point(23.0, 1.3 ),
                            point(24.0, 1.15),
                            point(24.5, 1.15),
                            point(26.2, 1.3 ),
                            point(27.0, 1.2 ) };

size_t number_coefficients = 5;
std::vector<double> coefficients(number_coefficients);
const OneDimensionalProblem sample_test_problem(
    [] (double x) -> double {
        std::vector<std::function<double(double, double)>> functions{
            [] (double t, double x) -> double { return 1.0; },
            [] (double t, double x) -> double { return t; },
            [] (double t, double x) -> double { return t * t; },
            [] (double t, double x) -> double { return std::sin(x * t); },
            [] (double t, double x) -> double { return std::cos(x * t); }
        };

        std::vector<std::vector<double>> A(number_coefficients,
                                           std::vector<double>(number_coefficients, 0));
        std::vector<double> b(number_coefficients, 0);

        for (size_t i = 0; i < number_coefficients; ++i) {
            for (size_t j = 0; j < number_coefficients; ++j) {
                for (size_t k = 0; k < number_test_points; ++k) {
                    A[i][j] += functions[j](q[k].x, x) * functions[i](q[k].x, x);
                }
            }
            for (size_t j = 0; j < number_test_points; ++j) {
                b[i] += q[j].y[0] * functions[i](q[j].x, x);
            }
        }

        mnk minimizer(A, b);
        minimizer.solve(coefficients);

        double result = 0.0, sum;
        for (size_t i = 0; i < number_test_points; ++i) {
            sum = 0.0;
            for (size_t j = 0; j < number_coefficients; ++j) {
                sum += coefficients[j] * functions[j](q[i].x, x);
            }
            result += (sum - q[i].y[0]) * (sum - q[i].y[0]);
        }

        return result;
    }, "Sample test problem", 1,
    opt::OneDimensionalSearchArea(0.5, 8.0), std::vector<double>{1.105}, 0.017, -1.0);

int main() {
    double accuracy = 0.001, reliability = 2.0;
    int maxTrials = 100000, maxFevals = 100000;
    GsaMethod<OneDimensionalProblem>::Parameters parameters(accuracy, 0.0, maxTrials, maxFevals, reliability);
    opt::Task task("Sample test task", sample_test_problem, parameters);

    size_t number_methods = 3;
    ScanningMethod<OneDimensionalProblem> scanning;
    PiyavskyMethod<OneDimensionalProblem> piyavsky;
    GsaMethod<OneDimensionalProblem> gsa;

    std::vector<ScanningMethod<OneDimensionalProblem>::GeneralNumericalMethod*> methods{ &scanning, &piyavsky, &gsa };

    ScanningMethod<OneDimensionalProblem>::Result scanning_result;
    PiyavskyMethod<OneDimensionalProblem>::Result piyavsky_result;
    GsaMethod<OneDimensionalProblem>::Result gsa_result;

    std::vector<ScanningMethod<OneDimensionalProblem>::GeneralNumericalMethod::Result*> method_results{ &scanning_result,
                                                                                                        &piyavsky_result,
                                                                                                        &gsa_result };

    ScanningMethod<OneDimensionalProblem>::Report scanning_report;
    PiyavskyMethod<OneDimensionalProblem>::Report piyavsky_report;
    GsaMethod<OneDimensionalProblem>::Report gsa_report;

    std::vector<ScanningMethod<OneDimensionalProblem>::GeneralNumericalMethod::IReport*> method_reports{ &scanning_report,
                                                                                                   &piyavsky_report,
                                                                                                   &gsa_report };

    double total_start_time = omp_get_wtime();
#if defined( CALC )
    OutputFile vars_file("output_data/sample_test_problem/vars.txt");
    if (!vars_file.isOpen()) std::cerr << "vars_file opening error\n";

    vars_file.setVariable("number_methods", number_methods, false);
    vars_file.setVariable("number_coefficients", number_coefficients, false);
    vars_file.initArray("method_names", number_methods);
    vars_file.initArray("x_opt", number_methods);
    vars_file.initArray("c", number_coefficients * number_methods);

    OutputFile trialsFile;
    std::vector<double> optimal_points;
    std::vector<Trial> trials;
    double start_time, end_time, work_time;

    for (size_t i = 0; i < number_methods; ++i) {
        trialsFile.open("output_data/sample_test_problem/" + method_names[i] + "_trials.txt");
        if (!trialsFile.isOpen()) std::cerr << "trialsFile opening error\n";

        methods[i]->setParameters(parameters);
        methods[i]->setProblem(sample_test_problem);

        start_time = omp_get_wtime();
        methods[i]->solve(*method_results[i]);
        end_time = omp_get_wtime();
        work_time = end_time - start_time;

        std::cout << "Method name: " << method_names[i] << '\n';
        method_reports[i]->print(std::cout, task, *method_results[i], work_time);

        sample_test_problem.computeObjFunction(method_results[i]->point);
        std::cout << "c = (" << coefficients[0];
        for (int i = 1; i < number_coefficients; ++i) {
            std::cout << "; " << coefficients[i];
        }
        std::cout << ")\n\n";

        vars_file.setValueInArray("method_names", i + 1, method_names[i]);
        vars_file.setValueInArray("x_opt", i + 1, method_results[i]->point, false);
        vars_file.setValuesInArray("c", i * number_coefficients + 1, coefficients, false);

        sample_test_problem.getOptimalPoints(optimal_points);
        trialsFile.addPoints(optimal_points, sample_test_problem.getOptimalValue());
        trialsFile.addPoint(method_results[i]->point, method_results[i]->value);

        methods[i]->getTrialPoints(trials);
        trialsFile.addPoints(trials);

        trialsFile.close();
    }

    vars_file.close();

    OutputFile test_points_file("output_data/sample_test_problem/test_points.txt");
    if (!test_points_file.isOpen()) std::cerr << "test_points_file opening error\n";

    for (int i = 0; i < number_test_points; ++i) {
        test_points_file.addPoint(q[i].x, q[i].y[0], false);
    }

    test_points_file.close();

    OutputFile function_points_file("output_data/sample_test_problem/function_points.txt");
    if (!function_points_file.isOpen()) std::cerr << "function_points_file opening error\n";

    double a = sample_test_problem.getSearchArea().lowerBound;
    double b = sample_test_problem.getSearchArea().upBound;
    for (double i = a; i < b; i += step) {
        function_points_file.addPoint(i, sample_test_problem.computeObjFunction(i), false);
    }

    function_points_file.close();
#endif
    double total_end_time = omp_get_wtime();
    std::cout << "Total time: " << total_end_time - total_start_time << "\n";

#if defined( DRAW )
    Script script("scripts/sample_test_problem.gp");
    script.addArgs(std::vector<int>{ display_type, graph_type, method_index });
    script.addArgs(std::vector<double>{ q[0].x, q[number_test_points - 1].x,
                                        sample_test_problem.getSearchArea().lowerBound,
                                        sample_test_problem.getSearchArea().upBound });
    script.start();
    if (script.isError() == 2) std::cerr << "Error gnuplot\n";
    if (script.isError() == 1) std::cerr << "Error chmod\n";
#endif

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}
