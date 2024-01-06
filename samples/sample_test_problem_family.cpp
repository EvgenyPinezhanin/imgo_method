#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <algorithm>

#include <opt_methods/MggsaMethod.h>
#include <opt_methods/GsaMethod.h>
#include <opt_problems/MultiDimensionalConstrainedProblem.h>
#include <opt_problems/OneDimensionalProblem.h>
#include <gnuplot/OutputFile.h>
#include <gnuplot/Script.h>
#include <MyMath.h>
#include <omp.h>

#define CALC
#define DRAW

const std::string method_name = "mggsa";

const int display_type = 0; // 0 - application, 1 - png, 2 - png(notitle)
const int graph_type = 0; // 0 - u(t), (1 - trial points) - ?
const int method_index = 0; // 0 - mggsa

const double step = 0.01;

size_t N = 3;
size_t number_coefficients = 2 * N + 2;
std::vector<double> coefficients(number_coefficients);
// std::vector<double> coefficients{ 34.93, -2.151, -8.075, -7.936, 19.461, 1.048, -7.375, 5.411 };
std::vector<double> omega(N + 1);
double alpha = 0.05, delta = 0.3, a = 1.0, b = 10.0;
size_t number_window_points = 9;
const std::vector<point> q{ point( a,      0.0 ),
                            point( 3.0,    0.0 ),
                            point( 5.0,    0.0 ),
                            point( 7.0,    0.0 ),
                            point( 9.0,    0.0 ),
                            point( b,      0.0 ),
                            point( 13.0,   7.65 ),
                            point( 16.65, -9.86 ),
                            point( 18.0,   0.0  ) };

bool is_min;
const OneDimensionalProblem u(
    [] (double t) -> double {
        double result = 0.0;

        for (size_t i = 0; i < N + 1; ++i) {
            result += coefficients[2 * i] * std::sin(omega[i] * t) +
                      coefficients[2 * i + 1] * std::cos(omega[i] * t) ;
        }

        return is_min ? result : -result;
    },
    "Sample test problem family", 1, opt::OneDimensionalSearchArea(a, b));
double accuracy_in = 0.01, reliability_in = 2.5;
int maxTrials_in = 10000, maxFevals_in = 10000;
GsaMethod<OneDimensionalProblem>::Parameters parameters(accuracy_in, 0.0, maxTrials_in, maxFevals_in, reliability_in);
GsaMethod<OneDimensionalProblem>::Result result;

double u_der(std::vector<double> x) {
    double result = 0.0;
    for (size_t i = 0; i < N + 1; ++i) {
        result += coefficients[2 * i] * x[i] * std::cos(x[i] * q[number_window_points - 1].x) -
                  coefficients[2 * i + 1] * x[i] * std::sin(x[i] * q[number_window_points - 1].x) ;
    }
    return result;
}

// const MultiDimensionalConstrainedProblem sample_test_problem_family(
//     [] (std::vector<double> x, int index) -> double {
double problem(std::vector<double> x, int index) { 
        std::vector<std::vector<double>> A(number_coefficients,
                                           std::vector<double>(number_coefficients, 0));
        std::vector<double> b(number_coefficients, 0);

        std::vector<std::function<double(double, double)>> functions{
            [] (double t, double x) -> double { return std::sin(x * t); },
            [] (double t, double x) -> double { return std::cos(x * t); }
        };

        mnk minimizer;

        GsaMethod<OneDimensionalProblem> gsa(u, parameters);
        double min_value, max_value;

        switch (index) {
            case 0: case 1: case 2:
                std::sort(x.begin(), x.end());
                return x[index] * (1.0 + alpha) - x[index + 1] * (1.0 - alpha);

            case 3: case 4: case 5:
                for (size_t i = 0; i < number_coefficients; ++i) {
                    for (size_t j = 0; j < number_coefficients; ++j) {
                        for (size_t k = 0; k < number_window_points; ++k) {
                            A[i][j] += functions[j % 2](q[k].x, x[j / 2]) * functions[i % 2](q[k].x, x[i / 2]);
                        }
                    }
                    for (size_t j = 0; j < number_window_points; ++j) {
                        b[i] += q[j].y[0] * functions[i % 2](q[j].x, x[i / 2]);
                    }
                }
                
                minimizer.set_A(A);
                minimizer.set_b(b);
                minimizer.solve(coefficients);

                omega = x;
                is_min = true;

                return std::abs(u.computeObjFunction(q[index - 3].x) - q[index - 3].y[0]) - delta;
            
            case 6:
                omega = x;

                is_min = true;
                gsa.setProblem(u);
                gsa.solve(result);
                min_value = result.value;

                is_min = false;
                gsa.setProblem(u);
                gsa.solve(result);
                max_value = -result.value;

                // std::cout << max_value - min_value - 2 * delta << "\n";

                return max_value - min_value - 2 * delta;

            case 7:
                return -std::abs(u_der(x));

            default: return std::numeric_limits<double>::quiet_NaN();
        }
    }
//    opt::MultiDimensionalSearchArea(N + 1, std::vector<double> { 0.01, 0.01, 0.01, 0.01 },
//                                           std::vector<double> { 2.0,  2.0,  2.0,  2.0 }), 2 * N + 1);
std::vector<double> A{ 0.01, 0.01, 0.01, 0.01 };
std::vector<double> B{ 2.0,  2.0,  2.0,  2.0 };

int main() {
    double accuracy = 0.01, reliability = 3.0, d = 0.0;
    int maxTrials = 10000, maxFevals = 10000;
    int den = 10, key = 3 , incr = 1;

    double total_start_time = omp_get_wtime();
#if defined( CALC )
    OutputFile vars_file("output_data/sample_test_problem_family/vars.txt");
    if (!vars_file.isOpen()) std::cerr << "vars_file opening error\n";

    size_t number_methods = 1;
    MggsaMethod mggsa;

    vars_file.setVariable("number_methods", number_methods, false);
    vars_file.setVariable("A", a, false);
    vars_file.setVariable("B", b, false);
    vars_file.setVariable("delta", delta, false);
    vars_file.setVariable("number_coefficients", number_coefficients, false);
    vars_file.initArray("method_names", number_methods);
    // vars_file.init_array("x_opt", number_methods);

    // output_file trialsFile;
    std::vector<double> optimal_points;
    std::vector<Trial> trials;
    double start_time, end_time, work_time;

    int numberTrials, numberFevals;
    std::vector<double> X, L;

    // trialsFile.open("output_data/sample_test_problem_family/" + method_name + "_trials.txt");
    // if (!trialsFile.is_open()) std::cerr << "trialsFile opening error\n";

    mggsa.setF(problem);
    mggsa.setN(N + 1);
    mggsa.setNumberConstraints(2 * N + 1);
    // mggsa.setNumberConstraints(3);
    mggsa.setAB(A, B);
    mggsa.setEps(accuracy);
    mggsa.setMaxTrials(maxTrials);
    mggsa.setMaxFevals(maxFevals);
    mggsa.setR(reliability);
    mggsa.setD(d);
    mggsa.setDen(den);
    mggsa.setKey(key);
    mggsa.setIncr(incr);

    start_time = omp_get_wtime();
    mggsa.solve(numberTrials, numberFevals, X);
    mggsa.getL(L);
    end_time = omp_get_wtime();
    work_time = end_time - start_time;

    std::cout << "Method name: " << method_name << '\n';
    printResultMggsa("sample_test_problem_family", N + 1, 2 * N + 1, A, B, L, std::vector<double>{0, 0, 0, 0}, -1,
                     maxTrials, maxFevals, accuracy, reliability, d, den, key, incr, numberTrials, numberFevals, L, X,
                     problem(X, 2 * N + 1));
    std::cout << "max |u_der| = " << std::abs(problem(X, 2 * N + 1)) << std::endl; 

    problem(X, 3);
    std::cout << "c = (" << coefficients[0];
    for (int i = 1; i < number_coefficients; ++i) {
        std::cout << "; " << coefficients[i];
    }
    std::cout << ")\n\n";

// 
            // trialsFile.add_point(taskArray[i].X_opt, taskArray[i].f(taskArray[i].X_opt, taskArray[i].numberConstraints));
            // trialsFile.add_point(X, taskArray[i].f(X, taskArray[i].numberConstraints));
// 
            // mggsa.getPoints(points);
            // trialsFile.add_points(points);
// 
            // trialsFile.close();

    vars_file.setValueInArray("method_names", 1, method_name);
    vars_file.initArray("x_opt", X, false);
    vars_file.initArray("c", coefficients, false);
    vars_file.initArray("T", 3);
    vars_file.initArray("Q", 3);
    for (int i = 0; i < 3; ++i) {
        vars_file.setValueInArray("T", i + 1, q[i + 6].x, false);
        vars_file.setValueInArray("Q", i + 1, q[i + 6].y[0], false);
    }

    // sample_test_problem.getOptimalPoints(optimal_points);
    // trialsFile.add_points(optimal_points, sample_test_problem.getOptimalValue());
    // trialsFile.add_point(result.point, result.value);

    // methods[i]->getTrialPoints(trials);
    // trialsFile.add_points(trials);

    // trialsFile.close();

    vars_file.close();

    OutputFile test_points_file("output_data/sample_test_problem_family/test_points.txt");
    if (!test_points_file.isOpen()) std::cerr << "test_points_file opening error\n";

    for (int i = 0; i < number_window_points; ++i) {
        test_points_file.addPoint(q[i].x, q[i].y[0], false);
    }

    test_points_file.close();

    // output_file function_points_file("output_data/sample_test_problem/function_points.txt");
    // if (!function_points_file.is_open()) std::cerr << "function_points_file opening error\n";
// 
    // double a = sample_test_problem.getSearchArea().lowerBound;
    // double b = sample_test_problem.getSearchArea().upBound;
    // for (double i = a; i < b; i += step) {
    //     function_points_file.add_point(i, sample_test_problem.computeObjFunction(i), false);
    // }

    // function_points_file.close();
#endif
    double total_end_time = omp_get_wtime();
    std::cout << "Total time: " << total_end_time - total_start_time << "\n";

#if defined( DRAW )
    Script script("scripts/sample_test_problem_family.gp");
    script.addArgs(std::vector<int>{ display_type, graph_type, method_index });
    script.addArgs(std::vector<double>{ 0, 20,
                                        -1, -1 });
    script.start();
    if (script.isError() == 2) std::cerr << "Error gnuplot\n";
    if (script.isError() == 1) std::cerr << "Error chmod\n";
#endif

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}
