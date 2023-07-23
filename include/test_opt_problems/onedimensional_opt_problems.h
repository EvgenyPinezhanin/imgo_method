#ifndef ONEDIMENSIONAL_OPT_PROBLEMS_H_
#define ONEDIMENSIONAL_OPT_PROBLEMS_H_

#include <vector>
#include <string>

#include <opt_problems/OneDimensionalProblem.h>
#include <base_structures/search_areas/OneDimensionalSearchArea.h>
#include <my_math.h>

using std::vector;
using std::string;
using std::sin;
using std::cos;
using std::pow;

using om::OneDimensionalSearchArea;

const int numberBlocks = 2;
const vector<int> numberFunctions = { 4, 20 };
const vector<string> blockNames = { "sample", "test" };

const vector<OneDimensionalProblem> sampleTasks{
    OneDimensionalProblem([] (double x) { return -4.0 * x + 1.0; },
                          OneDimensionalSearchArea(3.0, 4.0), vector<double>{ 4.0 }, -15.0, 4.0),

    OneDimensionalProblem([] (double x) { return 5.0 * x * x + 3.0 * x - 1.0; },
                          OneDimensionalSearchArea(-2.0, 2.0), vector<double>{ -0.3 }, -1.45, 23.0),

    OneDimensionalProblem([] (double x) { return x * sin(x); },
                          OneDimensionalSearchArea(0.0, 20.0), vector<double>{ 17.336 }, -17.307607, 18.955),

    OneDimensionalProblem([] (double x) { return x != 0 ? x * sin(1 / x) : 0.0; },
                          OneDimensionalSearchArea(-0.4, -0.05), vector<double>{ -0.2225 }, -0.21723353, 6.0 * M_PI)
};

const vector<OneDimensionalProblem> testTasks{
    OneDimensionalProblem(
        [] (double x) {
            return -(-1.0 / 6.0 * pow(x, 6) + 52.0 / 25.0 * pow(x, 5) - 39.0 / 80.0 * pow(x, 4) - 
                     71.0 / 10.0 * pow(x, 3) + 79.0 / 20.0 * x * x + x - 1.0 / 10.0);
        }, OneDimensionalSearchArea(-1.5, 11.0), vector<double>{ 10.0 }, -29763.233, 13870.0),

    OneDimensionalProblem([] (double x) { return -(-sin(x) - sin(10.0 / 3.0 * x)); },
                          OneDimensionalSearchArea(2.7, 7.5), vector<double>{ 5.145735 }, -1.8995993, 4.29),

    OneDimensionalProblem(
        [] (double x) {
            double sum = 0.0;
            for (int i = 1; i <= 5; i++) {
                sum += i * sin((i + 1.0) * x  + i);
            }
            return -sum;
        }, OneDimensionalSearchArea(-10.0, 10.0), vector<double>{ -6.7745761, -0.49139, 5.791785 }, -12.031249, 67.0),

    OneDimensionalProblem([] (double x) { return -(16.0 * x * x - 24.0 * x + 5.0) * exp(-x); },
                          OneDimensionalSearchArea(1.9, 3.9), vector<double>{ 2.868 }, -3.8504507, 3.0),

    OneDimensionalProblem([] (double x) { return -(-3.0 * x + 1.4) * sin(18.0 * x); },
                          OneDimensionalSearchArea(0.0, 1.2), vector<double>{ 0.96609 }, -1.4890725, 36.0),

    OneDimensionalProblem([] (double x) { return -((x + sin(x)) * exp(-x * x)); },
                          OneDimensionalSearchArea(-10.0, 10.0), vector<double>{ 0.67956 }, -0.8242394, 2.5),

    OneDimensionalProblem([] (double x) { return -(-sin(x) - sin(10.0 / 3.0 * x) - log(x) + 0.84 * x - 3.0); },
                          OneDimensionalSearchArea(2.7, 7.5), vector<double>{ 5.19978 }, -1.6013075, 6.0),

    OneDimensionalProblem(
        [] (double x) {
            double sum = 0.0;
            for (int i = 1; i <= 5; i++) {
                sum += i * cos((i + 1.0) * x  + i);
            }
            return -sum;
        }, OneDimensionalSearchArea(-10.0, 10.0), vector<double>{ -7.0835, -0.8003, 5.48286}, -14.508008, 67.0),

    OneDimensionalProblem([] (double x) { return -(-sin(x) - sin(2.0 / 3.0 * x)); },
                          OneDimensionalSearchArea(3.1, 20.4), vector<double>{ 17.039 }, -1.9059611, 1.7),

    OneDimensionalProblem([] (double x) { return -(x * sin(x)); },
                          OneDimensionalSearchArea(0.0, 10.0), vector<double>{ 7.9787 }, -7.9167274, 11.0),

    OneDimensionalProblem([] (double x) { return -(-2.0 * cos(x) - cos(2.0 * x)); },
                          OneDimensionalSearchArea(-1.57, 6.28), vector<double>{ 2.094, 4.189 }, -1.4999998, 3.0),

    OneDimensionalProblem([] (double x) { return -(-pow(sin(x), 3) - pow(cos(x), 3)); },
                          OneDimensionalSearchArea(0.0, 6.28), vector<double>{ 3.142, 4.712 }, -0.99999975, 2.2),

    OneDimensionalProblem(
        [] (double x) {
            return x * x - 1 < 0 ? -(pow(x, 2.0 / 3.0) + pow(-(x * x - 1), 1.0 / 3.0)) :
                           -(pow(x, 2.0 / 3.0) - pow(x * x - 1, 1.0 / 3.0));
        }, OneDimensionalSearchArea(0.001, 0.99), vector<double>{ 0.7071 }, -1.5874011, 8.5),

    OneDimensionalProblem([] (double x) { return -(exp(-x) * sin(2 * M_PI * x)); },
                          OneDimensionalSearchArea(0.0, 4.0), vector<double>{ 0.224885 }, -0.78868539, 6.5),

    OneDimensionalProblem([] (double x) { return -((-x * x + 5.0 * x - 6.0) / (x * x + 1)); },
                          OneDimensionalSearchArea(-5.0, 5.0), vector<double>{ 2.4142 }, -0.035533906, 6.5),

    OneDimensionalProblem([] (double x) { return -(-2.0 * (x - 3) * (x - 3) - exp(- x * x / 2)); },
                          OneDimensionalSearchArea(-3.0, 3.0), vector<double>{ 3.0 }, 0.011108997, 85.0),

    OneDimensionalProblem([] (double x) { return -(-pow(x, 6) + 15.0 * pow(x, 4) - 27.0 * x * x - 250.0); },
                          OneDimensionalSearchArea(-4.0, 4.0), vector<double>{ -3.0, 3.0 }, 7.0, 2520.0),

    OneDimensionalProblem([] (double x) { return x <= 3.0 ? (x - 2.0) * (x - 2.0) : -(-2.0 * log(x - 2.0) - 1.0); },
                          OneDimensionalSearchArea(0.0, 6.0), vector<double>{ 2.0 }, 0.0, 4.0),

    OneDimensionalProblem([] (double x) { return -(x - sin(3.0 * x) + 1.0); },
                          OneDimensionalSearchArea(0.0, 6.5), vector<double>{ 5.87287 }, -7.8156745, 4.0),

    OneDimensionalProblem([] (double x) { return -(x - sin(x)) * exp(- x * x); },
                          OneDimensionalSearchArea(-10.0, 10.0), vector<double>{ 1.195137 }, -0.063490529, 1.3)
};

#endif // ONEDIMENSIONAL_OPT_PROBLEMS_H_
