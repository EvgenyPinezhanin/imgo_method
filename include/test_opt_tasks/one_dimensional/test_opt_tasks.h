#ifndef TEST_OPT_TASKS_H_
#define TEST_OPT_TASKS_H_

#include <vector>

#include <opt_tasks/OneDimensionalTask.h>
#include <my_math.h>

using std::vector;
using std::sin;
using std::cos;
using std::pow;

const vector<OneDimensionalTask> testTasks{
    OneDimensionalTask(
        [] (double x) {
            return -(-1.0 / 6.0 * pow(x, 6) + 52.0 / 25.0 * pow(x, 5) - 39.0 / 80.0 * pow(x, 4) - 
                     71.0 / 10.0 * pow(x, 3) + 79.0 / 20.0 * x * x + x - 1.0 / 10.0);
        }, OneDimensionalSearchArea(-1.5, 11.0), 10.0, 13870.0),

    OneDimensionalTask([] (double x) { return -(-sin(x) - sin(10.0 / 3.0 * x)); },
                       OneDimensionalSearchArea(2.7, 7.5), 5.145735, 4.29),

    OneDimensionalTask(
        [] (double x) {
            double sum = 0.0;
            for (int i = 1; i <= 5; i++) {
                sum += i * sin((i + 1.0) * x  + i);
            }
            return -sum;
        }, OneDimensionalSearchArea(-10.0, 10.0), -0.49139, 67.0),

    OneDimensionalTask([] (double x) { return -(16.0 * x * x - 24.0 * x + 5.0) * exp(-x); },
                       OneDimensionalSearchArea(1.9, 3.9), 2.868, 3.0),

    OneDimensionalTask([] (double x) { return -(-3.0 * x + 1.4) * sin(18.0 * x); },
                       OneDimensionalSearchArea(0.0, 1.2), 0.96609, 36.0),

    OneDimensionalTask([] (double x) { return -((x + sin(x)) * exp(-x * x)); },
                       OneDimensionalSearchArea(-10.0, 10.0), 0.67956, 2.5),

    OneDimensionalTask([] (double x) { return -(-sin(x) - sin(10.0 / 3.0 * x) - log(x) + 0.84 * x - 3.0); },
                       OneDimensionalSearchArea(2.7, 7.5), 5.19978, 6.0),

    OneDimensionalTask(
        [] (double x) {
            double sum = 0.0;
            for (int i = 1; i <= 5; i++) {
                sum += i * cos((i + 1.0) * x  + i);
            }
            return -sum;
        }, OneDimensionalSearchArea(-10.0, 10.0), 5.48286, 67.0),

    OneDimensionalTask([] (double x) { return -(-sin(x) - sin(2.0 / 3.0 * x)); },
                       OneDimensionalSearchArea(3.1, 20.4), 17.039, 1.7),

    OneDimensionalTask([] (double x) { return -(x * sin(x)); },
                       OneDimensionalSearchArea(0.0, 10.0), 7.9787, 11.0),

    OneDimensionalTask([] (double x) { return -(-2.0 * cos(x) - cos(2.0 * x)); },
                       OneDimensionalSearchArea(-1.57, 6.28), 2.094, 3.0),

    OneDimensionalTask([] (double x) { return -(-pow(sin(x), 3) - pow(cos(x), 3)); },
                       OneDimensionalSearchArea(0.0, 6.28), 4.712, 2.2),

    OneDimensionalTask(
        [] (double x) {
            return x * x - 1 < 0 ? -(pow(x, 2.0 / 3.0) + pow(-(x * x - 1), 1.0 / 3.0)) :
                           -(pow(x, 2.0 / 3.0) - pow(x * x - 1, 1.0 / 3.0));
        }, OneDimensionalSearchArea(0.001, 0.99), 0.7071, 8.5),

    OneDimensionalTask([] (double x) { return -(exp(-x) * sin(2 * M_PI * x)); },
                       OneDimensionalSearchArea(0.0, 4.0), 0.224885, 6.5),

    OneDimensionalTask([] (double x) { return -((-x * x + 5.0 * x - 6.0) / (x * x + 1)); },
                       OneDimensionalSearchArea(-5.0, 5.0), 2.4142, 6.5),

    OneDimensionalTask([] (double x) { return -(-2.0 * (x - 3) * (x - 3) - exp(- x * x / 2)); },
                       OneDimensionalSearchArea(-3.0, 3.0), 3.0, 85.0),

    OneDimensionalTask([] (double x) { return -(-pow(x, 6) + 15.0 * pow(x, 4) - 27.0 * x * x - 250.0); },
                       OneDimensionalSearchArea(-4.0, 4.0), -3.0, 2520.0),

    OneDimensionalTask([] (double x) { return x <= 3.0 ? (x - 2.0) * (x - 2.0) : -(-2.0 * log(x - 2.0) - 1.0); },
                       OneDimensionalSearchArea(0.0, 6.0), 2.0, 4.0),

    OneDimensionalTask([] (double x) { return -(x - sin(3.0 * x) + 1.0); },
                       OneDimensionalSearchArea(0.0, 6.5), 5.87287, 4.0),

    OneDimensionalTask([] (double x) { return -(x - sin(x)) * exp(- x * x); },
                       OneDimensionalSearchArea(-10.0, 10.0), 1.195137, 1.3)
};

#endif // TEST_OPT_TASKS_H_
