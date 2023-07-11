#ifndef SAMPLE_OPT_TASKS_H_
#define SAMPLE_OPT_TASKS_H_

#include <vector>

#include <opt_tasks/OneDimensionalTask.h>
#include <myMath.h>

using std::vector;
using std::sin;

const vector<OneDimensionalTask> sampleTasks{
    OneDimensionalTask([] (double x) { return -4.0 * x + 1.0; }, OneDimensionalSearchArea(3.0, 4.0), 4.0, 4.0),
    OneDimensionalTask([] (double x) { return 5.0 * x * x + 3.0 * x - 1.0; },
                       OneDimensionalSearchArea(-2.0, 2.0), -0.3, 23.0),
    OneDimensionalTask([] (double x) { return x * sin(x); }, OneDimensionalSearchArea(0.0, 20.0), 17.336, 18.955),
    OneDimensionalTask([] (double x) { return x != 0 ? x * sin(1 / x) : 0.0; },
                       OneDimensionalSearchArea(-0.4, -0.05), -0.2225, 6.0 * M_PI)
};

#endif // SAMPLE_OPT_TASKS_H_
