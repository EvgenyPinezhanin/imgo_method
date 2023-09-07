#ifndef SCANNING_TASK_H_
#define SCANNING_TASK_H_

#include <tasks/Task.h>
#include <opt_problems/OneDimensionalProblem.h>
#include <opt_methods/ScanningMethod.h>

using ScanningTask = Task<OneDimensionalProblem, ScanningParameters>;

#endif // SCANNING_TASK_H_
