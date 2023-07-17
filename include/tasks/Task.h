#ifndef TASK_H_
#define TASK_H_

#include <base_structures/tasks/GeneralNumericalTask.h>

template <typename OptimizationProblemType>
struct Task : public GeneralNumericalTask<OptimizationProblemType> {
    int blockNumber, functionNumber;

    Task(string _name, const OptimizationProblemType &_optProblem, int _blockNumber,
         int _functionNumber, double _accuracy, double _error, int _maxTrials, int _maxFevals, bool _use = true)
        : GeneralNumericalTask<OptimizationProblemType>(_name, _optProblem, _accuracy, _error, _maxTrials,
        _maxFevals, _use), blockNumber(_blockNumber), functionNumber(_functionNumber) {};
};

#endif // TASK_H_
