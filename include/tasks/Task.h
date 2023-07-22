#ifndef TASK_H_
#define TASK_H_

#include <base_structures/tasks/GeneralNumericalTask.h>

template <typename OptProblemType>
struct Task : public om::GeneralNumericalTask<OptProblemType> {
    int blockNumber, functionNumber;

    Task(string _name, const OptProblemType &_optProblem, int _blockNumber, int _functionNumber,
         double _accuracy, double _error, int _maxTrials, int _maxFevals, bool _use = true)
        : om::GeneralNumericalTask<OptProblemType>(_name, _optProblem, _accuracy, _error,
        _maxTrials, _maxFevals, _use), blockNumber(_blockNumber), functionNumber(_functionNumber) {};
};

#endif // TASK_H_
