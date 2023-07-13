#ifndef GENERAL_NUMERICAL_TASK_H_
#define GENERAL_NUMERICAL_TASK_H_

#include <base_structures/tasks/GeneralTask.h>

template <typename OptimizationProblemType>
struct GeneralNumericalTask : public GeneralTask<OptimizationProblemType> {
    double accuracy, error;
    int maxTrials;

    GeneralNumericalTask(string _name, const OptimizationProblemType &_optTask, double _accuracy,
                         double _error, int _maxTrials, int _maxFevals, bool _use = true)
                         : GeneralTask<OptimizationProblemType>(_name, _optTask, _maxFevals, _use),
                         accuracy(_accuracy), error(_error), maxTrials(_maxTrials) {};
};

#endif // GENERAL_NUMERICAL_TASK_H_
