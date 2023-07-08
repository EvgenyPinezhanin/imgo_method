#ifndef GENERAL_NUMERICAL_TASK_H_
#define GENERAL_NUMERICAL_TASK_H_

#include <base_structures/tasks/GeneralTask.h>

template <typename OptimizationTaskType>
struct GeneralNumericalTask : public GeneralTask<OptimizationTaskType> {
    double accuracy;
    int maxTrials;

    GeneralNumericalTask(string _name, const OptimizationTaskType &_optTask, double _accuracy,
                         int _maxTrials, int _maxFevals, bool _use)
                         : GeneralTask<OptimizationTaskType>(_name, _optTask, _maxFevals, _use),
                         accuracy(_accuracy), maxTrials(_maxTrials) {};
};

#endif // GENERAL_NUMERICAL_TASK_H_
