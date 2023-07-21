#ifndef GENERAL_NUMERICAL_TASK_H_
#define GENERAL_NUMERICAL_TASK_H_

#include <base_structures/tasks/GeneralTask.h>

namespace om {
    template <typename OptProblemType>
    struct GeneralNumericalTask : public GeneralTask<OptProblemType> {
        double accuracy, error;
        int maxTrials, maxFevals;

        GeneralNumericalTask(const string &_name, const OptProblemType &_optProblem, double _accuracy,
                             double _error, int _maxTrials, int _maxFevals, bool _use = true)
                            : GeneralTask<OptimizationProblemType>(_name, _optProblem, _use),
                            accuracy(_accuracy), error(_error), maxTrials(_maxTrials), maxFevals(_maxFevals) {};
    };
}

#endif // GENERAL_NUMERICAL_TASK_H_
