#ifndef GENERAL_TASK_H_
#define GENERAL_TASK_H_

#include <string>

using std::string;

template <typename OptimizationProblemType>
struct GeneralTask {
    string name;

    OptimizationProblemType optProblem;
    int maxFevals;

    bool use;

    GeneralTask(string _name, const OptimizationProblemType &_optTask, int _maxFevals, bool _use)
                : name(_name), optProblem(_optTask), maxFevals(_maxFevals), use(_use) {};
};

#endif // GENERAL_TASK_H_
