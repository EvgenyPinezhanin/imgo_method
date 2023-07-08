#ifndef GENERAL_TASK_H_
#define GENERAL_TASK_H_

#include <string>

using std::string;

template <typename OptimizationTaskType>
struct GeneralTask {
    string name;

    OptimizationTaskType optTask;
    int maxFevals;

    bool use;

    GeneralTask(string _name, const OptimizationTaskType &_optTask, int _maxFevals, bool _use)
                : name(_name), optTask(_optTask), maxFevals(_maxFevals), use(_use) {};
};

#endif // GENERAL_TASK_H_
