#ifndef TASK_H_
#define TASK_H_

#include <string>

using std::string;

template <typename OptimizationTaskType>
struct Task {
    string name;

    OptimizationTaskType optTask;
    int maxFevals;

    bool use;

    Task(string _name, const OptimizationTaskType &_optTask, int _maxFevals, bool _use)
         : name(_name), optTask(_optTask), maxFevals(_maxFevals), use(_use) {};
};

#endif // TASK_H_
