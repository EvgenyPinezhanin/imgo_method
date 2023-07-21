#ifndef GENERAL_TASK_H_
#define GENERAL_TASK_H_

#include <string>

using std::string;

template <typename OptProblemType>
struct GeneralTask {
    string name;

    OptProblemType optProblem;

    bool use;

    GeneralTask(string _name, const OptProblemType &_optProblem, bool _use)
                : name(_name), optProblem(_optProblem), use(_use) {};
};

#endif // GENERAL_TASK_H_
