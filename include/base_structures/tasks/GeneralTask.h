#ifndef GENERAL_TASK_H_
#define GENERAL_TASK_H_

#include <string>

using std::string;

namespace opt {
    template <typename OptProblemType, typename ParametersMethodType>
    struct GeneralTask {
        string name;

        OptProblemType optProblem;
        ParametersMethodType parameters;

        bool use;

        GeneralTask(const string &_name, const OptProblemType &_optProblem,
                    const ParametersMethodType &_parameters, bool _use):
            name(_name),
            optProblem(_optProblem),
            parameters(_parameters),
            use(_use)
        {};
    };
}

#endif // GENERAL_TASK_H_
