#ifndef TASK_H_
#define TASK_H_

#include <base_structures/tasks/GeneralTask.h>

template <typename OptProblemType, typename ParametersMethodType>
struct Task : public opt::GeneralTask<OptProblemType, ParametersMethodType> {
    int blockNumber, functionNumber;

    Task(string _name, int _blockNumber, int _functionNumber, const OptProblemType &_optProblem,
         const ParametersMethodType &_parameters, bool _use = true):
        opt::GeneralTask<OptProblemType, ParametersMethodType>(_name, _optProblem, _parameters, _use),
        blockNumber(_blockNumber),
        functionNumber(_functionNumber)
    {};
};

#endif // TASK_H_
