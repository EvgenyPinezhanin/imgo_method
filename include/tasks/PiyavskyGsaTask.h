#ifndef PIYAVSKY_GSA_TASK_H_
#define PIYAVSKY_GSA_TASK_H_

#include <string>

#include <tasks/Task.h>
#include <opt_problems/OneDimensionalProblem.h>

using std::string;

struct PiyavskyGsaTask : public Task<OneDimensionalProblem> {
    double reliability;

    PiyavskyGsaTask(string _name, const OneDimensionalProblem &_optProblem, int _blockNumber,
                    int _functionNumber, double _reliability, double _accuracy, double _error,
                    int _maxTrials, int _maxFevals, bool _use = true)
                   : Task<OneDimensionalProblem>(_name, _optProblem, _blockNumber, _functionNumber,
                   _accuracy, _error, _maxTrials, _maxFevals, _use), reliability(_reliability) {};
};

#endif // PIYAVSKY_GSA_TASK_H_
