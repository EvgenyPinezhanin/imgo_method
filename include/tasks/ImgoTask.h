#ifndef IMGO_TASK_H_
#define IMGO_TASK_H_

#include <vector>
#include <string>

#include <tasks/Task.h>
// #include <opt_problems/OneDimensionalConstrainedProblem.h>

using std::vector;
using std::string;

/* template <typename OptimizationProblemType>
struct StronginTask : public GeneralNumericalTask<OptimizationProblemType> {
    vector<double> reliability, constantsEstimation, d;

    StronginTask(string _name, const OptimizationProblemType &_optTask, vector<double> _r, vector<double> constantEstimation,
                 vector<double> _d, double _accuracy, int _maxTrials, int _maxFevals, bool _use)
                 : GeneralNumericalTask<OptimizationProblemType>(_name, _optTask, _accuracy, _maxTrials, _maxFevals, use),
                 r(_r), constantEstimation(_constantEstimation), d(_d) {};
}; */

#endif // CONSTRAINED_STRONGIN_TASK_H_


/* struct TaskImgo : public Task {
    double (*f)(double, int);
    int numberConstraints;
    double r, d;

    TaskImgo(double (*_f)(double, int), string _name, int _numberConstraints, double _a, double _b, double _xOpt,
             vector<double> _L, double _eps, int _maxTrials, int _maxFevals, double _r, double _d, bool _used = true)
             : Task(_name, 1, vector<double>{_a}, vector<double>{_b}, vector<double>{_xOpt}, _L, _eps, _maxTrials,
             _maxFevals, _used), f(_f), numberConstraints(_numberConstraints), r(_r), d(_d) {};
}; */