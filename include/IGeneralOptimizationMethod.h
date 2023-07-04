#ifndef I_GENERAL_OPTIMIZATION_METHOD_H_
#define I_GENERAL_OPTIMIZATION_METHOD_H_

#include <vector>

#include <ResultMethod.h>

using std::vector;

template <typename TaskOptimizationMethodType>
class IGeneralOptimizationMethod {
protected:
    TaskOptimizationMethodType task;

    int maxFevals, numberFevals;

public:
    IGeneralOptimizationMethod(const TaskOptimizationMethodType &_task, int _maxFevals)
                               : task(_task), maxFevals(_maxFevals), numberFevals(0) {};

    void setTask(const TaskOptimizationMethodType &_task) { task = _task; };
    void getTask(TaskOptimizationMethodType &_task) const { _task = task; };

    void setMaxFevals(int _maxFevals) { maxFevals = _maxFevals; };
    int getMaxFevals() const { return maxFevals; };

    virtual void solve(ResultMethod *result) = 0;
    virtual bool solveTest(vector<double> XOpt, ResultMethod *result) = 0;
};

#endif // I_GENERAL_OPTIMIZATION_METHOD_H_
