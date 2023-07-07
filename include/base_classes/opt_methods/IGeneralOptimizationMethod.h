#ifndef I_GENERAL_OPTIMIZATION_METHOD_H_
#define I_GENERAL_OPTIMIZATION_METHOD_H_

template <typename OptimizationTaskType, typename ResultMethodType>
class IGeneralOptimizationMethod {
protected:
    OptimizationTaskType task;

    int maxFevals, numberFevals;

    virtual void setDataInResultMethod(ResultMethodType &result) const = 0;

public:
    IGeneralOptimizationMethod(const OptimizationTaskType &_task, int _maxFevals)
                               : task(_task), maxFevals(_maxFevals), numberFevals(0) {};

    void setTask(const OptimizationTaskType &_task) { task = _task; };
    OptimizationTaskType getTask() const { return task; };

    void setMaxFevals(int _maxFevals) { maxFevals = _maxFevals; };
    int getMaxFevals() const { return maxFevals; };

    virtual void solve(ResultMethodType &result) = 0;
    virtual bool solveTest(ResultMethodType &result) = 0;
};

#endif // I_GENERAL_OPTIMIZATION_METHOD_H_
